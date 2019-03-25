module mps_mod
use gen_util_mod
use types_mod
use umat_mod
use output_mod
use mps_util_mod
use sim_writing_util_mod
use sim_getincr_mod
use sim_error_mod

implicit none
    
    private
    public      ::  mp_simulate
    
    contains 
    
    
subroutine mp_simulate(error, props, f_data, simnr, evec)
implicit none
double precision, intent(out)           :: error
double precision, intent(in)            :: props(:)
type(fdata_typ), intent(inout)          :: f_data
integer, intent(in)                     :: simnr
double precision, optional, allocatable :: evec(:)

! Variables from global settings
logical                         :: save_results     ! Should results be saved to file?
procedure(umat_template),pointer:: umat_address     ! Addresss to umat subroutine
character(len=strl)             :: cmname           ! Material model name
logical                         :: nlgeom           ! Geometrically nonlinear analysis?
integer                         :: psize            ! Problem size (6 or 9), derived from nlgeom

! Variables derived from simulation settings
logical                         :: result_onlymain  ! Only include results from main increments, or from all
logical                         :: result_inclexp   ! Include experiment data in result file?
double precision, allocatable   :: result_steps(:)  ! Which steps to write results for (if = 0 then all, if =-1 then none)
double precision, allocatable   :: error_steps(:)   ! Which steps to save error for (if =0 then all)

! Loading
integer, allocatable            :: exp_info(:)      ! Columns in exp_data. If 0 then column does not exist (and will be set to zero). 
                                                    ! [Step, time, forc, astr, torq, rota, p_i, cstr_i, p_o, cstr_o, temp] (cstr=circumferential strain)
! Iteration settings
type(iter_typ)                  :: iter             ! Iteration settings

! Objective function scaling values
integer                         :: err_norm_met     ! Method for scaling errors (i.e. no scaling, max-min value in cycle, max-min value in simulation)

! Internal variables
 ! Error vector
logical                         :: eout             ! Logical if evec should be given
integer                         :: e_cnt            ! Number of times the error is evaluated
integer                         :: evec_cnt         ! Counter for current position in evec
integer                         :: evec_est_length  ! Estimated length for evec
 ! Time stepping
double precision                :: stp_time_incr(4) ! Current value for time_incr (dtmain, dtmin, dt0, dtmax)
integer, allocatable            :: stp_ctrl(:)      ! Current value for ctrl
double precision                :: stp_time_tot     ! Total time for current step
integer                         :: kstep, kinc, convincr
double precision                :: step             ! Double with step number (can be decimal)
 ! Error scaling
double precision, allocatable   :: stp_err_scale(:) ! Current value for err_scale
logical                         :: cyclic_error
double precision, allocatable   :: cyc_fsim(:,:), cyc_fexp(:,:), cyc_time(:)
integer                         :: cyc_ptcnt
double precision                :: last_cyc_start_step
 ! Experimental data
double precision, allocatable   :: expdata(:,:)     ! Experiment data
double precision, allocatable   :: stp_time_vec(:)  ! Temporary storage of time points in current step

 ! Other internals
! time stepping
double precision                :: time(2), pnewdt, dt, dt_old
double precision                :: dtmain, dtmin, dtmax, dtmain_true
double precision, allocatable   :: tmain(:)                             !main time increments (used for error calculation)
integer, allocatable            :: stprows(:)           ! Which row each step starts on
! 
double precision, allocatable   :: disp(:), disp_old(:), disp_exp(:)    !strain or dfgrd-I
double precision, allocatable   :: v(:), v_old(:), du(:)                !velocity of disp for use for qualified guessing of values at next increment
double precision, allocatable   :: load(:), load_old(:), load_exp(:)    !stress measure
double precision, allocatable   :: stat(:), stat_old(:)                 !state variables
double precision                :: temp, temp_old                       !temperature
double precision, allocatable   :: additional_output(:)                 !user chosen additional output variables
double precision                :: adjustment(2,4)      ! Adjustment variable for smooth following of experimental data if control mode is changed
integer, allocatable            :: free_dofs(:)         ! Free degrees of freedom in the problem
!
double precision, allocatable   :: evec_tmp(:)          !Temp variable used to resize the evec output
!
logical                         :: laststep, lconv, res_in_step, err_in_step, resize_evec
integer                         :: nstep, niter, kmain
integer                         :: exp_row    ! Row in experiment data
! program timing
double precision                :: starttime, stoptime
integer                         :: clock_count, clock_rate


! == Simulation timing
call system_clock ( clock_count, clock_rate)
starttime = real(clock_count)/real(clock_rate)

! == Material ==
umat_address => f_data%glob%umat_address
cmname       = f_data%glob%cmname
nlgeom       = f_data%glob%nlgeom

! == Problem size ==
! Depending on analysis type, stress (load) and disp (strain/defgrad) will have different sizes
if (nlgeom) then
    psize = 9
else
    psize = 6
endif

allocate(disp(psize), disp_old(psize), disp_exp(psize), v(psize), v_old(psize), du(psize))
allocate(load(psize), load_old(psize), load_exp(psize))
allocate(stat(f_data%glob%nstatv), stat_old(f_data%glob%nstatv))
allocate(additional_output(1)); additional_output = 0.d0
allocate(stp_ctrl(psize), stp_err_scale(psize))

! == Use shorter names for some variables == 
! slight performance loss, but increases readability
! Experiment
allocate(exp_info, source=f_data%sim(simnr)%exp%exp_info)
allocate(expdata, source=f_data%sim(simnr)%sim_setup%expdata_array)
allocate(stprows, source=f_data%sim(simnr)%sim_setup%stprows)

nstep = size(stprows)-1
! Iteration settings
iter = f_data%sim(simnr)%iter

! Result output
allocate(result_steps, source=f_data%sim(simnr)%outp%result_steps)
save_results    = (f_data%glob%resnr>0).and.(any(result_steps>-1.d-10))
result_onlymain = f_data%sim(simnr)%outp%result_onlymain
result_inclexp  = f_data%sim(simnr)%outp%result_inclexp

! == Get error vector ready if needed ==
resize_evec = .false.
eout = present(evec)
if (eout) then
    if (f_data%sim(simnr)%n_errors/=-1) then ! If a previous simulation already has found the correct number of errors
        evec_est_length = max(f_data%sim(simnr)%n_errors,1)
    else
        evec_est_length = 4*(int(expdata(size(expdata,1), exp_info(2))/minval(iter%time_incr(:,2)))+1) !Worst case size
        resize_evec = .true.
    endif
    allocate(evec(evec_est_length))
endif
evec_cnt = 1    ! Initiate evec_cnt (also use to get n_errors)
e_cnt = 0       ! Initiate e_cnt to zero (number of error evaluations)
error = 0.d0    ! Set error to zero

!Define initial values
time = expdata(1,exp_info(2))               ! Start time 
kstep = 0; kinc  = 0; niter = 0;            ! 
temp_old = f_data%sim(simnr)%init%temp_init     ! If not given in expdata, this value will be used throughout

exp_row = 2 ! Start at row 2 to allow interpolation (if time=t_start then interpolated value will be to row 1 anyway)
call get_incr(temp_old, load_exp, disp_exp, temp, time(2), expdata, exp_info, exp_row)

! Should add option to specify initial conditions to have continued analysis
load_old = 0.d0
disp_old = 0.d0
if (allocated(f_data%sim(simnr)%init%statev_init)) then
    stat_old = f_data%sim(simnr)%init%statev_init(1,:)
else
    stat_old = 0.d0
endif

disp = disp_old
load = load_old
stat = stat_old

! Error settings
call error_settings(f_data%sim(simnr)%err, psize, err_norm_met, error_steps, & 
    cyclic_error, cyc_fsim, cyc_fexp, cyc_time, cyc_ptcnt, last_cyc_start_step)

! Output file
if (save_results) then
    call setup_result_output(trim(f_data%glob%outname)//'_sim'//int2str(simnr)//'_'//int2str(f_data%glob%resnr)//'.txt', result_inclexp, f_data%sim(simnr)%stype, nlgeom, f_data%sim(simnr)%outp%dbl_format)
endif

STEP_LOOP: do kstep = 1,nstep
    step = f_data%sim(simnr)%sim_setup%steps(kstep)
    ! Settings for current step
    call get_step_data_dbl(stp_time_incr, iter%time_incr, step)
    call get_step_data_int(stp_ctrl, f_data%sim(simnr)%exp%ctrl, step)
    call get_step_data_dbl(stp_err_scale, f_data%sim(simnr)%err%err_scale, step)
    
    ! Initialize current step (reset variables and read in new)
    !load_step0 = load;  disp_step0 = disp;  temp_step0 = temp
    convincr    = 0;
    dtmain      = stp_time_incr(1); dtmin = stp_time_incr(2)
    dt          = stp_time_incr(3); dtmax = stp_time_incr(4)
    v           = 0.d0
    time(1)     = 0.d0; kinc       = 0
    
    
    ! Determine main time steps: tmain(nmain)
    allocate(stp_time_vec(stprows(kstep+1)+1-stprows(kstep)))
    stp_time_vec = expdata(stprows(kstep):stprows(kstep+1), exp_info(2))
    call get_tmain(tmain, stp_time_tot, stp_time_vec, dtmain)    
    deallocate(stp_time_vec)
    !call get_tmain(tmain, stp_time_tot, pack(expdata(stprows(kstep):stprows(kstep+1), exp_info(2)), .true.), dtmain)    
    
    ! Check if results should be written for current step
    res_in_step = (save_results).and.((all(dbl_comp_array(result_steps, 0.d0))).or.(any(dbl_comp_array(result_steps, step))))
    err_in_step = (all(dbl_comp_array(error_steps, 0.d0))).or.(any(dbl_comp_array(error_steps, step)))
    
    
    
    ! Find free degrees of freedom for current step
    call mps_gen_free_dofs(free_dofs, stp_ctrl)
    
    ! Find adjustments if any to the experiment data 
    ! E.g. to offset to match initial calculated value in step, 
    ! or offset but scale to get same end value as the experiment data
    adjustment = get_incr_adj(load, disp, expdata(stprows(kstep):stprows(kstep+1), :), exp_info, stp_ctrl)
    
    
    
    TMAIN_LOOP: do kmain=1,size(tmain)  ! Loop over all main increments
        laststep = .false.
        
        ! Calculate error after at the beginning of each main increment (the last is part of next step)
        if (err_in_step) then
            if (cyclic_error) then
                ! Add to cyc_fsim, cyc_fexp and cyc_time
                call cyclic_error_update(load_exp, load, disp_exp, disp, stp_ctrl, time, &
                                        cyc_ptcnt, cyc_fsim, cyc_fexp, cyc_time)
            else
                call error_update(load_exp, load, disp_exp, disp, stp_ctrl, stp_err_scale, &
                                    f_data%sim(simnr)%sim_setup%load_error_scale(:,kstep), &
                                    f_data%sim(simnr)%sim_setup%disp_error_scale(:,kstep), &
                                    error, e_cnt, evec_cnt, evec)
            endif
        endif
        
        ! Write out results for the main increment if it should for current step
        if (res_in_step) then
            call write_result(step, kinc, niter, time(2), temp, load, load_exp, disp, disp_exp, additional_output, result_inclexp)
        endif
        
        dtmain_true = tmain(kmain)-time(1)  ! True main time increment
        do while (.not.laststep)
            ! Set starting variables for current increment
            disp = disp_old
            load = load_old
            stat = stat_old
            pnewdt    = 1.d0
        
            ! Update time for this increment
            call timeupdate(kinc, time, dt, laststep, tmain(kmain), dtmin)
        
            ! Find load_exp, disp_exp and temp
            call get_incr(temp_old, load_exp, disp_exp, temp, time(2), expdata, exp_info, exp_row)
            
            ! Determine load and disp (adjusted if requested by ctrl values)
            call adjust_incr(load_exp, disp_exp, load, disp, stp_ctrl, adjustment, time(2))

            ! Make a guess for the displacement if not all displacements are prescribed
            if (free_dofs(1)/=-1) then
                call update_guess(du, v, dt, v_old, dt_old, kinc, free_dofs)
                disp(free_dofs) = disp(free_dofs) + du(free_dofs)
            endif
            
            ! Solve the increment
            call mps_solve_incr(load_old, disp_old, temp_old, stat_old, load, disp, temp, stat, time, dt, free_dofs, &
                                iter, niter, lconv, pnewdt, kinc, kstep, props, cmname, umat_address, nlgeom)
            
            ! Check results and if OK go to next step
            if (lconv) then
                v_old = v; v = (disp-disp_old)/dt; dt_old = dt                      ! Update values for better guessing
                load_old = load; disp_old = disp; stat_old = stat; temp_old = temp  ! Update actual result values
                
                ! Write out results if that is requested for current simulation/step/increment, and it is not the laststep (then it will be written later)
                if ((res_in_step).and.(.not.result_onlymain).and.(.not.laststep)) then
                    call write_result(step, kinc, niter, time(2), temp, load, load_exp, disp, disp_exp, additional_output, result_inclexp)
                endif
                
                ! Update time stepping
                convincr = convincr + 1
                if (convincr>iter%nconv_incr) then !Increase dt if enough converged iterations has precided this increment
                    dt = min(dt*iter%dt_incr, dtmax)
                endif
            else
                if (dt<=dtmin*(1.d0+1.d-5)) then ! Note that this tolerance must be larger than dt_last_tol in timeupdate!
                    call write_output('Time increment smaller than dtmin requested, exiting simulation at increment '//int2str(kinc), 'warning', 'sim:mps')
                    call write_output('Simulation nr '//int2str(simnr)//' in step '//dbl2str(step), 'warning', 'sim:mps', loc=.false.)
                    call write_output('Setting error to a huge number', 'warning', 'sim:mps', loc=.false.)
                    error = huge(1.d0)
                    exit STEP_LOOP
                endif
                convincr = 0            !Set convincr: dt will not be increased before enough converged increments
                kinc    = kinc - 1      !Reset increment counter
                time    = time - dt     !Reset times
                dt      = max(dt*pnewdt, dtmin)  !Update time step length
                laststep= .false.       !To allow for a smaller last step if required
            endif
        
        enddo
        
    enddo TMAIN_LOOP
    
    if (kstep==nstep) then
        ! Write out results for last main step in last step, if it should for this step
        if (res_in_step) then
            call write_result(step, kinc, niter, time(2), temp, load, load_exp, disp, disp_exp, additional_output, result_inclexp)
        endif
        
        if (err_in_step) then
            if (cyclic_error) then
                ! Add to cyc_fsim, cyc_fexp and cyc_time
                call cyclic_error_update(load_exp, load, disp_exp, disp, stp_ctrl, time, &
                                        cyc_ptcnt, cyc_fsim, cyc_fexp, cyc_time)
            else
                call error_update(load_exp, load, disp_exp, disp, stp_ctrl, stp_err_scale, &
                                    f_data%sim(simnr)%sim_setup%load_error_scale(:,kstep), &
                                    f_data%sim(simnr)%sim_setup%disp_error_scale(:,kstep), &
                                    error, e_cnt, evec_cnt, evec)
            endif ! End cyclic_error
        endif   ! End err_in_step
    endif   ! End kstep==nstep
    
    if (err_in_step.and.cyclic_error) then
        ! Export error and reset counters if last in cycle
        call get_cyclic_error(f_data%sim(simnr)%sim_setup%steps(kstep:min(kstep+1,nstep)), f_data%sim(simnr)%err, stp_ctrl, &
            f_data%sim(simnr)%sim_setup%load_error_scale(:,kstep), f_data%sim(simnr)%sim_setup%disp_error_scale(:,kstep), &
            stp_err_scale, last_cyc_start_step, cyc_ptcnt, cyc_fsim, cyc_fexp, cyc_time, error, e_cnt, evec_cnt, evec)
        
    endif
    
    deallocate(free_dofs)
    
    deallocate(tmain)
    
end do STEP_LOOP

if (error<(huge(1.d0)/10)) then ! Check that simulation succeeded before setting end values
    evec_cnt = evec_cnt - 1                 ! Correct to get n_errors
    f_data%sim(simnr)%n_errors = evec_cnt   ! Update n_errors
    
    error = error/(1.d-20+e_cnt)         ! Normalize with number of calculated increments (avoid division by zero if scale is zero for entire sim)
    
    if (resize_evec) then
        allocate(evec_tmp(max(evec_cnt,1)))
        evec_tmp = evec(1:max(evec_cnt,1))
        deallocate(evec)
        allocate(evec(max(evec_cnt,1)))
        evec = evec_tmp
        deallocate(evec_tmp)
    endif
endif

call system_clock ( clock_count, clock_rate)
stoptime = real(clock_count)/real(clock_rate)

call close_result_output(save_results, stoptime-starttime, error)


end subroutine mp_simulate
    
end module mps_mod  