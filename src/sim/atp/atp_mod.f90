
module atp_mod
use gen_util_mod
use atp_util_mod
use atp_element_mod
use atp_import_mod
use sim_error_mod
use sim_getincr_mod
use sim_writing_util_mod
use types_mod
use umat_mod
use output_mod
implicit none

    private
    public  :: atp_simulate         ! Main axial-torsion-pressure simulation routine
    
    contains
    
subroutine atp_simulate(error, props, f_data, simnr, evec)
implicit none
! Input parameters
double precision, intent(out)   :: error            ! Objective function value output
double precision, intent(in)    :: props(:)         ! Material parameters (true, full) input
type(fdata_typ)                 :: f_data           ! Settings input
integer                         :: simnr            ! Simulation number (=> f_data%sim(simnr) is pertinent to current simulation)
double precision, optional, allocatable :: evec(:)  ! Error vector


! Variables from global settings
logical                         :: save_results     ! Should results be saved to file?
procedure(umat_template),pointer:: umat_address     ! Addresss to umat subroutine
character(len=strl)             :: cmname           ! Material model name
logical                         :: nlgeom           ! Geometrically nonlinear analysis?

! Variables derived from simulation settings
logical                         :: result_onlymain  ! Only include results from main increments, or from all
logical                         :: result_inclexp   ! Include experiment data in result file?
double precision, allocatable   :: result_steps(:)  ! Which steps to write results for (if = 0 then all, if =-1 then none)
double precision, allocatable   :: error_steps(:)   ! Which steps to save error for (if =0 then all)
        
! Geometry and mesh
double precision, allocatable   :: rpos(:,:)        ! Radial position of each node for each element
double precision                :: h0               ! Height of considered gauge section
integer                         :: ngp              ! Number of gauss points in each element
logical                         :: bbar             ! Use the Abaqus bbar method or not

! Loading
integer                         :: exp_info(11)     ! Columns in exp_data. If 0 then column does not exist (and will be set to zero). 
                                                    ! [Step, time, forc, astr, torq, rota, p_i, cstr_i, p_o, cstr_o, temp] (cstr=circumferential strain)
! Iteration settings
double precision, allocatable   :: iter_err_norm(:) ! Geometry factors for converting equilibrium error to a stress measure
type(iter_typ)                  :: iter             ! Iteration settings

! Internal variables
 ! Error calculation
integer                         :: e_cnt            ! Number of times the error is evaluated
double precision, allocatable   :: err_tim_hist(:,:)! Step and time for use to calculate error
double precision, allocatable   :: err_exp_hist(:,:)! Experiment data for use to calculate error
double precision, allocatable   :: err_sim_hist(:,:)! Simulated data for use to calculate error
integer, allocatable            :: err_hist_comp(:) ! Colums in err_[exp/sim]_hist for disp values (load values at pos - 1)

 ! Time stepping
double precision                :: stp_time_incr(4) ! Current value for time_incr (dtmain, dtmin, dt0, dtmax)
integer                         :: stp_ctrl(4)      ! Current value for ctrl
double precision                :: stp_time_tot     ! Total time for current step
integer                         :: kstep, kinc, convincr
double precision                :: step


 ! Experimental data
double precision, allocatable   :: expdata(:,:)     ! Experiment data

 ! Other internals
double precision                :: time(2), pnewdt, dt, starttime, stoptime, dt_old
double precision                :: dtmain, dtmin, dtmax
double precision                :: disp_conv(4)    ! Convert from strain inputs to displacements used in FE solution
double precision                :: disp(4), disp_new(4), disp_exp(4)  ![eps_z, phi, eps_ci, eps_co] (c=circumferential)
double precision                :: load(4), load_new(4), load_exp(4)  ![Force, Torque, p_i, p_o]
double precision                :: temp, temp_new           !Temperature
double precision                :: denergy(4)       ! sse, spd, sce, rpl integrated over the volume
double precision                :: adjustment(2,4)  ! Adjustment variable for smooth following of experimental data if control mode is changed
double precision                :: len_adj          ! Adjustment of heigh (1) due to previous analysis (if cont_analysis<0) (otherwise 1.d0)
double precision, allocatable   :: gp_F(:,:,:), gp_F0(:,:,:), gp_s(:,:,:), gp_s0(:,:,:), gp_sv(:,:,:), gp_sv0(:,:,:)
double precision, allocatable   :: gp_strain(:,:,:), gp_strain0(:,:,:), additional_output(:)
double precision, allocatable   :: u(:), u0(:), du(:), v(:), v_old(:), tmain(:)
double precision, allocatable   :: stp_time_vec(:)
integer, allocatable            :: free_dofs(:), known_dofs(:), stprows(:)
logical                         :: laststep, lconv, res_in_step, err_in_step
logical                         :: adj_geom ! Should the geometry be adjusted due to previous analysis?
integer                         :: nprops, nstatv, nel, nnod, nstep, niter, ndof_tot, kmain
integer                         :: exp_row    ! Row in experiment data
integer                         :: clock_count, clock_rate


call system_clock ( clock_count, clock_rate)
starttime = real(clock_count)/real(clock_rate)

! == Material ==
umat_address => f_data%glob%umat_address
cmname       = f_data%glob%cmname
nlgeom       = f_data%glob%nlgeom
nprops       = size(props)

! == Use shorter names for some variables == 
! slight performance loss, but increases readability
! Experiment
exp_info = f_data%sim(simnr)%exp%exp_info
allocate(expdata, source=f_data%sim(simnr)%sim_setup%expdata_array)
allocate(stprows, source=f_data%sim(simnr)%sim_setup%stprows)
nstep = size(stprows)-1
! Iteration settings
iter = f_data%sim(simnr)%iter

! Error setup, should allocate error_hist and error_steps
call error_settings(f_data%sim(simnr)%err, err_tim_hist, err_exp_hist, err_sim_hist, err_hist_comp, &
                        error_steps, expdata(size(expdata,1), exp_info(2)), minval(iter%time_incr(2,:)))
e_cnt = 0       ! Initiate e_cnt to zero (number of error evaluations)
error = 0.d0    ! Set error to zero

! Result output
allocate(result_steps, source=f_data%sim(simnr)%outp%result_steps)
save_results    = (f_data%glob%resnr>0).and.(any(result_steps>-1.d-10))
result_onlymain = f_data%sim(simnr)%outp%result_onlymain
result_inclexp  = f_data%sim(simnr)%outp%result_inclexp

! == Import initial conditions (and adjust experiment data if continued analysis of separate experiments) ==
call atp_import_init(f_data, simnr, exp_info, nstatv, gp_s0, gp_strain0, gp_F0, gp_sv0, u0, temp, disp, load, disp_exp, load_exp, time, expdata, len_adj)

! == Mesh and geometry == 
call atp_import_mesh(f_data%sim(simnr)%mesh1d, len_adj, u0, h0, ngp, bbar, nnod, nel, ndof_tot, rpos, disp_conv, iter_err_norm)
call element_setup(ngp, nnod, bbar, simnr, set_full_output=(f_data%sim(simnr)%outp%log_output>=2))

!Allocate gp_, dispvars, additional_output
allocate(u(ndof_tot), du(ndof_tot), v(ndof_tot), v_old(ndof_tot))
allocate(gp_s(6, ngp, nel))
allocate(gp_strain(6, ngp, nel))
allocate(gp_F(9, ngp, nel))
allocate(gp_sv(nstatv, ngp, nel))
allocate(additional_output(1)); additional_output(1) = 0.d0

! RUN SIMULATION
! ============================================================================================================

!Define initial values
kstep = 0; kinc  = 0; niter = 0;
exp_row = 2 ! Start at row 2 to allow interpolation (if time=t_start then interpolated value will be to row 1 anyway)

! Output file
if (save_results) then
    call setup_result_output(trim(f_data%glob%outname)//'_sim'//int2str(simnr)//'_'//int2str(f_data%glob%resnr)//'.txt', 1, nlgeom, f_data%sim(simnr)%outp, ngp)
endif

STEP_LOOP: do kstep = 1,nstep
    step = f_data%sim(simnr)%sim_setup%steps(kstep)
    ! Settings for current step
    call get_step_data_dbl(stp_time_incr, iter%time_incr, step)
    call get_step_data_int(stp_ctrl, f_data%sim(simnr)%exp%ctrl, step)
    
    ! Initialize current step (reset variables and read in new)
    convincr   = 0;
    dtmain     = stp_time_incr(1); dtmin = stp_time_incr(2)
    dt         = stp_time_incr(3); dtmax = stp_time_incr(4)
    v          = 0.d0 !Should reset this guess for each new step
    time(1)    = 0.d0; kinc       = 0
    
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
    call gen_free_dofs(free_dofs, known_dofs, ndof_tot, stp_ctrl, f_data%sim(simnr)%exp%intpres_extstrn)
    
    ! Find adjustments if any to the experiment data 
    ! E.g. to offset to match initial calculated value in step, 
    ! or offset but scale to get same end value as the experiment data
    adjustment = get_incr_adj(load, disp, expdata(stprows(kstep):stprows(kstep+1), :), exp_info, stp_ctrl)
    
    TMAIN_LOOP: do kmain=1,size(tmain)  ! Loop over all main increments
        laststep = .false.
        
        if (err_in_step) then
            e_cnt = e_cnt + 1
            call update_error(load_exp, load, disp_exp, disp, step, time(2), err_hist_comp, e_cnt, err_tim_hist, err_exp_hist, err_sim_hist)
        endif
        
        ! Write out results for main step if it should for current step
        if (res_in_step) then
            call write_result(step, kinc, niter, time(2), temp, load, load_exp, disp, disp_exp, f_data%sim(simnr)%outp, gp_sv0, gp_s0, gp_strain0, gp_F0, u0(3:size(u)), denergy)
        endif
        
        do while (.not.laststep)
            ! Set starting variables for current increment
            gp_F = gp_F0; gp_s = gp_s0; gp_sv = gp_sv0; gp_strain = gp_strain0
            pnewdt    = 1.d0
            u = u0
        
            ! Update time and sought stress/displacement for this increment
            call timeupdate(kinc, time, dt, laststep, tmain(kmain), dtmin)
        
            ! Find load_exp, disp_exp and temp_new
            call get_incr(temp, load_exp, disp_exp, temp_new, time(2), expdata, exp_info, exp_row)
            ! Determine load_new and disp_new (adjusted if requested by ctrl values)
            call adjust_incr(load_exp, disp_exp, load_new, disp_new, stp_ctrl, adjustment, time(2))
            
            ! Apply boundary conditions
            call apply_bc(disp, disp_new, du, stp_ctrl, disp_conv)

            ! Make a guess for the displacement
            call update_guess(du, v, dt, v_old, dt_old, kinc, free_dofs)
        
            ! Solve increment
            call solve_incr(rpos, h0, load_new, disp_new, temp, temp_new-temp, time, dt, free_dofs, known_dofs, gp_F, gp_s, gp_sv, gp_strain, u, du, iter, niter, lconv, pnewdt, &
            kinc, kstep, props, cmname, umat_address, nlgeom, iter_err_norm, denergy)
            ! Check results and if OK go to next step
            if (lconv) then
                load = load_new; disp = disp_new; temp = temp_new
                gp_F0 = gp_F;   gp_s0 = gp_s;   gp_sv0 = gp_sv; gp_strain0 = gp_strain
                u0    = u;  v_old = v; v     = du/dt; dt_old= dt

                ! Write out results if that is requested for current simulation/step/increment for steps that are not main steps
                if ((res_in_step).and.(.not.result_onlymain).and.(.not.laststep)) then
                    call write_result(step, kinc, niter, time(2), temp, load, load_exp, disp, disp_exp, f_data%sim(simnr)%outp, gp_sv, gp_s, gp_strain, gp_F, u(3:size(u)), denergy)
                endif
                
                ! Update time stepping
                convincr = convincr + 1
                if (convincr>iter%nconv_incr) then !Increase dt if enough converged iterations has precided this increment
                    dt = min(dt*iter%dt_incr, dtmax)
                endif
            else
                if (dt<=dtmin*(1.d0+1.d-5)) then ! Note that this tolerance must be larger than dt_last_tol in timeupdate!
                    call write_output('Time increment smaller than dtmin requested, exiting simulation at increment '//int2str(kinc), 'warning', 'sim:atp')
                    call write_output('Simulation nr '//int2str(simnr)//' in step '//dbl2str(step), 'warning', 'sim:atp', loc=.false.)
                    call write_output('Setting error to a huge number', 'warning', 'atp', loc=.false.)
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
            call write_result(step, kinc, niter, time(2), temp, load, load_exp, disp, disp_exp, f_data%sim(simnr)%outp, gp_sv, gp_s, gp_strain, gp_F, u(3:size(u)), denergy)
        endif
        
        if (err_in_step) then
            e_cnt = e_cnt + 1
            call update_error(load_exp, load, disp_exp, disp, step, time(2), err_hist_comp, e_cnt, err_tim_hist, err_exp_hist, err_sim_hist)
        endif   ! End err_in_step
    endif   ! End kstep==nstep

    deallocate(free_dofs, known_dofs, tmain)
    
end do STEP_LOOP

if (error<(huge(1.d0)/10)) then ! Check that simulation succeeded before calculating error and setting end values
    if (e_cnt>0) then
        call calculate_error(f_data%sim(simnr)%err, f_data%sim(simnr)%exp%ctrl, &
                             err_tim_hist, err_exp_hist, err_sim_hist, err_hist_comp, e_cnt, error, evec)
    endif
    
    ! Export end values
    call atp_export_end(f_data, simnr, u, gp_s, gp_strain, gp_F, gp_sv, disp, load, h0)

endif

if (f_data%sim(simnr)%outp%log_output==1) then
    if (get_material_didnt_converge_count()>0) then
        call write_output('material routine requested smaller timestep '//int2str(get_material_didnt_converge_count())//' times during simulation nr '//int2str(simnr), 'status', 'sim:atp')
    endif
endif

call system_clock ( clock_count, clock_rate)
stoptime = real(clock_count)/real(clock_rate)

call close_result_output(save_results, stoptime-starttime, error)

end subroutine atp_simulate

end module 

