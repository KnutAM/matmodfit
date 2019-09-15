module sim_getincr_mod
    use output_mod
    implicit none
    
    private
    
    public  :: get_step_data_dbl, get_step_data_int     ! Obtain data of "step data type"
    public  :: get_incr_adj                             ! Get the adjustment for the load and displacement increments to smoothly work when control type changes
    public  :: get_incr                                 ! Get the experiment load, displacement and temperature values for current increment
    public  :: adjust_incr                              ! Adjust the simulated control variables in disp and load according to smoothly change control mode.
    public  :: get_tmain                                ! Determine the main time increments to match the experiment data points
    public  :: timeupdate                               ! Update kinc, time and potentially dt and laststep
    public  :: update_guess                             ! Make a qualified guess for the free degrees of freedom
    
    
    integer, parameter, dimension(4):: load_ind = (/1,3,5,7/)   !NOTE:  Does not suffice to change these here if convention updated,
    integer, parameter, dimension(4):: disp_ind = (/2,4,6,8/)   !       this will not result in a correct program...
    
    contains 
    
subroutine get_step_data_dbl(stp_var, glob_var, stpnr)
    implicit none
    double precision, intent(inout) :: stp_var(:)
    double precision, intent(in)    :: glob_var(:,:)
    double precision, intent(in)    :: stpnr
    integer                         :: k1
    
    k1 = size(glob_var,2)

    do while (k1 > 1)
        if (stpnr>(glob_var(1,k1)-1.d-10)) then
            exit
        endif
        k1 = k1 - 1
    enddo
    
    stp_var = glob_var(2:size(glob_var, 1),k1)
    
end subroutine

subroutine get_step_data_int(stp_var, glob_var, stpnr)
    implicit none
    integer, intent(inout)          :: stp_var(:)
    double precision, intent(in)    :: glob_var(:,:)
    double precision, intent(in)    :: stpnr
    integer                         :: k1

    k1 = size(glob_var,2)

    do while (k1 > 1)
        if (stpnr>(glob_var(1,k1)-1.d-10)) then
            exit
        endif
        k1 = k1 - 1
    enddo
    
    stp_var = nint(glob_var(2:size(glob_var, 1),k1))
    
end subroutine

function get_incr_adj(load, disp, expdata_step, exp_info, ctrl) result(adjustment)
! Find adjustment factors to ensure consistent continuation if control mode changes
implicit none
double precision, intent(in)    :: load(:), disp(:), expdata_step(:,:)
integer, intent(in)             :: exp_info(:), ctrl(:)
double precision, allocatable   :: adjustment(:,:)
integer                         :: k1
double precision                :: vc1, ve1     !Calculated and experiment values at start
double precision                :: times(2)     !(/time1, time2/) (Time at start and end of step)

times = expdata_step((/1, size(expdata_step,1)/), exp_info(2))
allocate(adjustment(2,size(ctrl)))
adjustment = 0.d0

do k1=1,size(ctrl)
    if (abs(ctrl(k1))>1) then ! Adjustment requested
        if (ctrl(k1)<0) then ! Load controlled
            vc1 = load(k1)  ! Calculated value at time 1
            ve1 = expdata_step(1, exp_info(load_ind(k1)+2)) ! Experiment value at time 1
        else                    ! Displacement controlled
            vc1 = disp(k1)  ! Calculated value at time 1
            ve1 = expdata_step(1, exp_info(disp_ind(k1)+2)) ! Experiment value at time 1
        endif
        
        if (abs(ctrl(k1))==2) then  !Shift exp data to match calculated start
            adjustment(1,k1) = vc1-ve1
        elseif (abs(ctrl(k1))==3) then  !Scale exp data to match calculated start and exp data end
            ! c(t) = e(t) + ((t2-t)/(t2-t1))*(s1-e1) = e(t) + (t2/(t2-t1))*(s1-e1) + t*(-1/(t2-t1))*(s1-e1)
            adjustment(1,k1) = (times(2)/(times(2)-times(1)))*(vc1-ve1)
            adjustment(2,k1) = (-1.d0/(times(2)-times(1)))*(vc1-ve1)
        else
            call write_output('abs(ctrl)>3 not supported', 'error', 'atp')
        endif
    endif
    
enddo


end function

subroutine get_incr(temp, load_exp, disp_exp, temp_new, total_time, expdata, exp_info, fpos)
implicit none
    double precision, intent(in)    :: expdata(:,:)     ! Columns according to exp_info
    double precision, intent(in)    :: total_time
    double precision, intent(in)    :: temp             ! Temperature
    integer, intent(in)             :: exp_info(:)      ! expdata info: (stp, time, chA_load, chA_disp, chB_load, ..., temp)
    double precision, intent(inout) :: load_exp(:)      ! Load for each channel
    double precision, intent(inout) :: disp_exp(:)      ! Disp for each channel
    double precision, intent(inout) :: temp_new         ! Updated temperature
    integer, intent(inout)          :: fpos             ! Current position in expdata matrix
    ! Internal variables
    integer                         :: k1, k2           ! Iteration variables
    integer                         :: nmaxcol          ! Length of exp_info
    double precision                :: val, interpdata(2,2)  ! Interpolation data: (:,1)=time, (:,2) yval
    logical                         :: use_fpos         ! Logical variable to check if current file position should be used
    double precision                :: time_extrap      ! Length of extrapolated time
    
    ! Initialize to zero for experiment load and disp (default values if not specified in experiment)
    load_exp = 0.d0
    disp_exp = 0.d0
    ! Keep old temperature
    temp_new = temp
    nmaxcol = size(exp_info)
    ! Go backwards if we have gone too far (this can happen if the time step is reduced!)
    ! We never allow fpos=1, at least two is needed for interpolation. 
    do while ((expdata(fpos,exp_info(2))>=total_time).and.(fpos>2))
        fpos = fpos - 1
    enddo
    
    ! Loop until we come to a point where the experiment time is larger than the current total time
    do k1=fpos,size(expdata,1)
        use_fpos = .false.
        if (expdata(k1,exp_info(2))>=total_time) then   ! Experiment time larger than current time
            use_fpos = .true.
        elseif (k1==size(expdata,1)) then               ! Last experiment time point (extrapolate, but warn if more than numerically small time difference)
            use_fpos = .true.
            time_extrap = total_time-expdata(k1,exp_info(2))
            if ((time_extrap/total_time)>1e-9) then
                call write_output('Extrapolating last time step by '//dbl2str(time_extrap)//' s', 'warning', 'atp')
            endif
        endif
        
        if (use_fpos) then
            do k2=3,nmaxcol
                if (exp_info(k2)>0) then
                    interpdata = expdata((/k1-1, k1/), (/exp_info(2), exp_info(k2)/))
                    val = interp(total_time, interpdata)                    
                    if (k2==nmaxcol) then
                        temp_new = val
                    elseif (isodd(k2)) then
                        load_exp((k2-1)/2) = val ! k2=[3, 5, ...]
                    else
                        disp_exp((k2-2)/2) = val ! k2=[4, 6, ...]
                    endif
                endif
            enddo
            exit    ! We found our values, and can now exit the outer loop
        endif
    enddo
    
    fpos = k1
    
end subroutine get_incr

subroutine adjust_incr(load_exp, disp_exp, load_new, disp_new, ctrl, adjustment, total_time)
! Adjust load and displacement increments according to previously determined adjustment and the adjustment type specified by ctrl for current step
implicit none
double precision, intent(in)    :: load_exp(:), disp_exp(:), adjustment(:,:), total_time
integer, intent(in)             :: ctrl(:)
double precision, intent(inout) :: load_new(:), disp_new(:)
integer                         :: k1

do k1=1,size(ctrl)
    if (ctrl(k1)==-1) then      ! No adjustment, load controlled
        load_new(k1) = load_exp(k1)
    elseif (ctrl(k1)==1) then   ! No adjustment, disp controlled
        disp_new(k1) = disp_exp(k1)
    elseif (ctrl(k1)<0) then    ! Load controlled
        load_new(k1) = adjust(load_exp(k1), total_time, abs(ctrl(k1)), adjustment(:,k1))
    else                        ! Disp controlled
        disp_new(k1) = adjust(disp_exp(k1), total_time, ctrl(k1), adjustment(:,k1))
    endif
enddo

end subroutine
    
function adjust(val, time, ctrl_check, adjustment) result(val_adjusted)
implicit none
double precision, intent(in)    :: val, time
integer, intent(in)             :: ctrl_check
double precision, intent(in)    :: adjustment(2)
double precision                :: val_adjusted

if (ctrl_check<=1) then
    val_adjusted = val
elseif (ctrl_check==2) then
    val_adjusted = val + adjustment(1)
elseif (ctrl_check==3) then
    val_adjusted = val + adjustment(1) + adjustment(2)*time
else
    call write_output('abs(ctrl)>3 not supported', 'error', 'atp')
endif

end function
        
subroutine get_tmain(tmain, stp_time_tot, texp, dtmain)
implicit none
    double precision, allocatable   :: tmain(:)     ! The actual main time increments in current step
    double precision, intent(out)   :: stp_time_tot ! Total step time
    double precision, intent(in)    :: texp(:)      ! Array of all experiment data time points in current step
    double precision, intent(in)    :: dtmain       ! The requested main time increment
    
    integer                         :: nmain        ! Number of main increments within current step
    double precision, allocatable   :: tmp_tmain(:) ! Temporary storage until final size decided.
    integer                         :: k1, k2, i    ! Counters
    integer                         :: minind, minind_old   ! Index of found minimum distance between main time increment and real.
    
    
    stp_time_tot = texp(size(texp))-texp(1)
    nmain = int(stp_time_tot/dtmain + 1 - 1.d-8)  ! Round upwards the expected number of main time increments
    allocate(tmp_tmain(nmain))
    do i = 1,nmain
        tmp_tmain(i) = i*(stp_time_tot/nmain) + texp(1)
    enddo
    !Modify last entry to avoid numerical issues when checking if experiment time larger than main time. 
    !Note that this doesn't change the final value as this is given by the experiment times
    tmp_tmain(nmain) = tmp_tmain(nmain) - dtmain*1.d-9
    !tmp_tmain =  (/((i*(stp_time_tot/nmain)+texp(1)),i=1,nmain)/)   ! Create vector of main times
    
    minind = 1  ! We do not accept the first index as a main step, because this is taken care of by the last. 
                ! If this is desired, one has to modify the timeupdate to avoid getting stuck in a zero increment...
    minind_old = 1
    
    ! We need to ensure that the tmain hits exactly at the experiment data points 
    ! (because we don't want to evaluate the error on interpolated data!)
    ! This could cause double entries (if the experiment data points are further appart than dtmain)
    k2 = 0
    do k1=1,nmain
        minind_old = minind_old+minind-1
        i=0
        do while(texp(minind_old+i)<tmp_tmain(k1))
            i = i + 1
        enddo
        minind = i-1+ minloc( abs( texp((minind_old+i-1):(minind_old+i))-tmp_tmain(k1)),1)
        if (minind>1) then ! Only add if the current time has not been added before
            k2 = k2 + 1
            tmp_tmain(k2) = texp(minind_old+minind-1)
        endif
    enddo
    
    allocate(tmain(k2))
    tmain = tmp_tmain(1:k2)-texp(1)
end subroutine

subroutine timeupdate(kinc, time, dt, laststep, main_time, dtmin)
!Update kinc, time and potentially dt and laststep
        implicit none
        double precision, parameter:: dt_last_tol   = 1.d-6 ! Scaled with dtmin
        double precision :: time(2), dt, main_time, dtmin
        integer          :: kinc
        logical          :: laststep

        kinc = kinc + 1
        
        if ((time(1)+dt)>(main_time-dtmin*dt_last_tol)) then    !Attempt to go further than main time, reduce dt to hit main_time. 
            dt = main_time-time(1)                              !A small tolerance scaled with the current time step is used to 
            time(1) = main_time                                 !avoid numerical problems with the comparison.
            time(2) = time(2) + dt
            laststep = .true.
        elseif ((main_time-(time(1)+dt))<dtmin) then    ! This means that the next time step will need to be smaller than dtmin.
            dt = (main_time-time(1))/2                  ! We therefore make these two last steps equal
            time = time + dt
        else    !Business as usual...
            time = time + dt
        endif

end subroutine timeupdate

subroutine update_guess(du, v, dt, v_old, dt_old, kinc, free_dofs)
    implicit none
    double precision, intent(in)    :: v(:), v_old(:), dt, dt_old
    integer, intent(in)             :: kinc, free_dofs(:)
    double precision                :: du(:)
! Have done a simple check if there is any advantage to the 2nd order guess, and there seems to be a slight advantage. 
    
    if (free_dofs(1)/=-1) then      ! If not all prescribed, assign guess
        if (kinc==1) then
            du(free_dofs) = 0.d0
        elseif (kinc==2) then
            du(free_dofs) = v(free_dofs)*dt
        else
            du(free_dofs) = v(free_dofs)*dt + ((v(free_dofs)-v_old(free_dofs))*dt**2)/(dt_old+dt)
        endif
    endif
    
end subroutine

! Internal utility functions
function interp(x, dmatrix) result (y)
implicit none
    double precision    :: x, dmatrix(2,2)
    double precision    :: y
    double precision    :: xval(2), yval(2)
    ! dmatrix contains x values in first column and y values in second column
    xval = dmatrix(:,1)
    yval = dmatrix(:,2)
    y = (yval(2)-yval(1))*(x - xval(1))/(xval(2)-xval(1)) + yval(1)
    
end function interp

function isodd(int)
! Return true if int is odd, and false if int is even (including zero)
implicit none
    integer :: int
    logical :: isodd
    
    isodd = .not.(mod(int,2)==0)
    
end function

end module