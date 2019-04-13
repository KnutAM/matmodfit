module sim_error_mod
use gen_util_mod
use types_mod
use output_mod
use sim_getincr_mod
use usr_interface_mod
    implicit none
    
    private
    public error_settings   ! Read in error settings and allocate error_hist variable
    public update_error     ! Append errors to the error_hist variable
    public calculate_error  ! Calculate the simulation error based on error_hist
    
    contains

subroutine error_settings(err, err_tim_hist, err_exp_hist, err_sim_hist, err_hist_comp, error_steps, tmax, dtmain_min)
implicit none
    type(err_typ)                   :: err                      ! Custom type containing error settings
    double precision, allocatable   :: err_tim_hist(:,:)        ! Keep track of step and time for each error contribution
    double precision, allocatable   :: err_exp_hist(:,:)        ! Experimental values used to calculate each error contribution
    double precision, allocatable   :: err_sim_hist(:,:)        ! Simulated values used to calculate each error contribution
    integer, allocatable            :: err_hist_comp(:)         ! Describes where in err_[exp/sim]_hist to put the displacement values (Loads are put at pos-1). 0=> do not put
    double precision, allocatable   :: error_steps(:)           ! Which steps to calculate the error for
    double precision, intent(in)    :: tmax, dtmain_min         ! Total simulation time and the smallest dtmain (used to calculate maximum number of rows in error_hist if unknown)
    
    integer                         :: hist_rows, hist_cols, k1
    
    ! Determine number of rows in error_hist (i.e. how many times will an error be calculated)
    if (err%hist_rows == -1) then            ! Unknown (first time simulation is run)
        hist_rows = int(tmax/dtmain_min)*2+1 ! Multiplying with ratio+1 with 2 should ensure that it can never go above. Usually this will be way too big, but only runs once so should be ok
    else                                     ! Known (simulation has been run before)
        hist_rows = err%hist_rows
    endif
    
    ! Determine number of columns in error_hist
    ! err_scale(chnl, stp_spec)
    allocate(err_hist_comp(size(err%err_scale, 1)-1))
    err_hist_comp = 0
    
    hist_cols = 0
    do k1=1,size(err%err_scale, 1)-1  ! Loop over each channel in err_scale
        if (any(err%err_scale(k1+1, :)>0.d0).or.any(err%err_scale_ctrl(k1+1, :)>0.d0)) then
            hist_cols = hist_cols + 2
            err_hist_comp(k1) = hist_cols
        endif
    enddo
    
    allocate(err_tim_hist(hist_rows, 2), err_exp_hist(hist_rows, hist_cols), err_sim_hist(hist_rows, hist_cols))
    err_tim_hist = 0.d0
    err_exp_hist = 0.d0
    err_sim_hist = 0.d0
    
    allocate(error_steps,  source=err%err_steps)
    
end subroutine

    
subroutine update_error(load_exp, load, disp_exp, disp, step, time, err_hist_comp, e_cnt, err_tim_hist, err_exp_hist, err_sim_hist)
implicit none
    double precision, intent(in)    :: load_exp(:), load(:), disp_exp(:), disp(:)
    double precision, intent(in)    :: step, time
    integer, intent(in)             :: err_hist_comp(:) ! Which column to put each channel in [L1, D1, L2, D2, ...]. 0=> Do not put anywhere
    integer, intent(in)             :: e_cnt ! Row in error_hist
    double precision, intent(inout) :: err_tim_hist(:,:), err_exp_hist(:,:), err_sim_hist(:,:)
    
    integer                         :: k1, pos(2)
    
    err_tim_hist(e_cnt, 1) = step
    err_tim_hist(e_cnt, 2) = time
    
    do k1=1,size(err_hist_comp)
        if (err_hist_comp(k1).ne.0) then
            pos = err_hist_comp(k1) + [-1, 0]
            err_exp_hist(e_cnt, pos) = [load_exp(k1), disp_exp(k1)]
            err_sim_hist(e_cnt, pos) = [load(k1), disp(k1)]
        endif
    enddo
    
end subroutine

subroutine calculate_error(err, ctrl, err_tim_hist, err_exp_hist, err_sim_hist, err_hist_comp, e_cnt, error, evec)
implicit none
    type(err_typ)                           :: err                      ! Custom type containing error settings
    double precision, intent(in)            :: ctrl(:,:)                ! Control as step variable type (i.e. first row represents which step it the remaining starts being valid from)
    double precision, intent(inout)         :: err_tim_hist(:,:), err_exp_hist(:,:), err_sim_hist(:,:)
    integer                                 :: err_hist_comp(:)         ! Describes where in err_[sim/exp]_hist the disp values are put (load values are put at pos-1)
    integer                                 :: e_cnt                    ! Counter for how many errors have been calculated
    double precision, intent(out)           :: error                    ! Calculated error
    double precision, allocatable, optional :: evec(:)                  ! Error vector
    
    ! Update size of error_hist
    if (err%hist_rows == -1) then            ! Unknown (first time simulation is run)
        err%hist_rows = e_cnt
    elseif (err%hist_rows.ne.e_cnt) then
        call write_output('Error count doesn''t match from previous simulation (please save input files and report to developers)', 'warning', 'sim:err')
        call write_output('Error calculation may be incorrect', 'warning', 'sim:err', loc=.false.)
        err%hist_rows = e_cnt
    endif
    
    ! Update scalar error
    select case (err%error_type)
    case (0)  ! User written subroutine for calculating error and potentially evec
        call user_error(err, ctrl, e_cnt, err_tim_hist, err_exp_hist, err_sim_hist, err_hist_comp, error, evec)
    case (1)  ! Default type, square sum of error
        call square_error(err, ctrl, e_cnt, err_tim_hist, err_exp_hist, err_sim_hist, err_hist_comp, error, evec)
    case (2)  ! Cyclic error
        call cyclic_error(err, ctrl, e_cnt, err_tim_hist, err_exp_hist, err_sim_hist, err_hist_comp, error, evec)
    case default
        call write_output('error_type = '//int2str(err%error_type)//' is not supported', 'error', 'sim:err')
    end select
    
end subroutine
    
subroutine user_error(err, ctrl, e_cnt, err_tim_hist, err_exp_hist, err_sim_hist, err_hist_comp, error, evec)
implicit none
    type(err_typ)                           :: err                      ! Custom type containing error settings
    double precision, intent(in)            :: ctrl(:,:)                ! Control information
    double precision, intent(inout)         :: err_tim_hist(:,:), err_exp_hist(:,:), err_sim_hist(:,:)
    integer                                 :: err_hist_comp(:)         ! Describes where in err_[sim/exp]_hist the disp values are put (load values are put at pos-1)
    integer                                 :: e_cnt                    ! Counter for how many errors have been calculated
    double precision, intent(out)           :: error                    ! Calculated error
    double precision, allocatable, optional :: evec(:)                  ! Error vector
    procedure(user_error_template),pointer  :: usr_error ! Address to user subroutine
    
    usr_error => err%error_address
    
    if (present(evec)) then
        call usr_error(err%user_settings, ctrl, err_tim_hist(1:e_cnt,:), err_exp_hist(1:e_cnt,:), err_sim_hist(1:e_cnt,:), err_hist_comp, error, evec)
    else
        call usr_error(err%user_settings, ctrl, err_tim_hist(1:e_cnt,:), err_exp_hist(1:e_cnt,:), err_sim_hist(1:e_cnt,:), err_hist_comp, error)
    endif
    
end subroutine

subroutine square_error(err, ctrl, e_cnt, err_tim_hist, err_exp_hist, err_sim_hist, err_hist_comp, error, evec)
implicit none
    type(err_typ)                           :: err                      ! Custom type containing error settings
    double precision, intent(in)            :: ctrl(:,:)                ! Control information
    double precision, intent(inout)         :: err_tim_hist(:,:), err_exp_hist(:,:), err_sim_hist(:,:)
    integer                                 :: err_hist_comp(:)         ! Describes where in err_[sim/exp]_hist the disp values are put (load values are put at pos-1)
    integer                                 :: e_cnt                    ! Counter for how many errors have been calculated
    double precision, intent(out)           :: error                    ! Calculated error
    double precision, allocatable, optional :: evec(:)                  ! Error vector
    
    double precision                        :: step, error_scale_factor, ch_ind
    integer                                 :: k1, k2, thisind, nextind, tempind, num_channels
    double precision, allocatable           :: scaling(:), stp_err_scale(:), stp_err_scale_ctrl(:)
    integer, allocatable                    :: stp_ctrl(:)
    
    num_channels = size(ctrl,1)-1    ! Number of channels
    allocate(stp_ctrl(num_channels), stp_err_scale(num_channels), stp_err_scale_ctrl(num_channels))
    
    ! Calculate global scaling factors if scaling with max-min in experiment
    if (err%err_norm_met==2) then 
        allocate(scaling(size(err_exp_hist,2)))
        do k1=1,size(scaling)
            scaling(k1) = maxval(err_exp_hist(:,k1)) - minval(err_exp_hist(:,k1))
            if (scaling(k1)<10.d0/huge(1.d0)) then
                scaling(k1) = 0.d0  ! Have to avoid division by zero in cases where the experiment data is numerically equal to zero
            else
                scaling(k1) = 1.d0/scaling(k1)
            endif
        enddo
    endif
    
    ! Loop over each step in the error history saved. Treat each of these bulks by themselves as different steps may have different error scaling
    thisind = 1
    tempind = 1 ! Needs to be set to some value >0 to enter loop
    do while (tempind>0)
        step = err_tim_hist(thisind, 1)
        
        tempind = findfirst(err_tim_hist(thisind:, 1)>step)  ! Gives zero if no true values
        if (tempind==0) then
            nextind = size(err_tim_hist,1)+1
        else
            nextind = thisind + tempind - 1
        endif
        
        call get_step_data_int(stp_ctrl, ctrl, step)
        call get_step_data_dbl(stp_err_scale, err%err_scale, step)
        call get_step_data_dbl(stp_err_scale_ctrl, err%err_scale_ctrl, step)
        
        
        do ch_ind=1,num_channels
            if (err_hist_comp(ch_ind)>0) then
                do k2=0,1   ! Load first, then disp. (2*k2-1) = [-1, 1]
                    k1=err_hist_comp(ch_ind)-1+k2   ! First err_hist_comp(ch_ind)-1 (load), then err_hist_comp(ch_ind) (disp)
                    if (stp_ctrl(ch_ind)*(2*k2-1)>0) then    ! Controlled part ctrl and (2*k2-1) have the same sign
                        error_scale_factor = stp_err_scale_ctrl(ch_ind)
                    else    ! The part of the channel not controlled
                        error_scale_factor = stp_err_scale(ch_ind)
                    endif
                        
                    select case (err%err_norm_met)
                    case (0) ! No scaling
                        err_sim_hist(thisind:(nextind-1), k1) = error_scale_factor*(err_sim_hist(thisind:(nextind-1), k1) - err_exp_hist(thisind:(nextind-1), k1))
                    case (1) ! Scale with max - min in step
                        err_sim_hist(thisind:(nextind-1), k1) = error_scale_factor*(err_sim_hist(thisind:(nextind-1), k1) - err_exp_hist(thisind:(nextind-1), k1))/&
                                                                (maxval(err_exp_hist(thisind:(nextind-1), k1)) - minval(err_exp_hist(thisind:(nextind-1), k1)))
                    case (2) ! Scale with max - min in experiment
                        err_sim_hist(thisind:(nextind-1), k1) = error_scale_factor*(err_sim_hist(thisind:(nextind-1), k1) - err_exp_hist(thisind:(nextind-1), k1))*scaling(k1)
                    case default
                        call write_output('err_norm_met='//int2str(err%err_norm_met)//' is not supported', 'error', 'sim:err')
                    end select
                enddo
            endif
        enddo
        
        thisind = nextind
    enddo
    
    
    ! Calculate actual error
    error = sum(err_sim_hist**2)/e_cnt
        
    if (present(evec)) then
        if (.not.allocated(evec)) then
            allocate(evec(e_cnt*size(err_sim_hist,2)))
        endif
        do k1=1,size(err_sim_hist,2)
            evec(((k1-1)*e_cnt+1):(k1*e_cnt)) = err_sim_hist(1:e_cnt, k1)
        enddo
    endif
    

end subroutine
    



subroutine cyclic_error(err, ctrl, e_cnt, err_tim_hist, err_exp_hist, err_sim_hist, err_hist_comp, error, evec)
implicit none
    type(err_typ)                           :: err                      ! Custom type containing error settings
    double precision, intent(in)            :: ctrl(:,:)                ! Control information
    double precision, intent(inout)         :: err_tim_hist(:,:), err_exp_hist(:,:), err_sim_hist(:,:)
    integer                                 :: err_hist_comp(:)         ! Describes where in err_[sim/exp]_hist the disp values are put (load values are put at pos-1)
    integer                                 :: e_cnt                    ! Counter for how many errors have been calculated
    double precision, intent(out)           :: error                    ! Calculated error
    double precision, allocatable, optional :: evec(:)                  ! Error vector
    
    integer, allocatable                    :: cycle_inds(:)            ! Indices in error_hist (rows) pertaining to start and stop of cycles
    double precision                        :: mean_exp_cycle, mean_sim_cycle ! Mean values during each cycle
    integer                                 :: c1, c2                   ! Indicies for cycle start and end
    integer                                 :: k1, k2, k3, ch_ind, evec_cnt
    double precision, allocatable           :: dt(:)
    double precision                        :: cycle_time, step, error_scale_factor, tmp_error
    double precision, allocatable           :: stp_err_scale(:), stp_err_scale_ctrl(:)
    integer, allocatable                    :: stp_ctrl(:)
        
    allocate(stp_err_scale(size(ctrl,1)-1), stp_err_scale_ctrl(size(ctrl,1)-1), stp_ctrl(size(ctrl,1)-1))
    ! Step 1: Identify the cycles
    ! nstep_cyc_err_calc=1.d0  ! What difference in the step column builds up one cycle?
    ! nstep_cyc_initial=0.d0   ! What value in the step column before the first cycle?
    ! err_type_scale=0.5       ! Weight of cycle average error versus cycle shape error
    call determine_cycles(err_tim_hist(1:e_cnt,1), err%nstep_cyc_err_calc, err%nstep_cyc_initial, cycle_inds)
    
    if (present(evec)) then
        allocate(evec(size(err_sim_hist) + size(err_sim_hist,2)*(size(cycle_inds)-1)))
        evec = 0.d0
        evec_cnt = 1
    endif
    
    ! Step 2: Calculate the error for each cycle
    error = 0.d0
    do k1=1,size(cycle_inds)-1
        c1 = cycle_inds(k1)+1   ! First index of cycle
        c2 = cycle_inds(k1+1)   ! Last index of cycle
        step = err_tim_hist(c2,1)
        allocate(dt(c2-c1+1))
        dt(2:(c2-c1)) = (err_tim_hist((c1+2):(c2),2)-err_tim_hist((c1):(c2-2),2))/2.d0
        dt(1) = (err_tim_hist(c1+1,2)-err_tim_hist(c1,2))/2.d0  ! Only half first cycle (as the remaining is in the previous cycle
        dt(c2-c1+1) = (err_tim_hist(c2,2)-err_tim_hist(c2-1,2)) ! The full last cycle, because the final point that chould be used is not in this cycle (in total we thus are missing a half time step)
        cycle_time = sum(dt)
        
        call get_step_data_int(stp_ctrl, ctrl, step)
        call get_step_data_dbl(stp_err_scale, err%err_scale, step)
        call get_step_data_dbl(stp_err_scale_ctrl, err%err_scale_ctrl, step)
        
        do ch_ind=1,size(ctrl,1)-1
            if (err_hist_comp(ch_ind)>0) then
                do k2=0,1   ! Load first, then disp. (2*k2-1) = [-1, 1]
                    k3=err_hist_comp(ch_ind)-1+k2   ! First err_hist_comp(ch_ind)-1 (load), then err_hist_comp(ch_ind) (disp)
                    if (stp_ctrl(ch_ind)*(2*k2-1)>0) then    ! Controlled part ctrl and (2*k2-1) have the same sign
                        error_scale_factor = stp_err_scale_ctrl(ch_ind)
                    else    ! The part of the channel not controlled
                        error_scale_factor = stp_err_scale(ch_ind)
                    endif
                    
                    mean_exp_cycle = sum(err_exp_hist(c1:c2, k3)*dt, 1)/cycle_time
                    mean_sim_cycle = sum(err_sim_hist(c1:c2, k3)*dt, 1)/cycle_time
                        
                    err_sim_hist(c1:c2, k3) = error_scale_factor*(dt/cycle_time)* &
                        ( (err_sim_hist(c1:c2, k3) - mean_sim_cycle) - (err_exp_hist(c1:c2, k3) - mean_exp_cycle) )/&
                        (abs(maxval(err_exp_hist(c1:c2, k3)) - minval(err_exp_hist(c1:c2, k3)))+10.d0/huge(1.d0))
                    
                    tmp_error = (mean_sim_cycle-mean_exp_cycle)*error_scale_factor/(maxval(abs(err_exp_hist(c1:c2,k3)))+10.d0/huge(1.d0))
                    error = error + tmp_error**2
                    if (present(evec)) then
                        evec(evec_cnt:(evec_cnt+c2-c1)) = err_sim_hist(c1:c2, k3)
                        evec_cnt = evec_cnt + c2-c1+1
                        evec(evec_cnt) = tmp_error
                        evec_cnt = evec_cnt + 1
                    endif
                    
                enddo ! End of load/disp "loop"
            endif
        enddo ! End of ch_ind loop
        deallocate(dt)
    enddo   ! End of cycle_inds loop
    
    error = (err%err_type_scale**2)*error + ((1.d0-err%err_type_scale)**2)*sum(err_sim_hist**2)
    
end subroutine

subroutine determine_cycles(steps, step_per_cycle, step_initial, cycle_inds)
implicit none
    double precision, intent(in)    :: steps(:)
    double precision, intent(in)    :: step_per_cycle, step_initial
    integer, allocatable            :: cycle_inds(:), cycle_inds_tmp(:)
    
    integer                         :: num_cycles, k1
    double precision                :: end_of_cycle_step
    logical                         :: notlast
    
    
    
    num_cycles = count(steps(2:)>steps(1:size(steps)-1)) ! Calculate number of step changes and use this as upper bound for number of cycles
    allocate(cycle_inds_tmp(num_cycles+1))
    
    notlast = .true.
    k1=1
    cycle_inds_tmp(k1) = findfirst(steps >= step_initial)-1
    do while (notlast)
        k1 = k1 + 1
        
        end_of_cycle_step = steps(cycle_inds_tmp(k1-1)+1) + step_per_cycle
        cycle_inds_tmp(k1) = cycle_inds_tmp(k1-1) + findfirst(steps(cycle_inds_tmp(k1-1)+1:)>end_of_cycle_step)-1
        
        if (cycle_inds_tmp(k1)<cycle_inds_tmp(k1-1)) then
            notlast = .false.
            cycle_inds_tmp(k1) = size(steps)
        endif

    enddo
    
    allocate(cycle_inds(k1))
    cycle_inds = cycle_inds_tmp(1:k1)
    
end subroutine

function findfirst(log_arr) result(ind)
! Replacement for intrinsic function in f2008, but not the same inputs!
implicit none
    logical :: log_arr(:)
    integer :: ind
    integer :: k1
    
    ind = 0
    do k1=1,size(log_arr)
        if (log_arr(k1)) then
            ind = k1
            return
        endif
    enddo
    
end function

! OLD
!subroutine error_settings(err, psize, err_norm_met, error_steps, cyclic_error, cyc_fsim, cyc_fexp, cyc_time, cyc_ptcnt, start_step)
!implicit none
!    type(err_typ)                   :: err
!    integer, intent(in)             :: psize
!    integer                         :: err_norm_met
!    double precision, allocatable   :: error_steps(:)
!    logical                         :: cyclic_error
!    double precision, allocatable   :: cyc_fsim(:,:), cyc_fexp(:,:), cyc_time(:)
!    integer, intent(out)            :: cyc_ptcnt
!    double precision, intent(out)   :: start_step
!    
!    err_norm_met = err%err_norm_met
!    allocate(error_steps,  source=err%err_steps)
!    
!    cyclic_error = err%cyclic_error
!    
!    if (cyclic_error) then
!        allocate(cyc_fsim(err%n_cyc_err, psize))
!        allocate(cyc_fexp(err%n_cyc_err, psize))
!        allocate(cyc_time(err%n_cyc_err+2))
!        cyc_fsim = 0.d0
!        cyc_fexp = 0.d0
!        cyc_time = 0.d0
!        start_step = err%nstep_cyc_initial
!        cyc_ptcnt = 0
!    endif
!    
!
!end subroutine
!    
!subroutine error_update(load_exp, load, disp_exp, disp, ctrl, err_scale, load_scale, disp_scale, error, e_cnt, evec_cnt, evec)
!implicit none
!    double precision, intent(in)    :: load_exp(:), load(:), disp_exp(:), disp(:)
!    integer, intent(in)             :: ctrl(:)
!    double precision, intent(in)    :: err_scale(:), load_scale(:), disp_scale(:)
!    double precision, intent(inout) :: error
!    integer, intent(inout)          :: e_cnt, evec_cnt
!    double precision, allocatable, optional      :: evec(:)
!    
!    double precision, allocatable   :: error_inc(:)
!    integer                         :: evec_cnt_tmp
!    
!    allocate(error_inc(size(ctrl)))
!    
!    error_inc = error_function(load_exp, load, disp_exp, disp, ctrl, err_scale, load_scale, disp_scale)
!    error = error + sum(error_inc**2) ! Add the sum of squares
!    evec_cnt_tmp = evec_cnt + count(err_scale>1.d-20) - 1
!    if (present(evec)) then
!        evec(evec_cnt:evec_cnt_tmp) = pack(error_inc, .not.dbl_comp_array(err_scale,0.d0))
!    endif
!    evec_cnt = evec_cnt_tmp + 1
!    e_cnt = e_cnt + 1
!    
!end subroutine
!
!function error_function(load_exp, load, disp_exp, disp, ctrl, stp_scale, load_scale, disp_scale) result(error)
!implicit none
!    double precision, intent(in)    :: load_exp(:), load(:), disp_exp(:), disp(:)
!    integer, intent(in)             :: ctrl(:)
!    double precision, intent(in)    :: stp_scale(:), load_scale(:), disp_scale(:)
!    double precision, allocatable   :: error(:)
!    ! Give out error as (e_channel1, e_channel2, e_channel3, ...)
!    integer                         :: k1
!        
!    allocate(error(size(ctrl)))
!    error = 0.d0
!    do k1=1,size(ctrl)
!        if (ctrl(k1)>0) then    ! Disp controlled, error on load
!            error(k1) = ( load_exp(k1) - load(k1) )*stp_scale(k1)*load_scale(k1)
!        else                    ! Load controlled, error on disp
!            error(k1) = ( disp_exp(k1) - disp(k1) )*stp_scale(k1)*disp_scale(k1)
!        endif
!    enddo
!end function
!    
!subroutine cyclic_error_update(load_exp, load, disp_exp, disp, ctrl, time, cnt, fsim, fexp, time_vec)
!implicit none
!    double precision, intent(in)                 :: load_exp(:), load(:), disp_exp(:), disp(:)
!    integer, intent(in)                          :: ctrl(:)
!    double precision, intent(in)                 :: time(:)
!    integer, intent(inout)                       :: cnt
!    double precision, allocatable, intent(inout) :: fsim(:,:), fexp(:,:), time_vec(:)
!    ! Internal variables
!    double precision, allocatable                :: tmp1(:), tmp2(:,:)
!    integer                                      :: k1
!    
!    ! If needed, increase size of fsim, fexp and time_vec
!    cnt = cnt + 1
!    if (cnt>size(fsim,1)) then
!        allocate(tmp2(2*size(fsim,1), size(fsim,2)))
!        tmp2(1:size(fsim,1),:) = fsim
!        deallocate(fsim); allocate(fsim, source=tmp2)
!        tmp2(1:size(fexp,1),:) = fexp
!        deallocate(fexp); allocate(fexp, source=tmp2)
!        allocate(tmp1(2*size(fsim,1)+2))
!        tmp1(1:size(time_vec)) = time_vec
!        deallocate(time_vec); allocate(time_vec, source=tmp1)
!    endif
!    
!    ! Add data to saved arrays for later error evaluation (after cycle has completed)
!    do k1=1,size(ctrl)
!        if (ctrl(k1)<0) then    ! Load ctrl - error on disp
!            fsim(cnt,k1) = disp(k1)
!            fexp(cnt,k1) = disp_exp(k1)
!        else                    ! Disp ctrl - error on load
!            fsim(cnt,k1) = load(k1)
!            fexp(cnt,k1) = load_exp(k1)
!        endif    
!    enddo
!    time_vec(cnt+1) = time(2)
!    if (cnt==1) then
!        time_vec(1) = time(2)-time(1)   ! If first in new cycle, set start time of cycle as start time of current step
!    endif
!    
!end subroutine
!                               
!subroutine get_cyclic_error(steps, err, ctrl, load_scale, disp_scale, stp_err_scale, & 
!            last_start, n, cyc_fsim, cyc_fexp, cyc_time, error, e_cnt, evec_cnt, evec)
!implicit none
!    double precision, intent(in)                :: steps(:) ! The current and the next load step, only the current if it is the last step in the experiment.
!    type(err_typ), intent(in)                   :: err
!    integer                                     :: ctrl(:)
!    double precision, intent(in)                :: load_scale(:), disp_scale(:)
!    double precision, intent(in)                :: stp_err_scale(:)
!    double precision, intent(inout)             :: last_start
!    integer, intent(inout)                      :: n ! cyc_ptcnt
!    double precision, intent(inout)             :: cyc_fsim(:,:), cyc_fexp(:,:), cyc_time(:)
!    double precision, intent(inout)             :: error
!    integer, intent(inout)                      :: e_cnt, evec_cnt
!    double precision, allocatable, optional     :: evec(:)
!    
!    integer                                     :: k1, pdim
!    double precision                            :: shape_scale(size(cyc_fsim,2))
!    logical                                     :: dim_use(size(cyc_fsim,2))
!    double precision                            :: cyc_nr_since_last 
!    double precision                            :: fsim_hat, fexp_hat, dt(n), tot_time
!    double precision                            :: mean_error, shape_error
!    double precision                            :: tmp_evec(n)
!    
!    ! Figure out if this is the last step in a cycle
!    cyc_nr_since_last = (steps(size(steps))-last_start)/err%nstep_cyc_err_calc    ! Should normally be 1 if current step is last in cycle and error is evaluated every cycle
!    ! Consider as last step if next step is the first in next cycle, or if this is the last step in the experiment
!    ! The latter could potientially be added as a user setting if the last cycle should be counted to the error measure or not.
!    if (dbl_compare(cyc_nr_since_last, real(int(cyc_nr_since_last+0.5), kind(cyc_nr_since_last)), tol_setting=1.d-6).or.(size(steps)==1)) then
!        ! Pre-calculations
!        pdim = size(cyc_fsim,2) 
!        dim_use = .not.dbl_comp_array(stp_err_scale, 0.d0)
!        where (ctrl>0) ! I.e. true if disp control
!            shape_scale = stp_err_scale*load_scale
!        elsewhere
!            shape_scale = stp_err_scale*disp_scale
!        end where
!        
!        ! Errors should be calculated and we reset all values
!        cyc_time(n+2) = cyc_time(n+1)
!        ! Calculate custom scale for shape (using entire cycle and not just the current step)
!        if (err%err_norm_met==1) then
!            do k1=1,pdim
!                if (dim_use(k1)) then
!                    shape_scale(k1) = stp_err_scale(k1)/(maxval(cyc_fexp(1:n,k1))-minval(cyc_fexp(1:n,k1)))
!                endif
!            enddo
!        endif
!        
!        dt = (cyc_time(3:(n+2))-cyc_time(1:n))/2
!        tot_time = cyc_time(n)+dt(n)/2 - (cyc_time(2) - dt(1)/2)
!
!        
!        ! Calculate errors
!        do k1=1,pdim
!            if (dim_use(k1)) then
!                ! Calculate mean values for the cycle
!                fsim_hat = sum(cyc_fsim(1:n,k1)*dt)/tot_time
!                fexp_hat = sum(cyc_fexp(1:n,k1)*dt)/tot_time
!                ! Calculate the mean error for the cycle
!                mean_error = stp_err_scale(k1)*(fsim_hat - fexp_hat)/maxval(abs(cyc_fexp(1:n,k1)))    ! Mean error scaled by maximum absolute value
!                ! Calculate the shape errors for the cycle
!                tmp_evec = (cyc_fsim(1:n,k1) - fsim_hat - (cyc_fexp(1:n,k1) - fexp_hat))*shape_scale(k1)
!                shape_error = sum( ((tmp_evec)**2)*dt )/tot_time
!                if (present(evec)) then
!                    evec(evec_cnt) = err%err_type_scale*mean_error
!                    evec((evec_cnt+1):(evec_cnt+n)) = (1.d0-err%err_type_scale)*tmp_evec*dt/tot_time
!                endif
!                evec_cnt = evec_cnt + n + 1
!                error = error + err%err_type_scale*mean_error**2 + (1.d0-err%err_type_scale)*shape_error
!            endif
!        enddo
!        e_cnt = e_cnt + 1
!        ! Reset values
!        n = 0
!    endif
!    
!
!end subroutine
!
!                               
end module
