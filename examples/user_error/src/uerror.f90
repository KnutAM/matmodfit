module uerror_mod
    implicit none
    
    contains
    
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
    
    subroutine get_step_data_int(stp_var, glob_var, stpnr)
    implicit none
        integer, intent(inout)          :: stp_var(:)
        double precision, intent(in)    :: glob_var(:,:)
        double precision, intent(in)    :: stpnr
        integer                         :: k1

        do k1 = size(glob_var, 2),1,-1
            if (stpnr>(glob_var(1,k1)-1.d-10)) then
                stp_var = nint(glob_var(2:size(glob_var, 1),k1))
                exit
            endif
        enddo
    
        if (stpnr<(glob_var(1,1)-1.d-10)) then
            write(*,*) 'First entry of step data types must be equal to first step of experiment (or 1 if no step column)'
        endif
    
    end subroutine
    
end module
    
    
    
    
subroutine uerror(settings, ctrl, e_cnt, err_tim_hist, err_exp_hist, err_sim_hist, err_hist_comp, error, evec)
!DEC$ ATTRIBUTES DLLEXPORT :: uerror
use uerror_mod
implicit none
    double precision, intent(in)            :: settings(:)              ! Settings given by err%user_settings
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
    integer                                 :: err_norm_met
    double precision                        :: sfac
    
    err_norm_met = int(settings(1))
    num_channels = size(ctrl,1)-1    ! Number of channels
    allocate(stp_ctrl(num_channels), stp_err_scale(num_channels), stp_err_scale_ctrl(num_channels))
    
    ! Calculate global scaling factors if scaling with max-min in experiment
    if (err_norm_met==2) then 
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
        if (size(settings)>=num_channels+1) then
            stp_err_scale = settings(2:(num_channels+1))
        else
            stp_err_scale = 1.d0
        endif
        stp_err_scale_ctrl = 0.d0       
        
        do ch_ind=1,num_channels
            if (err_hist_comp(ch_ind)>0) then
                do k2=0,1   ! Load first, then disp. (2*k2-1) = [-1, 1]
                    k1=err_hist_comp(ch_ind)-1+k2   ! First err_hist_comp(ch_ind)-1 (load), then err_hist_comp(ch_ind) (disp)
                    if (stp_ctrl(ch_ind)*(2*k2-1)>0) then    ! Controlled part ctrl and (2*k2-1) have the same sign
                        error_scale_factor = stp_err_scale_ctrl(ch_ind)
                    else    ! The part of the channel not controlled
                        error_scale_factor = stp_err_scale(ch_ind)
                    endif
                        
                    select case (err_norm_met)
                    case (0) ! No scaling
                        err_sim_hist(thisind:(nextind-1), k1) = error_scale_factor*(err_sim_hist(thisind:(nextind-1), k1) - err_exp_hist(thisind:(nextind-1), k1))
                    case (1) ! Scale with max - min in step
                        sfac = maxval(err_exp_hist(thisind:(nextind-1), k1)) - minval(err_exp_hist(thisind:(nextind-1), k1))
                        if (sfac>10.d0/huge(1.d0)) then
                            err_sim_hist(thisind:(nextind-1), k1) = error_scale_factor*(err_sim_hist(thisind:(nextind-1), k1) - err_exp_hist(thisind:(nextind-1), k1))/sfac
                        else
                            err_sim_hist(thisind:(nextind-1), k1) = 0.d0
                        endif
                    case (2) ! Scale with max - min in experiment
                        err_sim_hist(thisind:(nextind-1), k1) = error_scale_factor*(err_sim_hist(thisind:(nextind-1), k1) - err_exp_hist(thisind:(nextind-1), k1))*scaling(k1)
                    case default
                        write(*,*) 'err_norm_met=', err_norm_met, 'is not supported'
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
