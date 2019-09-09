! Example routine for user error
! Detailed interface description
! settings
!   User settings as defined by the *user_settings input (double vector)
!
! ctrl(:,:)
!   Channel control, each column = (stepnr, channel1, channel2, ...) (As specified in input, but transpose)
!
! err_tim_hist(:,:)
!   First column contains the steps, 2nd column contains the time
!
! err_exp_hist(:,:)
!   Each column contains the experimental data for each channel (both control modes). Column number given by err_hist_comp.
!
! err_sim_hist(:,:)
!   On input each column contains the simulated data for each channel (both control modes). Column number given by err_hist_comp.
!   This is intent inout to allow it to be used to changed to error instead inside uerror. What it is on output doesn't matter.
!
! err_hist_comp(:)
!   Describes where in err_[sim/exp]_hist the disp values are put. Load values are put at this position - 1. 
!   
! error
!   Scalar error calculated by uerror. This error is returned from the present simulation and linearly added to the other errors. 
!
! evec(:)
!   If (present(evec)) then return a (signed) error vector. This is used to e.g. calculate the correlation matrix. 
!   These errors are appended to the error vector for all simulations.
!
! End of detailed interface description
!
! In this example the error is defined as the square difference in the relative maximum values for the non-controlled variables
! settings(1) gives the scale factor for the errors
! The first row of ctrl is used
    
subroutine uerror(settings, ctrl, err_tim_hist, err_exp_hist, err_sim_hist, err_hist_comp, error, evec)
!DEC$ ATTRIBUTES DLLEXPORT :: uerror
implicit none
    double precision, intent(in)            :: settings(:)              ! Settings given by err%user_settings
    double precision, intent(in)            :: ctrl(:,:)                ! Control information
    double precision, intent(inout)         :: err_tim_hist(:,:), err_exp_hist(:,:), err_sim_hist(:,:)
    integer                                 :: err_hist_comp(:)         ! Describes where in err_[sim/exp]_hist the disp values are put (load values are put at pos-1)
    double precision, intent(out)           :: error                    ! Calculated error
    double precision, allocatable, optional :: evec(:)                  ! Error vector
    
    integer                                 :: num_errors, k1, k2, nr
    double precision, allocatable           :: error_vector(:)
    double precision                        :: max_sim, max_exp
    
    num_errors = count(err_hist_comp>0) ! This is the number of channels the error is calculated for (that have scale factor > 0)
    allocate(error_vector(num_errors))
    
    k2 = 0
    do k1=1,size(err_hist_comp) ! Loop over all components and check which ones that have data
        nr = err_hist_comp(k1)
        if (nr>0) then  ! Data saved for this component, load at nr-1 and disp at nr
            k2 = k2 + 1
            if (ctrl(k1+1,1)>0) then ! Displacement control, error on load
                max_sim = maxval(err_sim_hist(:, nr-1))
                max_exp = maxval(err_exp_hist(:, nr-1))
            else                    ! Load control, error on disp
                max_sim = maxval(err_sim_hist(:, nr))
                max_exp = maxval(err_exp_hist(:, nr))
            endif
            error_vector(k2) = settings(1)*(max_exp - max_sim)/max(max_exp, 10.d0/huge(1.d0))   ! To avoid 0/0 if max_exp=0
        endif
    enddo
    
    error = sum(error_vector**2)/num_errors
    
    if (present(evec)) then
        if (.not.allocated(evec)) then
            allocate(evec, source=error_vector)
        else
            evec = error_vector
        endif
    endif
    
end subroutine
