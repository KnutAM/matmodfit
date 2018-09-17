module sim_error_mod
use gen_util_mod
use types_mod
use output_mod
    implicit none
    
    private
    
    public error_settings       ! Read in error settings and allocate error variables
    public error_update         ! Standard error updating
    public cyclic_error_update  ! Special error update for cyclic error, saves data into cyc_fsim, cyc_fexp and cyc_time
    public get_cyclic_error     ! Routine for calculating the error measures for one cycle
    
    contains

subroutine error_settings(err, psize, err_norm_met, error_steps, cyclic_error, cyc_fsim, cyc_fexp, cyc_time, cyc_ptcnt, start_step)
implicit none
    type(err_typ)                   :: err
    integer, intent(in)             :: psize
    integer                         :: err_norm_met
    double precision, allocatable   :: error_steps(:)
    logical                         :: cyclic_error
    double precision, allocatable   :: cyc_fsim(:,:), cyc_fexp(:,:), cyc_time(:)
    integer, intent(out)            :: cyc_ptcnt
    double precision, intent(out)   :: start_step
    
    err_norm_met = err%err_norm_met
    allocate(error_steps,  source=err%err_steps)
    
    cyclic_error = err%cyclic_error
    
    if (cyclic_error) then
        allocate(cyc_fsim(err%n_cyc_err, psize))
        allocate(cyc_fexp(err%n_cyc_err, psize))
        allocate(cyc_time(err%n_cyc_err+2))
        cyc_fsim = 0.d0
        cyc_fexp = 0.d0
        cyc_time = 0.d0
        start_step = err%nstep_cyc_initial
        cyc_ptcnt = 0
    endif
    

end subroutine
    
subroutine error_update(load_exp, load, disp_exp, disp, ctrl, err_scale, load_scale, disp_scale, error, e_cnt, evec_cnt, evec)
implicit none
    double precision, intent(in)    :: load_exp(:), load(:), disp_exp(:), disp(:)
    integer, intent(in)             :: ctrl(:)
    double precision, intent(in)    :: err_scale(:), load_scale(:), disp_scale(:)
    double precision, intent(inout) :: error
    integer, intent(inout)          :: e_cnt, evec_cnt
    double precision, allocatable, optional      :: evec(:)
    
    double precision, allocatable   :: error_inc(:)
    integer                         :: evec_cnt_tmp
    
    allocate(error_inc(size(ctrl)))
    
    error_inc = error_function(load_exp, load, disp_exp, disp, ctrl, err_scale, load_scale, disp_scale)
    error = error + sum(error_inc**2) ! Add the sum of squares
    evec_cnt_tmp = evec_cnt + count(err_scale>1.d-20) - 1
    if (present(evec)) then
        evec(evec_cnt:evec_cnt_tmp) = pack(error_inc, .not.dbl_comp_array(err_scale,0.d0))
    endif
    evec_cnt = evec_cnt_tmp + 1
    e_cnt = e_cnt + 1
    
end subroutine
                        
function error_function(load_exp, load, disp_exp, disp, ctrl, stp_scale, load_scale, disp_scale) result(error)
implicit none
    double precision, intent(in)    :: load_exp(:), load(:), disp_exp(:), disp(:)
    integer, intent(in)             :: ctrl(:)
    double precision, intent(in)    :: stp_scale(:), load_scale(:), disp_scale(:)
    double precision, allocatable   :: error(:)
    ! Give out error as (e_channel1, e_channel2, e_channel3, ...)
    integer                         :: k1
        
    allocate(error(size(ctrl)))
    error = 0.d0
    do k1=1,size(ctrl)
        if (ctrl(k1)>0) then    ! Disp controlled, error on load
            error(k1) = ( load_exp(k1) - load(k1) )*stp_scale(k1)*load_scale(k1)
        else                    ! Load controlled, error on disp
            error(k1) = ( disp_exp(k1) - disp(k1) )*stp_scale(k1)*disp_scale(k1)
        endif
    enddo
end function
    
subroutine cyclic_error_update(load_exp, load, disp_exp, disp, ctrl, time, cnt, fsim, fexp, time_vec)
implicit none
    double precision, intent(in)                 :: load_exp(:), load(:), disp_exp(:), disp(:)
    integer, intent(in)                          :: ctrl(:)
    double precision, intent(in)                 :: time(:)
    integer, intent(inout)                       :: cnt
    double precision, allocatable, intent(inout) :: fsim(:,:), fexp(:,:), time_vec(:)
    ! Internal variables
    double precision, allocatable                :: tmp1(:), tmp2(:,:)
    integer                                      :: k1
    
    ! If needed, increase size of fsim, fexp and time_vec
    cnt = cnt + 1
    if (cnt>size(fsim,1)) then
        allocate(tmp2(2*size(fsim,1), size(fsim,2)))
        tmp2(1:size(fsim,1),:) = fsim
        deallocate(fsim); allocate(fsim, source=tmp2)
        tmp2(1:size(fexp,1),:) = fexp
        deallocate(fexp); allocate(fexp, source=tmp2)
        allocate(tmp1(2*size(fsim,1)+2))
        tmp1(1:size(time_vec)) = time_vec
        deallocate(time_vec); allocate(time_vec, source=tmp1)
    endif
    
    ! Add data to saved arrays for later error evaluation (after cycle has completed)
    do k1=1,size(ctrl)
        if (ctrl(k1)<0) then    ! Load ctrl - error on disp
            fsim(cnt,k1) = disp(k1)
            fexp(cnt,k1) = disp_exp(k1)
        else                    ! Disp ctrl - error on load
            fsim(cnt,k1) = load(k1)
            fexp(cnt,k1) = load_exp(k1)
        endif    
    enddo
    time_vec(cnt+1) = time(2)
    if (cnt==1) then
        time_vec(1) = time(2)-time(1)   ! If first in new cycle, set start time of cycle as start time of current step
    endif
    
end subroutine
                               
subroutine get_cyclic_error(steps, err, ctrl, load_scale, disp_scale, stp_err_scale, & 
            last_start, n, cyc_fsim, cyc_fexp, cyc_time, error, e_cnt, evec_cnt, evec)
implicit none
    double precision, intent(in)                :: steps(:) ! The current and the next load step, only the current if it is the last step in the experiment.
    type(err_typ), intent(in)                   :: err
    integer                                     :: ctrl(:)
    double precision, intent(in)                :: load_scale(:), disp_scale(:)
    double precision, intent(in)                :: stp_err_scale(:)
    double precision, intent(inout)             :: last_start
    integer, intent(inout)                      :: n ! cyc_ptcnt
    double precision, intent(inout)             :: cyc_fsim(:,:), cyc_fexp(:,:), cyc_time(:)
    double precision, intent(inout)             :: error
    integer, intent(inout)                      :: e_cnt, evec_cnt
    double precision, allocatable, optional     :: evec(:)
    
    integer                                     :: k1, pdim
    double precision                            :: shape_scale(size(cyc_fsim,2))
    logical                                     :: dim_use(size(cyc_fsim,2))
    double precision                            :: cyc_nr_since_last 
    double precision                            :: fsim_hat, fexp_hat, dt(n), tot_time
    double precision                            :: mean_error, shape_error
    double precision                            :: tmp_evec(n)
    
    ! Figure out if this is the last step in a cycle
    cyc_nr_since_last = (steps(size(steps))-last_start)/err%nstep_cyc_err_calc    ! Should normally be 1 if current step is last in cycle and error is evaluated every cycle
    ! Consider as last step if next step is the first in next cycle, or if this is the last step in the experiment
    ! The latter could potientially be added as a user setting if the last cycle should be counted to the error measure or not.
    if (dbl_compare(cyc_nr_since_last, real(int(cyc_nr_since_last+0.5), kind(cyc_nr_since_last)), tol_setting=1.d-6).or.(size(steps)==1)) then
        ! Pre-calculations
        pdim = size(cyc_fsim,2) 
        dim_use = .not.dbl_comp_array(stp_err_scale, 0.d0)
        where (ctrl>0) ! I.e. true if disp control
            shape_scale = stp_err_scale*load_scale
        elsewhere
            shape_scale = stp_err_scale*disp_scale
        end where
        
        ! Errors should be calculated and we reset all values
        cyc_time(n+2) = cyc_time(n+1)
        ! Calculate custom scale for shape (using entire cycle and not just the current step)
        if (err%err_norm_met==1) then
            do k1=1,pdim
                if (dim_use(k1)) then
                    shape_scale(k1) = stp_err_scale(k1)/(maxval(cyc_fexp(1:n,k1))-minval(cyc_fexp(1:n,k1)))
                endif
            enddo
        endif
        
        dt = (cyc_time(3:(n+2))-cyc_time(1:n))/2
        tot_time = cyc_time(n)+dt(n)/2 - (cyc_time(2) - dt(1)/2)

        
        ! Calculate errors
        do k1=1,pdim
            if (dim_use(k1)) then
                ! Calculate mean values for the cycle
                fsim_hat = sum(cyc_fsim(1:n,k1)*dt)/tot_time
                fexp_hat = sum(cyc_fexp(1:n,k1)*dt)/tot_time
                ! Calculate the mean error for the cycle
                mean_error = stp_err_scale(k1)*(fsim_hat - fexp_hat)/maxval(abs(cyc_fexp(1:n,k1)))    ! Mean error scaled by maximum absolute value
                ! Calculate the shape errors for the cycle
                tmp_evec = (cyc_fsim(1:n,k1) - fsim_hat - (cyc_fexp(1:n,k1) - fexp_hat))*shape_scale(k1)
                shape_error = sum( ((tmp_evec)**2)*dt )/tot_time
                if (present(evec)) then
                    evec(evec_cnt) = err%err_type_scale*mean_error
                    evec((evec_cnt+1):(evec_cnt+n)) = (1.d0-err%err_type_scale)*tmp_evec*dt/tot_time
                endif
                evec_cnt = evec_cnt + n + 1
                error = error + err%err_type_scale*mean_error**2 + (1.d0-err%err_type_scale)*shape_error
            endif
        enddo
        e_cnt = e_cnt + 1
        ! Reset values
        n = 0
    endif
    

end subroutine

                               
end module
