! subroutine that matches the required format by nlopt package
! Since this is called via the external keyword it cannot be part of a module
subroutine nlopt_obj_fun(error, nvar, xvar, grad, need_grad, f_data)
    use types_mod
    use simulation_mod
    use output_mod
    use convert_mpar_mod
    implicit none
    integer(kind=4)     :: nvar, need_grad
    double precision    :: error, xvar(nvar), grad(nvar), time
    type(fdata_typ)        :: f_data
    integer             :: clock_count, clock_rate
    double precision, allocatable   :: ipar(:)
    
    if (need_grad/=0) then
        call simulate_gradient(error, grad, xvar, f_data)
    else
        call simulate(error, xvar, f_data)
    endif
    f_data%optdata%n_objfun_calls = f_data%optdata%n_objfun_calls + 1 !Increment number of calls counter

    
	if (f_data%glob%error_history) then
        allocate(ipar(size(f_data%glob%ipar_init)))
        call xvar_to_ipar(ipar, xvar, f_data%glob)
		call write_errhist(f_data%optdata%n_objfun_calls, error, ipar)
	endif
		
    ! Write status update if improved objective function
    if (error<f_data%optdata%min_objfun) then
        ! Objective function value
        call write_output('Objective function value = '//dbl2str(error), 'status', 'opt')
        call write_output('Improvement in obj value = '//dbl2str(f_data%optdata%min_objfun - error), 'status', 'opt', loc=.false.)
        f_data%optdata%min_objfun = error   ! Update the saved best objective function value
        
        ! Number of simulations
        call write_output('Number of simulation runs: '//int2str(f_data%optdata%n_objfun_calls), 'status', 'opt', loc=.false.)
        
        ! Time spent on current optimization stage
        call system_clock ( clock_count, clock_rate)
        time = real(clock_count)/real(clock_rate)
        time = time - f_data%optdata%start_time
        if (time<0.d0) then !Ensure that time counter doesn't overflow... (implies negative time)
            time = time + f_data%optdata%start_time
            call system_clock(count_max=clock_count)
            f_data%optdata%start_time = f_data%optdata%start_time - real(clock_count)/real(clock_rate)
            time = time - f_data%optdata%start_time
        endif
        
        call write_output('Optimization stage time usage: '//tim2str(time), 'status', 'opt', loc=.false.)
        
        ! Write a blank line just to make the message separation more clear
        call write_output(' ', 'status', 'opt', loc=.false.)
    endif
    
end subroutine nlopt_obj_fun
    
