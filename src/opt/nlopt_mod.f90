    
module nlopt_mod
    use types_mod
    use output_mod
    use simulation_mod
    use convert_mpar_mod
    use nlopt_constants_mod
    use gen_util_mod
    implicit none
    
    private
    public nlopt_optimization

    ! Use this to ensure compatibility with the nlopt library, if matmodfit is compiled using 64 bit integers (i.e. -fdefault-integer-8)
    ! This may occur if the mkl library guide for linking is used, with the assumption that 64 bit integers go together with 64 bit system...
    ! The default project settings should use kind=4 by default, however.
    integer, parameter :: nlopt_ityp = 4
    
    contains

subroutine nlopt_optimization(error_out, xvar_out, error_in, xvar_in, opt, f_data, opt_step)                   
! Use nlopt: https://nlopt.readthedocs.io/en/latest/ to perform the optimization
! The subroutine optimize sets up the optimization using this optimization library
implicit none
    double precision                :: error_in(:)
    double precision                :: xvar_in(:,:)
    double precision, allocatable   :: error_out(:)
    double precision, allocatable   :: xvar_out(:,:)
    type(fdata_typ)                 :: f_data
    type(opt_typ)     		        :: opt
    
    integer                         :: nvar
    integer                         :: nset
    integer, allocatable            :: set_pos(:)
    integer(kind=nlopt_ityp)        :: ires
    integer 		                :: k1, opt_step
    
    !Pointer type of variable, should be large enough to have address on both 32 and 64 bit system:
    integer(kind=8)         :: nl_opt 	!See nlopt fortran example
    external                :: nlopt_obj_fun
    integer                 :: clock_count, clock_rate
    
    nvar = size(xvar_in,1)
    nset = min(size(error_in), opt%start%num_sets)
    allocate(set_pos(nset))
    set_pos = get_n_smallest_inds(error_in, nset)
!   allocate(xvar_out, source=xvar_in(:,set_pos)) Seems to give issues with gfortran?
    allocate(xvar_out(nvar,nset))
    xvar_out = xvar_in(:,set_pos)
    allocate(error_out(nset))
    error_out = error_in(set_pos)
    

    
    
    ! == Mandatory setup ==
    ! Optimization object
    call nlo_create(nl_opt, int(opt%start%algorithm, kind=nlopt_ityp), int(nvar, kind=nlopt_ityp))
        
    ! Objective function
    call nlo_set_min_objective(ires, nl_opt, nlopt_obj_fun, f_data) ! f_data is passed as reference
    if (ires<0) then
        call write_output('Could not set objective function for optimization step '//int2str(opt_step), 'error', 'opt')
    endif
        
    ! Lower bounds
    call nlo_set_lower_bounds(ires, nl_opt, get_lower_bounds(f_data%glob))
    if (ires<0) then
        call write_output('Could not set lower bounds in optimization step '//int2str(opt_step), 'error', 'opt')
    endif
        
    ! Upper bounds
    call nlo_set_upper_bounds(ires, nl_opt, get_upper_bounds(f_data%glob))
    if (ires<0) then
        call write_output('Could not set upper bounds in optimization step '//int2str(opt_step), 'error', 'opt')
    endif
    
    ! == Non-mandatory options ==
    ! User option setup (more options, but non-mandatory)
    call nlopt_set_user_pref(nl_opt, opt, opt_step, nvar)
    
    ! == Run optimization for each set == 
    do k1=1,nset
        ! Prepare optimization status input
        f_data%optdata%min_objfun = error_out(k1)   ! Set initial smallest value as "infinity"
        f_data%optdata%n_objfun_calls = 0           ! Set number of objective function calls to zero
        call system_clock ( count=clock_count, count_rate=clock_rate)
        f_data%optdata%start_time = real(clock_count)/real(clock_rate)
        call nlo_optimize(ires, nl_opt, xvar_out(:,k1), error_out(k1))  ! Run optimization
        call nlopt_termination(ires, k1)                        ! Write out termination message
    enddo
    
    ! Delete optimization object from memory
    call nlo_destroy(nl_opt)
    
end subroutine nlopt_optimization

subroutine nlopt_set_user_pref(nl_opt, opt, k1, nvar)
implicit none
    integer(kind=8)          :: nl_opt   ! "Pointer" to nlopt optimization object
    type(opt_typ)       :: opt      ! Optimization settings for stage k1
    integer                  :: k1       ! Optimization stage number (used here only for error output messages)
    integer                  :: nvar     ! Number of optimized parameters, used for checking that input is ok
    integer(kind=nlopt_ityp) :: ires     ! Result from setting nlo options (error checking not implemented yet)
    double precision         :: tmp(nvar)! Temp double for putting all variables to single value 
    
        if (any(.not.dbl_comp_array(opt%start%initial_step,0.d0,1.d-20))) then
            if (size(opt%start%initial_step)==nvar) then
                call nlo_set_initial_step(ires, nl_opt, opt%start%initial_step)
            elseif (size(opt%start%initial_step)==1) then
                call nlo_set_initial_step1(ires, nl_opt, opt%start%initial_step(1))
            else
                call write_output('Wrong number of initial step parameters in optimization step '//int2str(k1), 'error', 'opt')
            endif
            if (ires<0) then
                call write_output('Could not set initial_step in optimization step '//int2str(k1), 'error', 'opt')
            endif
            
        endif
        
        if (opt%end_cond%stopval>0.d0) then
            call nlo_set_stopval(ires, nl_opt, opt%end_cond%stopval)
            if (ires<0) then
                call write_output('Could not set stop value in optimization step '//int2str(k1), 'error', 'opt')
            endif
        endif
        
        if (opt%end_cond%ftol_rel>0.d0) then
            call nlo_set_ftol_rel(ires, nl_opt, opt%end_cond%ftol_rel)
            if (ires<0) then
                call write_output('Could not set relative objective tolerance in optimization step '//int2str(k1), 'error', 'opt')
            endif
        endif
        
        if (opt%end_cond%ftol_abs>0.d0) then
            call nlo_set_ftol_rel(ires, nl_opt, opt%end_cond%ftol_abs)
            if (ires<0) then
                call write_output('Could not set absolute objective tolerance in optimization step '//int2str(k1), 'error', 'opt')
            endif
        endif
        
        if (any(.not.dbl_comp_array(opt%end_cond%xtol_rel,0.d0,1.d-20))) then
            if (size(opt%end_cond%xtol_rel)==nvar) then
                call nlo_set_xtol_rel(ires, nl_opt, opt%end_cond%xtol_rel)
            elseif (size(opt%end_cond%xtol_rel)==1) then
                tmp = opt%end_cond%xtol_rel(1)
                call nlo_set_xtol_rel(ires, nl_opt, tmp)
            else
                call write_output('Wrong number of xtol_rel parameters in optimization step '//int2str(k1), 'error', 'opt')
            endif
            if (ires<0) then
                call write_output('Could not set relative variable tolerance in optimization step '//int2str(k1), 'error', 'opt')
            endif
        endif
        
        if (any(.not.dbl_comp_array(opt%end_cond%xtol_abs,0.d0,1.d-20))) then
            if (size(opt%end_cond%xtol_abs)==nvar) then
                call nlo_set_xtol_abs(ires, nl_opt, opt%end_cond%xtol_abs)
            elseif (size(opt%end_cond%xtol_abs)==1) then
                call nlo_set_xtol_abs1(ires, nl_opt, opt%end_cond%xtol_abs(1))
            else
                call write_output('Wrong number of xtol_abs parameters in optimization step '//int2str(k1), 'error', 'opt')
            endif
            if (ires<0) then
                call write_output('Could not set absolute variable tolerance in optimization step '//int2str(k1), 'error', 'opt')
                stop
            endif
        endif
        
        if (opt%end_cond%maxeval>0) then
            call nlo_set_maxeval(ires, nl_opt, int(opt%end_cond%maxeval, kind=nlopt_ityp))
            if (ires<0) then
                call write_output('Could not set max number of simulations in optimization step '//int2str(k1), 'error', 'opt')
            endif
        endif
        
        if (opt%end_cond%maxtime>0.d0) then
            call nlo_set_maxtime(ires, nl_opt, opt%end_cond%maxtime)
            if (ires<0) then
                call write_output('Could not set max optimization time in optimization step '//int2str(k1), 'error', 'opt')
            endif
        endif
        
        ! This does not work, for some reason nlosr ans nlosrt is not found in nlopt, even though it is exported by the dll file...
        if (opt%start%use_fixed_seed) then
            ! call nlosr(opt%start%fixed_seed)    ! Set fixed seed
        else
            ! call nlosrt                         ! Set auto-seed based on cpu-time
        endif
        
        
    
end subroutine
    
subroutine nlopt_termination(ires, opt_stage)
implicit none
    integer(kind=nlopt_ityp) :: ires
    integer                  :: opt_stage
    
    call write_output('Optimization of parameter set '//int2str(opt_stage)//' ended (output code '//int2str(ires)//')', 'status', 'opt')
    select case (ires)
    case (nlopt_success)
        call write_output('Generic success of optimization', 'status', 'opt')
    case (nlopt_stopval_reached) 
        call write_output('Optimization ended because sufficiently good fit was obtained', 'status', 'opt')
    case (nlopt_ftol_reached)
        call write_output('Change of objective function less than specified tolerance', 'status', 'opt')
    case (nlopt_xtol_reached) 
        call write_output('Change of all parameters are less than specified tolerance', 'status', 'opt')
    case (nlopt_maxeval_reached)
        call write_output('The maximum number of function evaluations reached', 'status', 'opt')
    case (nlopt_maxtime_reached)
        call write_output('The maximum amount of time allowed was reached', 'status', 'opt')
    case (nlopt_failure)
        call write_output('Generic failure of nlopt', 'status', 'opt')
    case (nlopt_invalid_args)
        call write_output('Invalid arguments to nlopt', 'status', 'opt')
        call write_output('Check if lower bounds are bigger than upper bounds', 'status', 'opt')
    case (nlopt_out_of_memory)
        call write_output('Ran out of memory', 'status', 'opt')
    case (nlopt_roundoff_limited)
        call write_output('Halted because roundoff errors limited progress (Result usually still useful)', 'status', 'opt')
    case (nlopt_forced_stop)
        call write_output('Forced termination of nlopt from inside simulation', 'status', 'opt')
    case default
        call write_output('Unknown error code', 'status', 'opt')
    end select
    
end subroutine

! Should probably be moved to a utility module

end module nlopt_mod
    
