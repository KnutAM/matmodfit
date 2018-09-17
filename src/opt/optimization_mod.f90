    
module optimization_mod
    use types_mod
    use output_mod
    use simulation_mod
    use nlopt_mod
    use convert_mpar_mod
    use gen_util_mod
    use usr_interface_mod
    use iso_c_binding
    implicit none
    
    private
    public optimize
    
    contains
    
! main optimization subroutine
subroutine optimize(errors, xvars, f_data, opt)
implicit none
    double precision, allocatable   :: errors(:)    ! Error for best set of material parameters
    double precision, allocatable   :: xvars(:,:)   ! Best set of scaled and optimized material parameters
    type(fdata_typ)                 :: f_data
    type(opt_typ)     		        :: opt(:)
    
    double precision, allocatable   :: error_out(:)
    double precision, allocatable   :: xvar_out(:,:)
    integer                         :: otype
    integer                         :: k1
    
    ! Each optimization function (at stage k) takes a number of sets as input. 
    ! Its opt settings (opt(k)) determines the maximum number to be used. 
    ! It returns the set of optimized parameters and their corresponding errors, 
    ! these lists doesn't need to be sorted. Each optimization function should do this 
    ! if desired in the beginning when selecting which sets to work on.
    
    allocate(errors(size(xvars,2)))             ! The error for each variable set
    errors = huge(1.d0)
    
    do k1=1,size(opt)
        otype = opt(k1)%otype
        if (otype==-1) then
            call usr_optimization(error_out, xvar_out, errors, xvars, opt(k1), f_data)
        elseif (otype==1) then
            call nlopt_optimization(error_out, xvar_out, errors, xvars, opt(k1), f_data, k1)
        elseif (otype==2) then
            call pspace_optimization(error_out, xvar_out, opt(k1), f_data)
        else
            call write_output('otype = '//int2str(otype)//' is not supported', 'opt', 'error')
        endif
        
        ! Should write some results file after each optimization stage
        call optstage_output(xvar_out, error_out, f_data%glob, k1)
        
        ! Reallocate in and out (this could have been improved performance wise, 
        deallocate(errors, xvars)             !but no point as it is only called 
        allocate(errors, source=error_out)    !a few times for each full run...
        allocate(xvars, source=xvar_out)
        deallocate(error_out, xvar_out)
    enddo
    
end subroutine optimize

subroutine optstage_output(xvar, error, glob, optnr)
implicit none
    double precision                :: xvar(:,:)    ! Resulting scaled optimized parameters
    double precision                :: error(:)     ! Error corresponding to xvar
    type(glob_typ)                  :: glob         ! Global settings (used to scale back to regular mpar)
    integer                         :: optnr        ! Optnumber (good to write out)
    integer                         :: k1
    double precision, allocatable   :: mpar(:)
    
    ! For now, we write to the log file. This should be changed to a more appropriate file later
    call write_output('Results from optimization stage nr. '//int2str(optnr)//':', 'status', 'opt')
    call write_output('set nr, error, material parameters ...', loc=.false.)
    allocate(mpar, source=glob%ipar_init)
    do k1=1,size(error)
        call xvar_to_mpar(mpar, xvar(:,k1), glob)
        call write_optstage_output(k1, error(k1), mpar)
    enddo    
    
end subroutine

  
subroutine pspace_optimization(error_out, xvar_out, opt, f_data)
implicit none
    double precision, allocatable   :: error_out(:)
    double precision, allocatable   :: xvar_out(:,:)
    type(fdata_typ)                 :: f_data
    type(opt_typ)     		        :: opt
    integer                         :: nset
    integer                         :: k1
    double precision, allocatable   :: ipar(:)
    
    
    nset = opt%start%num_sets   ! May be updated during call to distribute_xvar
    
    call distribute_xvar(xvar_out, get_lower_bounds(f_data%glob), get_upper_bounds(f_data%glob), &
                         nset, opt%start%algorithm)
    
    allocate(error_out(nset))
    allocate(ipar(size(f_data%glob%ipar_init)))
    
    ! Ideal loop for parallelization?
    do k1=1,nset
        call simulate(error_out(k1), xvar_out(:,k1), f_data)
        call write_output('parameter set '//int2str(k1)//'/'//int2str(nset)//' completed', 'status', 'opt')
        if (f_data%glob%error_history) then
            call xvar_to_ipar(ipar, xvar_out(:,k1), f_data%glob)
		    call write_errhist(k1, error_out(k1), ipar)
	    endif
    enddo
    
    
end subroutine pspace_optimization

subroutine distribute_xvar(xvar, lb, ub, nset, algorithm)
implicit none
    double precision, allocatable   :: xvar(:,:)    ! Variables to be distributed
    double precision                :: lb(:), ub(:) ! Lower and upper bounds for variables
    integer                         :: nset         ! Number of variable sets
    integer                         :: algorithm    ! Distribution algorithm
    
    if (algorithm==1) then
        call full_factorial(xvar, lb, ub, nset)
    elseif (algorithm==2) then
        call lhs_basic(xvar, lb, ub, nset)
    else
        call write_output('algorithm '//int2str(algorithm)//' is not supported for pspace optimization', 'error', 'opt')
    endif
    
end subroutine

subroutine full_factorial(xvar, lb, ub, nset)
implicit none
    double precision, allocatable   :: xvar(:,:)    ! Variables to be distributed
    double precision                :: lb(:), ub(:) ! Lower and upper bounds for variables
    integer                         :: nset         ! Number of variable sets
    integer                         :: nvar         ! Number of variables in each set
    integer                         :: nlev         ! Number of variable levels
    integer                         :: k1, k2       ! Iterators
    integer, allocatable            :: xpos(:,:)    ! Choice of level for each material parameter
    
    nvar = size(lb)
    nlev = int(real(nset)**(1.d0/real(nvar)))       ! Calculate the number of levels of each parameter, rounded down
    nset = nlev**nvar                               ! Given the rounded number of levels, get the appropriate number of sets
    
    allocate(xvar(nvar,nset), xpos(nvar,nset))
    
    do k1=1,nlev
        xpos(1,((k1-1)*nset/nlev + 1):(k1*nset/nlev)) = k1
        xvar(1,((k1-1)*nset/nlev + 1):(k1*nset/nlev)) = lb(1) + k1*(ub(1)-lb(1))/(1.d0 + nlev)
    enddo
    
    do k2=2,nvar
        do k1=1,nlev
            xpos(k2,((k1-1)*nset/nlev + 1):(k1*nset/nlev)) = xpos(k2-1, ::nlev)
            xvar(k2,((k1-1)*nset/nlev + 1):(k1*nset/nlev)) = &
                        lb(1) + xpos(k2-1, ::nlev)*(ub(1)-lb(1))/(1.d0 + nlev)
        enddo
    enddo
    
end subroutine

subroutine lhs_basic(xvar, lb, ub, nset)
implicit none
    double precision, allocatable   :: xvar(:,:)    ! Variables to be distributed
    double precision                :: lb(:), ub(:) ! Lower and upper bounds for variables
    integer                         :: nset         ! Number of variable sets
    integer                         :: k1           ! Iterator
    double precision, allocatable   :: tmp_rand(:)  ! Temporary container for a set of random numbers      
    integer, allocatable            :: sortind(:)   ! Indicies for sorted column of xvar
    integer, allocatable            :: assign_vec(:)! Vector to be assigned (containing the end position of the intervals)
    
    
    allocate(xvar(size(lb),nset), tmp_rand(nset), assign_vec(nset), sortind(nset))
    
    ! Prevent the same numbers to be initiated every time
    call random_seed()
    
    call random_number(xvar)
    
    do k1=1,nset
        assign_vec(k1) = k1
    enddo
    
    do k1=1,size(lb)
        sortind = get_n_smallest_inds(xvar(k1,:), nset)
        call random_number(tmp_rand)
        xvar(k1,sortind) = lb(k1) + (ub(k1)-lb(k1))*(assign_vec-tmp_rand)/nset
    enddo
    
    
end subroutine

subroutine usr_optimization(error_out, xvar_out, error_in, xvar_in, opt, f_data)                  
implicit none
    double precision                :: error_in(:)
    double precision                :: xvar_in(:,:)
    double precision, allocatable   :: error_out(:)
    double precision, allocatable   :: xvar_out(:,:)
    type(fdata_typ), target         :: f_data
    type(fdata_typ), pointer        :: f_data_ptr
    type(c_ptr)                     :: f_data_c_ptr
    type(opt_typ)     		        :: opt
    
    procedure(user_opt_template),pointer :: usr_opt
        integer                         :: nvar
    integer                         :: nset
    integer, allocatable            :: set_pos(:)
    integer 		                :: k1
    double precision                :: lb(size(xvar_in,1)), ub(size(xvar_in,1))
    
    nvar = size(xvar_in,1)
    nset = min(size(error_in), opt%start%num_sets)
    allocate(set_pos(nset))
    set_pos = get_n_smallest_inds(error_in, nset)
    allocate(xvar_out(nvar,nset))
    xvar_out = xvar_in(:,set_pos)
    allocate(error_out(nset))
    error_out = error_in(set_pos)
    
    usr_opt => opt%usr_opt%opt_addr
    f_data_ptr => f_data
    f_data_c_ptr = c_loc(f_data)
    
    lb = get_lower_bounds(f_data%glob)
    ub = get_upper_bounds(f_data%glob)
    
    do k1=1,nset
        call usr_opt(usr_call_sim, xvar_out(:,k1), lb, ub, opt%usr_opt%user_data, error_out(k1), f_data_c_ptr) ! Run optimization
    enddo
    
end subroutine

subroutine usr_call_sim(error, xvar, f_data_c_ptr, evec)
use iso_c_binding
    double precision, intent(out)               :: error                        ! Error from simulation
    double precision, intent(in)                :: xvar(:)                      ! Variables to be simulated 
    type(c_ptr)                                 :: f_data_c_ptr                 ! Data to be used in simulation
    double precision, intent(out), allocatable, optional     :: evec(:)         ! If requested, supplies vector with errors e_i
    type(fdata_typ), pointer                    :: f_data_ptr
    
    call c_f_pointer(f_data_c_ptr, f_data_ptr)
    if (present(evec)) then
        call simulate(error, xvar, f_data_ptr, evec)
    else
        call simulate(error, xvar, f_data_ptr)
    endif
    
end subroutine


end module optimization_mod
    
