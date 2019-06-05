module opt_inp_mod
use types_mod
use output_mod
use reading_utilities_mod
use check_input_mod
use usr_interface_mod
use load_dll_mod
use nlopt_constants_mod
implicit none
    private
    public  :: read_optimization_settings
    
    contains
    
subroutine read_optimization_settings(opt)
    implicit none
    type(opt_typ), allocatable, intent(out) :: opt(:)
    type(opt_typ), allocatable      :: opt_tmp(:)
    type(usr_opt_typ)               :: usr_opt
    integer                         :: otype, k1

    k1 = 0
    
    allocate(opt_tmp(maxopt))
    
    do while(.not.end_of_category())
        if (adjustl(textline)=='<opt>') then
            k1 = k1 + 1
            call read_int(otype)
            opt_tmp(k1)%otype = otype
            if      (otype==-1) then !User supplied optimization routine
                call get_opt_usr(opt_tmp(k1)%usr_opt)
            elseif  (otype==1) then !NLOPT optimization
                call get_opt_nlopt(opt_tmp(k1))
            elseif  (otype==2) then !Parameter space search
                call get_opt_pspace(opt_tmp(k1))
            else                    !otype not defined
                call close_input()
                call write_output('otype = '//int2str(otype)//' not defined', 'error', 'inp')
            endif    
        else
            call readline()
        endif
    enddo
    
    if (.not.allocated(opt)) allocate(opt(k1))
    opt = opt_tmp(1:k1)
    
end subroutine read_optimization_settings

subroutine get_opt_nlopt(opt)
implicit none
    type(opt_typ)   :: opt
    logical         :: read_next_line

    read_next_line = .true.
    
    do while(status==0)
        if (read_next_line) then
            call readline()
        endif
        read_next_line = .false.
        if (end_of_category('<opt>')) then
            exit
        elseif (adjustl(textline)=='<<start>>') then
            call read_opt_start(opt%start, 1)
        elseif (adjustl(textline)=='<<end_cond>>') then
            call read_opt_end_cond(opt%end_cond)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
            read_next_line = .true.
        endif
    enddo
    
    ! Check input for nlopt optimization
    call check_opt_nlopt(opt)
end subroutine get_opt_nlopt

subroutine get_opt_pspace(opt)
implicit none
    type(opt_typ) :: opt
    logical         :: read_next_line

    read_next_line = .true.
    do while(status==0)
        if (read_next_line) then
            call readline()
        endif
        read_next_line = .false.
        if (end_of_category('<opt>')) then
            exit
        elseif (adjustl(textline)=='<<start>>') then
            call read_opt_start(opt%start, 2)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
            read_next_line = .true.
        endif
    enddo
    
    ! Check input for parameter space "optimization"
    call check_opt_pspace(opt)
end subroutine get_opt_pspace

subroutine read_opt_start(start, otype)
implicit none
type(start_typ)     :: start
integer             :: otype

    do while(status==0)
        call readline()
        if (end_of_subcategory()) then
            exit
        elseif  (adjustl(textline)=='*algorithm') then
            call read_alg(start%algorithm, otype)
        elseif  (adjustl(textline)=='*initial_step') then
            call read_dbl_vector_allocate(start%initial_step)
        elseif  (adjustl(textline)=='*num_sets') then
            call read_int(start%num_sets)
        elseif  (adjustl(textline)=='*fixed_seed') then
            start%use_fixed_seed = .true.
            call read_int(start%fixed_seed)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
        endif
    enddo
end subroutine read_opt_start

subroutine read_opt_end_cond(end_cond)
implicit none
type(end_cond_typ)  :: end_cond

    do while(status==0)
        call readline()
        if (end_of_subcategory()) then
            exit
        elseif  (adjustl(textline)=='*stopval') then
            call read_dbl(end_cond%stopval)
        elseif  (adjustl(textline)=='*ftol_rel') then
            call read_dbl(end_cond%ftol_rel)
        elseif  (adjustl(textline)=='*ftol_abs') then
            call read_dbl(end_cond%ftol_abs)
        elseif  (adjustl(textline)=='*xtol_rel') then
            call read_dbl_vector_allocate(end_cond%xtol_rel)
        elseif  (adjustl(textline)=='*xtol_abs') then
            call read_dbl_vector_allocate(end_cond%xtol_abs)
        elseif  (adjustl(textline)=='*maxeval') then
            call read_int(end_cond%maxeval)
        elseif  (adjustl(textline)=='*maxtime') then
            call read_dbl(end_cond%maxtime)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
        endif
    enddo
end subroutine read_opt_end_cond

subroutine get_opt_usr(usr_opt)
implicit none
    type(usr_opt_typ)                       :: usr_opt

    do while(status==0)
        call readline()
        if (end_of_category('<opt>')) then
            exit
        elseif (adjustl(textline)=='*lib') then
            call read_str(usr_opt%lib, strl)
            call load_user_opt(usr_opt%lib, usr_opt%opt_addr)
        elseif (adjustl(textline)=='*user_data') then
            call read_dbl_vector_allocate(usr_opt%user_data)
        elseif (adjustl(textline)=='<<usr_opt>>') then
            ! Do nothing, as we currently only read in this category there is no need to add another subprocedure...
            ! This just allows the default structure with the category
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
        endif
    enddo
    
    ! Check input for user optimization
    call check_opt_usr(usr_opt)
    
end subroutine get_opt_usr
    
subroutine read_alg(algnr, otype)
implicit none
integer                 :: algnr, otype, stat
character(len=strl)     :: algstring

call read_str(algstring, strl)

read(algstring,*, iostat=stat) algnr

if (stat.ne.0) then
    if (otype==1) then      ! nlopt type
        algnr = get_nlopt_alg(algstring)
    elseif (otype==2) then  ! parameter space search
        algnr = get_pspace_alg(algstring)
    else
        call write_output('Algorithm text option is not supported for optimization type '//int2str(otype), 'error', 'inp')
    endif
endif

end subroutine

function get_pspace_alg(algstring) result(algnr)
implicit none
    character(len=strl)         :: algstring
    integer                     :: algnr
    
    select case (algstring)
    case ('full_factorial')
        algnr = 1
    case ('latin_hypercube')
        algnr = 2
    case default
        algnr = 0
        call write_output('Optimization algorithm "'//algstring//'" is not recognized for parameter space search', 'error', 'inp')
    end select

end function

function get_nlopt_alg(algstring) result(algnr)
    implicit none
    character(len=strl)             :: algstring
    integer                         :: algnr
    
        select case (algstring)
            ! Code below automatically generated from nlopt.f using matlab
        case ('NLOPT_GN_DIRECT')
	        algnr = NLOPT_GN_DIRECT
        case ('NLOPT_GN_DIRECT_L')
	        algnr = NLOPT_GN_DIRECT_L
        case ('NLOPT_GN_DIRECT_L_RAND')
	        algnr = NLOPT_GN_DIRECT_L_RAND
        case ('NLOPT_GN_DIRECT_NOSCAL')
	        algnr = NLOPT_GN_DIRECT_NOSCAL
        case ('NLOPT_GN_DIRECT_L_NOSCAL')
	        algnr = NLOPT_GN_DIRECT_L_NOSCAL
        case ('NLOPT_GN_DIRECT_L_RAND_NOSCAL')
	        algnr = NLOPT_GN_DIRECT_L_RAND_NOSCAL
        case ('NLOPT_GN_ORIG_DIRECT')
	        algnr = NLOPT_GN_ORIG_DIRECT
        case ('NLOPT_GN_ORIG_DIRECT_L')
	        algnr = NLOPT_GN_ORIG_DIRECT_L
        case ('NLOPT_GD_STOGO')
	        algnr = NLOPT_GD_STOGO
        case ('NLOPT_GD_STOGO_RAND')
	        algnr = NLOPT_GD_STOGO_RAND
        case ('NLOPT_LD_LBFGS_NOCEDAL')
	        algnr = NLOPT_LD_LBFGS_NOCEDAL
        case ('NLOPT_LD_LBFGS')
	        algnr = NLOPT_LD_LBFGS
        case ('NLOPT_LN_PRAXIS')
	        algnr = NLOPT_LN_PRAXIS
        case ('NLOPT_LD_VAR1')
	        algnr = NLOPT_LD_VAR1
        case ('NLOPT_LD_VAR2')
	        algnr = NLOPT_LD_VAR2
        case ('NLOPT_LD_TNEWTON')
	        algnr = NLOPT_LD_TNEWTON
        case ('NLOPT_LD_TNEWTON_RESTART')
	        algnr = NLOPT_LD_TNEWTON_RESTART
        case ('NLOPT_LD_TNEWTON_PRECOND')
	        algnr = NLOPT_LD_TNEWTON_PRECOND
        case ('NLOPT_LD_TNEWTON_PRECOND_RESTART')
	        algnr = NLOPT_LD_TNEWTON_PRECOND_RESTART
        case ('NLOPT_GN_CRS2_LM')
	        algnr = NLOPT_GN_CRS2_LM
        case ('NLOPT_GN_MLSL')
	        algnr = NLOPT_GN_MLSL
        case ('NLOPT_GD_MLSL')
	        algnr = NLOPT_GD_MLSL
        case ('NLOPT_GN_MLSL_LDS')
	        algnr = NLOPT_GN_MLSL_LDS
        case ('NLOPT_GD_MLSL_LDS')
	        algnr = NLOPT_GD_MLSL_LDS
        case ('NLOPT_LD_MMA')
	        algnr = NLOPT_LD_MMA
        case ('NLOPT_LN_COBYLA')
	        algnr = NLOPT_LN_COBYLA
        case ('NLOPT_LN_NEWUOA')
	        algnr = NLOPT_LN_NEWUOA
        case ('NLOPT_LN_NEWUOA_BOUND')
	        algnr = NLOPT_LN_NEWUOA_BOUND
        case ('NLOPT_LN_NELDERMEAD')
	        algnr = NLOPT_LN_NELDERMEAD
        case ('NLOPT_LN_SBPLX')
	        algnr = NLOPT_LN_SBPLX
        case ('NLOPT_LN_AUGLAG')
	        algnr = NLOPT_LN_AUGLAG
        case ('NLOPT_LD_AUGLAG')
	        algnr = NLOPT_LD_AUGLAG
        case ('NLOPT_LN_AUGLAG_EQ')
	        algnr = NLOPT_LN_AUGLAG_EQ
        case ('NLOPT_LD_AUGLAG_EQ')
	        algnr = NLOPT_LD_AUGLAG_EQ
        case ('NLOPT_LN_BOBYQA')
	        algnr = NLOPT_LN_BOBYQA
        case ('NLOPT_GN_ISRES')
	        algnr = NLOPT_GN_ISRES
        case ('NLOPT_AUGLAG')
	        algnr = NLOPT_AUGLAG
        case ('NLOPT_AUGLAG_EQ')
	        algnr = NLOPT_AUGLAG_EQ
        case ('NLOPT_G_MLSL')
	        algnr = NLOPT_G_MLSL
        case ('NLOPT_G_MLSL_LDS')
	        algnr = NLOPT_G_MLSL_LDS
        case ('NLOPT_LD_SLSQP')
	        algnr = NLOPT_LD_SLSQP
        case ('NLOPT_LD_CCSAQ')
	        algnr = NLOPT_LD_CCSAQ
        case ('NLOPT_GN_ESCH')
	        algnr = NLOPT_GN_ESCH
            ! End of automatically generated code from nlopt.f using matlab
        case default
            call close_input()
            call write_output('The algorithm '//trim(algstring)//' is not recognized', 'error', 'inp')
        end select
    
end function

end module