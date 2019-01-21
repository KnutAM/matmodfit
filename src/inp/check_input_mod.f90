! Main reading routine, with the read_inputdata subroutine used to obtain input data    
module check_input_mod
    use reading_utilities_mod
    use output_mod
    use types_mod
    use convert_mpar_mod
    use gen_util_mod
    implicit none
    
    private                         ! All private by default
    
    ! Assign public parts
    public  :: bad_input            ! Routine to notify that bad input entered, and end the program
    public  :: check_glob           ! Check the global settings input
    
    public  :: check_sim_ext        ! Check the simulation settings input for external simulation
    public  :: check_sim_atp        ! Check the simulation settings input for atp simulation
    public  :: check_sim_mps        ! Check the simulation settings input for mp simulation (mps)
    public  :: check_sim_atp_elem_rem!Check the simulation settings input for atp element removal simulation
    public  :: check_sim_usr        ! Check the settings for user simulation
    
    public  :: check_opt_nlopt      ! Check the optimization settings input for nlopt
    public  :: check_opt_pspace     ! Check the optimization settings input for parameter space simulations
    public  :: check_opt_usr
    
    public  :: check_overall        ! Check stuff that can only be checked after all info exists
    
    ! Values for comparisons
    double precision, parameter     :: numtol = 1.d-20
    
    ! Default values for glob arrays
    integer, parameter              ::  DEF_ipar_sctype = 1
    
    ! Default values for sim:atp arrays
    double precision, parameter     :: DEF_err_scale = 1.d0
    integer, parameter, dimension(1):: DEF_result_steps = (/0/)
    integer, parameter, dimension(1):: DEF_error_steps = (/0/)
    double precision, parameter     :: DEF_statev_init = 0.d0
    
    ! Default values for opt arrays
    double precision, parameter, dimension(1)   :: DEF_initial_step = (/0.d0/)  !0=> nlopt default
    double precision, parameter, dimension(1)   :: DEF_xtol_rel = (/0.d0/)      !0=> no tolerance
    double precision, parameter, dimension(1)   :: DEF_xtol_abs = (/0.d0/)      !0=> no tolerance
    
    ! Variables 
    logical     :: nlgeom
    integer     :: nstatv, nchannels
    
    
    contains
    
subroutine bad_input(message)
implicit none
    character(len=*)    :: message
    
    call close_input()
    call write_output(message, 'error', 'inp')
    
end subroutine

subroutine check_glob(glob)
implicit none
    type(glob_typ)                  :: glob
    integer                         :: nparam, noptim, nsets   ! Number of specified parameters
    integer                         :: k1
    logical                         :: lower_bound_nok, upper_bound_nok !Logical variable for checking bounds
    double precision, allocatable   :: mpar(:)
    ! Subroutine tasks:
    ! 0) Assign default run_type if no assigned (needed for later stages)
    ! 1) Check that all mandatory variables have been defined correctly
    ! 2) Check that input is reasonable
    ! 3) Assign default values to optional allocatable variables.
    ! 4) Set variables for module to be used later
    
    if (.not.allocated(glob%run_type)) then
        allocate(glob%run_type(3))
        glob%run_type = [1,2,3]
    endif
    
    ! 1) Check that all mandatory variables have been defined correctly
    if (allocated(glob%ipar_init)) then
        nparam = size(glob%ipar_init)
    else
        call bad_input('ipar_init must be defined')
    endif
    
    if (allocated(glob%ipar_optim)) then
        noptim = size(glob%ipar_optim)
        ! Check if the indicies are within the bounds of the ipar_init array:
        if ((maxval(glob%ipar_optim)>nparam).or.(minval(glob%ipar_optim)<1)) then
            call bad_input('ipar_optim must contain indicies of the ipar_init array')
        endif
    elseif (all(glob%run_type/=1)) then  ! If not a single simulation
        call bad_input('ipar_optim must be defined unless for a single run')
    else                            ! Single simulation
        allocate(glob%ipar_optim(1))
        glob%ipar_optim(1) = 1
        noptim = 1
    endif
    
    if (allocated(glob%ipar_max)) then
        if (noptim/=size(glob%ipar_max)) then
            call bad_input('ipar_max must have the same number of elements as ipar_optim')
        endif
    elseif (all(glob%run_type/=1)) then  ! If not a single simulation
        call bad_input('ipar_max must be defined unless for a single run')
    else                            ! Single simulation
        allocate(glob%ipar_max(1))
        glob%ipar_max(1) = abs(glob%ipar_init(1))*2.d0 + 1.d0 ! Ensure larger (multiply and add to avoid numerical issues when scaling)
    endif
    
    if (allocated(glob%ipar_min)) then
        if (noptim/=size(glob%ipar_min)) then
            call bad_input('ipar_min must have the same number of elements as ipar_optim')
        endif
    elseif (all(glob%run_type/=1)) then  ! If not a single simulation
        call bad_input('ipar_min must be defined unless for a single run')
    else                            ! Single simulation
        allocate(glob%ipar_min(1))
        glob%ipar_min(1) = glob%ipar_init(1)  ! Don't decrement in case log or inverse scaling chosen and initial parameter is close to zero. max increased so they are not equal.
    endif
    
    ! 2) Check that input is reasonable
    ! Check that parameters are within bounds
    lower_bound_nok = any(glob%ipar_init(glob%ipar_optim)<(glob%ipar_min-numtol))
    upper_bound_nok = any(glob%ipar_init(glob%ipar_optim)>glob%ipar_max+numtol)
    if (lower_bound_nok.or.upper_bound_nok) then
        call bad_input('initial parameters must be within range')
    endif
    
    if (any(dbl_comp_array2(glob%ipar_max, glob%ipar_min))) then
        call bad_input('ipar_max cannot be equal to ipar_min, use ipar_optim to not optimize this/these parameters')
    endif
    
    
    if (glob%nstatv<1) then
        call bad_input('Number of state variables must be greater than 0')
    endif
    nstatv = glob%nstatv
    
    ! 3) Assign default values to optional allocatable variables. (Other 
    !    optional variables are set to default values in the type definition)
    
    if (.not.allocated(glob%ipar_sctype)) then
        allocate(glob%ipar_sctype(size(glob%ipar_optim)))
        glob%ipar_sctype = DEF_ipar_sctype
    endif
    
    
    if (.not.allocated(glob%xvar_sets)) then
        allocate(mpar, source=glob%ipar_init)
        if (allocated(glob%ipar_sets)) then
            nsets = size(glob%ipar_sets, 2)
            allocate(glob%xvar_sets(noptim, nsets))
            do k1 = 1,size(glob%xvar_sets,2)
                mpar(glob%ipar_optim) = glob%ipar_sets(:,k1)
                call ipar_to_xvar(mpar, glob%xvar_sets(:,k1), glob)
            enddo
        else
            allocate(glob%ipar_sets(nparam, 1), glob%xvar_sets(noptim, 1))
            call ipar_to_xvar(mpar, glob%xvar_sets(:,1), glob)
        endif
    endif
    
    
    ! 4) Set variables for the module to be used later
    nlgeom = glob%nlgeom
    
    ! 5) Check user subroutines
    if (glob%user_scale) then
        call check_user_scaling(glob)
    endif
    
    
end subroutine
    
subroutine check_user_scaling(glob)
implicit none
    type(glob_typ)                  :: glob
    double precision, allocatable   :: ipar(:), xvar(:), mpar(:)
    
    
    allocate(ipar, source=glob%ipar_init)
    allocate(mpar(size(ipar)))
    allocate(xvar(size(glob%ipar_optim)))
    
    
    call ipar_to_xvar(ipar, xvar, glob)
    call xvar_to_ipar(ipar, xvar, glob)
    
    if (any(.not.dbl_comp_array2(ipar, glob%ipar_init))) then
        call bad_input('User scaling is not reversible')
    endif
    
    ipar(glob%ipar_optim) = glob%ipar_max
    call ipar_to_xvar(ipar, xvar, glob)
    if (any(.not.dbl_comp_array(xvar, 1.d0))) then
        call bad_input('User scaling from ipar_max doesn''t result in 1.d0')
    endif
    call xvar_to_mpar(mpar, xvar, glob)
    if (any(abs(mpar)>huge(1.d0)/10)) then
        call bad_input('User scaling at upper bound gives infinite parameters to umat')
    endif
    
    
    ipar(glob%ipar_optim) = glob%ipar_min
    call ipar_to_xvar(ipar, xvar, glob)
    if (any(.not.dbl_comp_array(xvar, 0.d0))) then
        call bad_input('User scaling from ipar_min doesn''t result in 0.d0')
    endif
    call xvar_to_mpar(mpar, xvar, glob)
    if (any(abs(mpar)>huge(1.d0)/10)) then
        call bad_input('User scaling at lower bound gives infinite parameters to umat')
    endif
    
    
    
    
end subroutine

subroutine check_sim_usr(usr_sim)
implicit none
    type(usr_sim_typ)   :: usr_sim
    
    if (usr_sim%lib=='') then
        call bad_input('A user library path name (*lib) must be specified when running user simulation')
    endif
    
    if (.not.allocated(usr_sim%user_data)) then
        allocate(usr_sim%user_data(1))
        usr_sim%user_data = 1.d0
    endif
    
end subroutine


subroutine check_sim_ext(sim)
implicit none
    type(sim_typ)   :: sim

    if (sim%ext_cmd%script=='') then
        call bad_input('A script command is required when using external simulation')
    endif
    
end subroutine

subroutine check_sim_atp(sim)
implicit none
    type(sim_typ)       :: sim
    
    nchannels = 4

    call check_sim_exp(sim%exp)
    
    call check_sim_mesh1d(sim%mesh1d, sim%exp%ctrl)

    call check_sim_iter(sim%iter)

    call check_sim_err(sim%err)
    
    call check_sim_outp(sim%outp)

end subroutine        

subroutine check_sim_atp_elem_rem(sim)
implicit none
    type(sim_typ)       :: sim
    double precision    :: dummy_ctrl(5,1)
    
    nchannels = 4
    
    call check_sim_iter(sim%iter)
    
    dummy_ctrl(:,1) = (/1.d0, 1.d0, 1.d0, 1.d0, 1.d0/)
    call check_sim_mesh1d(sim%mesh1d, dummy_ctrl)
    
    ! 2) Check that input is reasonable
    if (sim%init%cont_analysis==0) then
        call bad_input('cont_analysis input is mandatory for element removal analysis type. It can not be zero.')
    endif
    
    if (.not.allocated(sim%iter%time_incr)) then
        allocate(sim%iter%time_incr(1,5))
        sim%iter%time_incr(:,1) = (/1.d0, sim%atp_er%time_relx, sim%atp_er%time_relx, sim%atp_er%time_relx, sim%atp_er%time_relx/)
    endif
    
    call check_sim_err(sim%err)
    
    call check_sim_outp(sim%outp)
    
    if (.not.allocated(sim%exp%exp_scale)) then
        allocate(sim%exp%exp_scale(1))
        sim%exp%exp_scale = 1.d0
    endif
    
end subroutine

subroutine check_sim_mps(sim)
implicit none
type(sim_typ)   :: sim
    
    if (nlgeom) then
        nchannels = 9
    else
        nchannels = 6
    endif
    
    call check_sim_exp(sim%exp)
    call check_sim_iter(sim%iter)
    
    ! 2) Check that input is reasonable
    ! Check initial state variables
    if (allocated(sim%init%statev_init)) then
        if (size(sim%init%statev_init,1)/=1) then
            call bad_input('statev_init can only have one row for mps simulation')
        elseif (size(sim%init%statev_init,2)/=nstatv) then
            call bad_input('statev_init must have the same number of elements as specified by nstatv')
        endif
    endif
    
    ! 3) Assign default values to optional allocatable variables.
    call check_sim_err(sim%err)
    
    call check_sim_outp(sim%outp)
    
end subroutine


subroutine check_sim_exp(exp)
implicit none
type(exp_typ)   :: exp
integer         :: k1


    ! Check *ctrl
    if (allocated(exp%ctrl)) then
        if (size(exp%ctrl, 1)/=(nchannels+1)) then
            call bad_input('ctrl should have '//int2str(nchannels+1)//' columns')
        elseif (any(exp%ctrl(1,:)<0)) then
            call bad_input('The first column in ctrl represents steps and must be strictly positive')
        elseif (.not.dbl_compare(exp%ctrl(1,1),1.d0)) then
            if (exp%exp_info(1)==0) then
                call bad_input('The first row in ctrl should give the settings for the first load step')
            endif
        elseif (size(exp%ctrl, 2)>1) then
            k1 = size(exp%ctrl, 2)
            if (any( exp%ctrl(1,2:k1)<=exp%ctrl(1,1:(k1-1)) )) then
                call bad_input('The steps described by ctrl must be in increasing order')
            endif
        endif
    else
        call bad_input('ctrl must be defined')
    endif
    
    ! Check *exp_info
    if (allocated(exp%exp_info)) then
        if (size(exp%exp_info)/=2+2*nchannels+1) then
            call bad_input('exp_info should have '//int2str(3+2*nchannels)//' columns')
        endif
    else
        call bad_input('exp_info must be defined')
    endif
    
    ! Check *exp_data
    if (exp%exp_file=='') then
        call bad_input('exp_data file must be defined')
    endif
    
    ! Check *exp_scale
    if (.not.allocated(exp%exp_scale)) then
        allocate(exp%exp_scale(1))
        exp%exp_scale = 1.d0
    endif

end subroutine

subroutine check_sim_mesh1d(mesh1d, ctrl)
implicit none
type(mesh1d_typ)    :: mesh1d
double precision    :: ctrl(:,:)

    if (allocated(mesh1d%node_pos)) then
        if (any(mesh1d%node_pos<0.d0)) then
            call bad_input('node_pos must be positive')
        elseif (size(mesh1d%node_pos)<2) then
            call bad_input('node_pos: At least two nodes must be specified')
        elseif (any(mesh1d%node_pos(2:size(mesh1d%node_pos))<=mesh1d%node_pos(1:(size(mesh1d%node_pos)-1)))) then
            call bad_input('node_pos must be monotonically strictly increasing')
        endif
    else
        call bad_input('node_pos must be defined')
    endif

    if (dbl_compare(mesh1d%node_pos(1),0.d0,1.d-20)) then
        if (.not.all(dbl_comp_array(ctrl(4,:), 1.d0))) then
            call write_output('Displacement control to zero should be specified on inner radius if it is zero', 'warning', 'inp')
        endif
    endif
    
    ! Full gauss integration requires element_order=2*ngp-1
    if (mesh1d%element_order>(2*mesh1d%ngp - 1)) then
        call write_output('Reduced element integration is specified', 'warning', 'inp')
    elseif (mesh1d%element_order<(2*mesh1d%ngp - 1)) then
        call write_output('Over-integration of element is specified', 'status', 'inp')
    endif
    
end subroutine

subroutine check_sim_iter(iter)
implicit none
type(iter_typ)  :: iter
integer         :: k1

    if (allocated(iter%time_incr)) then
        do k1=1,size(iter%time_incr, 2)
            if (any(iter%time_incr(2:,k1)<1/huge(1.d0))) then
                call bad_input('time_incr: time increment must be greater than zero')
            elseif (any(iter%time_incr(3:4,k1)>iter%time_incr(5,k1))) then
                call bad_input('time_incr: dtmax cannot be smaller than dtmin or dt0')
            elseif (any(iter%time_incr(2:5,k1)<iter%time_incr(3,k1))) then
                call bad_input('time_incr: dtmin cannot be larger than any other increment input')
            endif
        enddo
    else
        call bad_input('time_incr must be defined')
    endif
    
end subroutine

subroutine check_sim_err(err)
implicit none
type(err_typ)   :: err
    
    if (allocated(err%err_scale)) then
        if (size(err%err_scale,1)/=(nchannels+1)) then
            call bad_input('*err_scale should have '//int2str(nchannels+1)//' columns')
        endif
    else
        allocate(err%err_scale((nchannels+1),1))
        err%err_scale(1,1) = 1.d0
        err%err_scale(2:,1) = DEF_err_scale
    endif
    
    if (.not.allocated(err%err_steps)) then
        allocate(err%err_steps(size(DEF_error_steps)))
        err%err_steps = DEF_error_steps
    endif
    
end subroutine

subroutine check_sim_outp(outp)
implicit none
    type(outp_typ)      :: outp
    character(len=50)   :: tmpstr
    integer             :: outpstatus

    if (.not.allocated(outp%result_steps)) then
        allocate(outp%result_steps(size(DEF_result_steps)))
        outp%result_steps = DEF_result_steps
    endif
    
    if (outp%dbl_format=='') then
        outp%dbl_format = '(ES25.15E3)'
    else
        write(tmpstr,outp%dbl_format, iostat=outpstatus) 1.0
        if (outpstatus/=0) then
            call write_output('The dbl_format specified "'//trim(outp%dbl_format)//'" has incorrect syntax.', 'error', 'inp', halt=.false.)
            call write_output('It should follow fortran output formats, e.g. "(ES15.5)"', 'error', 'inp', loc=.false.)
        endif
        call write_output('The dbl_format specified gives the number 1.0 written as:', 'status', 'inp')
        call write_output('"'//trim(tmpstr)//'"', 'status', 'inp', loc=.false.)
    endif

end subroutine

subroutine check_opt_nlopt(opt)
implicit none
    type(opt_typ) :: opt
    
        ! Subroutine tasks:
    ! 1) Check that all mandatory variables have been defined correctly: no mandatory variables
    ! 2) Check that input is reasonable: Nothing to check
    ! 3) Assign default values to optional allocatable variables:
    
    if (.not.allocated(opt%start%initial_step)) then
        allocate(opt%start%initial_step(1))
        opt%start%initial_step = DEF_initial_step !Set to default 0 value to use nlopt default
    endif
    
    if (.not.allocated(opt%end_cond%xtol_rel)) then
        allocate(opt%end_cond%xtol_rel(1))
        opt%end_cond%xtol_rel = DEF_xtol_rel         
    endif
    
    if (.not.allocated(opt%end_cond%xtol_abs)) then
        allocate(opt%end_cond%xtol_abs(1))
        opt%end_cond%xtol_abs = DEF_xtol_abs
    endif
end subroutine

subroutine check_opt_pspace(opt)
implicit none
type(opt_typ)       :: opt

end subroutine check_opt_pspace

subroutine check_opt_usr(usr_opt)
implicit none
    type(usr_opt_typ)   :: usr_opt
    
    if (usr_opt%lib=='') then
        call bad_input('A user library path name (*lib) must be specified when running user optimization')
    endif
    
    if (.not.allocated(usr_opt%user_data)) then
        allocate(usr_opt%user_data(1))
        usr_opt%user_data = 1.d0
    endif
    
end subroutine

subroutine check_overall(f_data, opt)
implicit none
    type(fdata_typ) :: f_data
    type(opt_typ)   :: opt(:)
    ! Internal variables
    integer         :: k1
    
    ! Check if required categories has been read
    if (any(f_data%glob%run_type==2)) then
        if (size(opt)==0) then
            call bad_input('Optimization settings must be defined when a run_type=2 is specified')
        endif
    endif
    
    ! Check if umat is read, if required (i.e. any stype >0)
    if (f_data%glob%umat_lib=='') then
        do k1=1,size(f_data%sim)
            if (f_data%sim(k1)%stype>0) then
                call bad_input('umat must be defined for simulation with stype = ' // int2str(f_data%sim(k1)%stype))
            endif
        enddo
    endif
    
end subroutine

end module    
