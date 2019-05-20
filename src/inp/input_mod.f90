! Main reading routine, with the read_inputdata subroutine used to obtain input data    

module input_mod
    use reading_utilities_mod   ! Defines status, textline and readcount as well as functions for reading the different datatypes
    use check_input_mod     ! Module for checking that the input is correct
    use output_mod          ! Module for producing output log messages
    use types_mod
    use sim_inp_mod
    use opt_inp_mod
    use load_dll_mod
    implicit none
    
    private                         ! All private by default
    public  :: read_inputdata       ! Only the main routine, to read the input data, is public
    
    contains

subroutine read_inputdata(inputfile, f_data, opt)
    implicit none
    character(len=*)            :: inputfile    ! Full or relative path to input file
    type(fdata_typ)             :: f_data       ! 
    type(opt_typ),allocatable   :: opt(:)
    ! End of header
    
    call write_output('Reading of input from "'//trim(inputfile)//'" is starting', 'status', 'inp')
    
    call setup_input(inputfile)
    
    ! File should always start with the global settings
    call read_global_settings(f_data%glob)
    
    ! Then comes the simulations
    call read_simulation_settings(f_data%sim)

    
    ! Finally the optimization input
    call read_optimization_settings(opt)
    
    ! Check complete input data (things that cannot be checked locally)
    call check_overall(f_data, opt)
    
    call close_input()
    
    call write_output('Reading of input from "'//trim(inputfile)//'" has completed', 'status', 'inp')
    call write_output(' ', 'status', 'inp', loc=.false.) ! Write out empty line after input reading finished.
    
end subroutine read_inputdata

! Global settings
subroutine read_global_settings(glob)
    implicit none
    type(glob_typ)                  :: glob

    do while(status==0)
        call readline()
        if (end_of_category('<sim>')) then
            exit
        elseif (adjustl(textline)=='*run_type') then
            call read_int_vector_allocate(glob%run_type)
        elseif (adjustl(textline)=='*umat_lib') then
            call read_str(glob%umat_lib, strl)
            call load_umat(glob%umat_lib, glob%umat_address)
        elseif (adjustl(textline)=='*cmname') then
            call read_str(glob%cmname, strl)
        elseif (adjustl(textline)=='*nlgeom') then
            call read_logical(glob%nlgeom)
        elseif (adjustl(textline)=='*nstatv') then
            call read_int(glob%nstatv)
        elseif (adjustl(textline)=='*ipar_init') then
            call read_dbl_vector_allocate(glob%ipar_init)
        elseif (adjustl(textline)=='*ipar_optim') then
            call read_int_vector_allocate(glob%ipar_optim)
        elseif (adjustl(textline)=='*ipar_sctype') then
            call read_int_vector_allocate(glob%ipar_sctype)
        elseif (adjustl(textline)=='*ipar_max') then
            call read_dbl_vector_allocate(glob%ipar_max)
        elseif (adjustl(textline)=='*ipar_min') then
            call read_dbl_vector_allocate(glob%ipar_min)
        elseif (adjustl(textline)=='*ipar_sets') then
            call read_dbl_mvector(glob%ipar_sets)
        elseif (adjustl(textline)=='*xvar_sets') then
            call read_dbl_mvector(glob%xvar_sets)
        elseif (adjustl(textline)=='*num_grad_pert') then
            call read_dbl(glob%num_grad_pert)
        elseif (adjustl(textline)=='*error_history') then
            call read_logical(glob%error_history)
        elseif (adjustl(textline)=='*user_scale_lib') then
            call read_str(glob%user_scale_lib, strl)
            glob%user_scale = .true.
            call load_user_scaling( glob%user_scale_lib, &
                                    glob%xvar_to_mpar_addr, &
                                    glob%ipar_to_xvar_addr, &
                                    glob%xvar_to_ipar_addr)
        elseif (adjustl(textline)=='*opt_resnr') then
            call read_int(glob%opt_resnr)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
        endif
    enddo
    
    call check_glob(glob)   ! Check that everything is assigned, and if neccessary assign default allocatable variables to defaults
end subroutine

end module input_mod
