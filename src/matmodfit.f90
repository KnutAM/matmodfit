 
program matmodfit 
    use types_mod
    use input_mod
    use output_mod
    use run_types_mod
    use simulation_mod
    implicit none

    character(len=80), parameter    :: version_string = 'matmodfit version 1.1'
    type(fdata_typ)                 :: f_data
    type(opt_typ), allocatable      :: opt(:)
    character(len=80)               :: inputfile
    double precision, allocatable   :: xvars(:,:)
    integer                         :: k1
    
    ! Get input file name
    call getarg(1, inputfile)
    if (inputfile==' ') then
        inputfile = 'matmodfit.inp' !Set default input file
    elseif (inputfile=='--version') then
        write(*,"(A)") trim(version_string)
        stop
    elseif (inputfile=='--help') then
        write(*,"(A)") trim(version_string)
        write(*,"(A)") 'To run an input file: "matmodfit inputfilename.inp"'
        write(*,"(A)") 'It is also possible give a full path if inputfilename'
        write(*,"(A)") 'is not in current directory'
        write(*,"(A)") 'For more information please see the full manual'
        stop
    endif
    
    
    ! Setup log output (must be done first)
    call setup_logout(inputfile, f_data)        

    ! Get input data from inputfile
    call read_inputdata(inputfile, f_data, opt)

    ! Setup output
    call setup_output(f_data)

    ! Setup simulations
    call setup_simulations(f_data)
    
    ! Extract initial material parameters
    allocate(xvars, source=f_data%glob%xvar_sets)
    
    do k1=1,size(f_data%glob%run_type)
        ! Choose run type
        select case (f_data%glob%run_type(k1))
        case (1)
            call run_simulation(f_data, xvars)
        case (2)
            call run_optimization(f_data, opt, xvars)
        case (3)
            call run_optanalyzer(f_data, xvars)
        case default
            call write_output('Run type '//int2str(f_data%glob%run_type(k1))//' not supported', 'error')
        end select
    enddo


end program