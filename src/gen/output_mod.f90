module output_mod
use types_mod
implicit none
    private
    
    public      :: setup_logout, setup_output   ! Setup routines
    public      :: int2str, dbl2str, tim2str    ! Convert to string routines
    public      :: write_output                 ! Main output routine
    public      :: write_sim_results            ! Routine to write out end results after single simulation
    public      :: write_opt_results            ! Routine to write out end results after optimizations
    public      :: write_opt_ana_results        ! Routine to write out results for optimum analysis (error and correlation matrix)
    public      :: write_matrix                 ! Write out a formatted matrix
    public	    :: write_errhist                ! Write out error history to file if requested
    public      :: write_optstage_output        ! Write out the stage results after one optimization (currently in log file)
    
    
    integer  :: matmodfit_log      ! Log file fid
    integer  :: matmodfit_res      ! Result file fid
    integer  :: matmodfit_errhist  ! Error history fid
    
    contains

subroutine set_outname(inputfile, f_data)
implicit none
	character(len=*)				:: inputfile
	type(fdata_typ)                :: f_data
	integer	             			:: pos
	character(len=100)				:: outname
    
    outname = inputfile
    
    ! Remove unix type folder deliminator if there are any
    pos = index(outname, '/', .true.)+1 ! If not found, index returns 0 hence pos is 1
    outname = outname(pos:len(outname)) ! and we do 1:len(simname) which is the entire outname
    
    ! Remove windows type of folder deliminator if there are any
    pos = index(outname, '\', .true.)+1 ! Same as for "/" above
    outname = outname(pos:len(outname))
    
    ! Remove file ending if there is any
    pos = index(outname, '.', .true.)
    if (pos/=0) then
        pos = pos - 1
        outname = outname(1:pos)
    endif
    
    ! Output file name
    write(f_data%glob%outname, "(A)") trim(outname)

end subroutine

subroutine setup_logout(inputfile, f_data)
implicit none
    character(len=*)	:: inputfile
    type(fdata_typ)    :: f_data
    integer             :: status
    character(len=100)	:: outname

    call set_outname(inputfile, f_data)

    outname = f_data%glob%outname
    
    ! Log file
    open(newunit=matmodfit_log, file=trim(outname)//'.log', status='replace', IOSTAT=status)
    if (status/=0) then
        write(*,*) 'matmodfit:error: Could not create logfile "', trim(outname), '.log', '".'
        write(*,*) 'Please check that you have write permissions to the current working directory'
        error stop 1
    endif

end subroutine


subroutine setup_output(f_data)
implicit none
    type(fdata_typ)    :: f_data
    integer             :: status
    character(len=100)	:: outname

    outname = f_data%glob%outname
    
    ! Result file
    open(newunit=matmodfit_res, file=trim(outname)//'.res', status='replace', IOSTAT=status)
    if (status/=0) then
        call write_output('Could not create resultfile "'//trim(outname)//'.res"', 'error')
    endif
    
    ! Error history file if requested
    if (f_data%glob%error_history) then
    	open(newunit=matmodfit_errhist, file=trim(outname)//'.err', status='replace', IOSTAT=status)
    	if (status/=0) then
            call write_output('Could not create error history file "'//trim(outname)//'.err"', 'error')
        endif
        call write_errhist_header()
    endif
end subroutine

function int2str(int) result(str)
implicit none
    integer                         :: int
    character(len=:),allocatable    :: str
    character(len=20)               :: strtmp
    
    write(strtmp,*) int
    str = trim(adjustl(strtmp))
    
end function

function dbl2str(dbl, formatspec) result(str)
implicit none
    double precision                :: dbl
    character(len=*), optional      :: formatspec
    character(len=:), allocatable   :: str
    character(len=20)               :: formatspec_int
    character(len=50)               :: strtmp
    integer                         :: status
    
    if (present(formatspec)) then
        formatspec_int = formatspec
    else
        formatspec_int = 'E15.5E3'
    endif
    
    formatspec_int = '('//trim(formatspec_int)//')'
    write(strtmp,formatspec_int, IOSTAT=status) dbl
    if (status/=0) then
        write(*,*) 'formatspec = ', formatspec_int, ' not ok'
        error stop 1
    endif
    
    str = trim(adjustl(strtmp))
    
end function

function tim2str(dbltime) result(strtime)
implicit none
    double precision                :: dbltime
    character(len=:), allocatable   :: strtime
    integer                         :: hour, minute, second
    double precision                :: time_internal
    
    !Calculate number of hours. Add 0.1 inside int() to avoid numerical errors creating rounding down to wrong value
    hour   = int( dbltime - mod(dbltime, 3600.d0) + 0.1)/3600   
    
    !Calculate number of minutes. Add 0.1 inside int() to avoid numerical errors creating rounding down to wrong value
    time_internal = dbltime - hour*3600 ! Remove already included hours from time
    minute = int( time_internal - mod(time_internal, 60.d0)  + 0.1)/60 
    
    !Calculate number of seconds (round to nearest integer by adding 0.5 inside int())
    second = int(time_internal - minute*60 + 0.5)
    
    !Create time reporting string
    strtime = int2str(hour)//'h:'//int2str(minute)//'m:'//int2str(second)//'s'
    
end function
    
    
subroutine write_output(message_string, output_type, location, loc, halt)
! write_output writes output to both stdout and the <simname>.log file
! Example usages:
! call write_output('something bad happened at line'//int2str(30)//'. Exiting', 'error', 'inp')
! call write_output('did you do something wrong when setting x = '//dbl2str(x)//'?', 'warning')
! Standards for which message types to use
! error:    Program cannot continue and will be terminated
! warning:  Program can maybe continue, but there might be unexpected behavior
! status:   Everything ok, just giving some feedback about what's going on to user

implicit none
    character(len=*)            :: message_string    
    character(len=*), optional  :: output_type, location
    logical, optional           :: loc  !Should location be shown, or just whitespace
    logical, optional           :: halt !Should program stop (default for error, and not for other message types)
    character(len=25)           :: output_type_str, location_str
    logical                     :: show_loc, stop_program
    character(len=50)           :: location_full
    integer                     :: loclen
    
    if (present(output_type)) then
        output_type_str = output_type//':'
    else
        output_type_str = ''
    endif
    
    if (present(location)) then
        location_str = location//':'
    else
        location_str = ''
    endif
    
    if (present(loc)) then
        show_loc = loc
    else
        show_loc = .true.
    endif
    
    if (present(halt)) then
        stop_program = halt
    elseif (output_type_str=='error:') then
        stop_program = .true.
    else
        stop_program = .false.
    endif
    
    
    location_full = 'matmodfit:'//trim(location_str)//trim(output_type_str)
    
    if (show_loc) then
        write(*,'(A)')          trim(location_full)//' '//trim(message_string)
        write(matmodfit_log,'(A)')  trim(location_full)//' '//trim(message_string)
    else
        loclen = len(trim(location_full))
        location_full = ''
        write(*,'(A)')          location_full(1:loclen)//' '//trim(message_string)
        write(matmodfit_log,'(A)')  location_full(1:loclen)//' '//trim(message_string)
    endif
    
    if (stop_program) then
        error stop 1
    endif
    
    
end subroutine

subroutine write_matrix(filenum, matrix, formatspec, zerotol)
    implicit none
    integer, optional           :: filenum
    integer                     :: k1, k2
    double precision, intent(in):: matrix(:,:)
    character(len=*), optional  :: formatspec
    double precision, optional  :: zerotol
    character(len=80)           :: spec
    double precision            :: tol
    logical                     :: default_output
    
    if (.not.present(filenum)) then
        default_output = .true.
    else
        default_output = .false.
    endif
    
    if (present(formatspec)) then
        spec = formatspec
    else
        spec = "(ES15.5E3)"
    endif
    if (present(zerotol)) then
        tol = zerotol
    else
        tol = -1.d0 !Never write zero in text unless specifically requested for some non-negative tol
    endif
    
    
    do k1=1,size(matrix, dim=1)
        do k2=1,size(matrix, dim=2)
            if ((.not.present(formatspec)).and.(abs(matrix(k1,k2)).le.tol)) then
                if (default_output) then
                    write(*, '(A15)', advance="no") 'ZERO'
                else
                    write(filenum, '(A15)', advance="no") 'ZERO'
                endif
                
            else
                if (default_output) then
                    write(*, spec, advance="no") matrix(k1,k2)
                else
                    write(filenum, spec, advance="no") matrix(k1,k2)
                endif
                
            endif
            
        enddo
        if (default_output) then
            write(*, *) ''
        else
            write(filenum, *) ''
        endif
        
    enddo
    
end subroutine write_matrix

subroutine write_vector(filenum, vector, formatspec, zerotol)
    implicit none
    integer, optional           :: filenum
    integer                     :: k1
    double precision, intent(in):: vector(:)
    character(len=*), optional  :: formatspec
    double precision, optional  :: zerotol
    character(len=80)           :: spec
    double precision            :: tol
    logical                     :: default_output
    
    if (.not.present(filenum)) then
        default_output = .true.
    else
        default_output = .false.
    endif
    
    if (present(formatspec)) then
        spec = formatspec
    else
        spec = "(ES15.5E3)"
    endif
    if (present(zerotol)) then
        tol = zerotol
    else
        tol = -1.d0 !Never write zero in text unless specifically requested for some non-negative tol
    endif
    
    
    do k1=1,size(vector, dim=1)
            if ((.not.present(formatspec)).and.(abs(vector(k1)).le.tol)) then
                if (default_output) then
                    write(*, '(A15)', advance="no") 'ZERO'
                else
                    write(filenum, '(A15)', advance="no") 'ZERO'
                endif
                
            else
                if (default_output) then
                    write(*, spec, advance="no") vector(k1)
                else
                    write(filenum, spec, advance="no") vector(k1)
                endif
                
            endif
    enddo
    if (default_output) then
        write(*, *) ''
    else
        write(filenum, *) ''
    endif
    
end subroutine write_vector

subroutine write_sim_results(error, mpar)
implicit none
    double precision            :: error, mpar(:)
    
    write(matmodfit_res,*) '== SIMULATION RESULTS =='
    write(matmodfit_res,*) 'Error after simulation:'
    write(matmodfit_res,"(ES15.5E3)") error
    write(matmodfit_res,*) 'Simulation material parameters:'
    call write_vector(matmodfit_res, mpar)
    write(matmodfit_res,*) '== END OF SIMULATION RESULTS =='

end subroutine

subroutine write_opt_results(errors, ipars)
implicit none
    double precision                :: errors(:), ipars(:,:)
    integer                         :: k1
    character(len=12)               :: dbl_format
    
    dbl_format = '(ES25.15E3)'
    
    write(matmodfit_res,*) '== OPTIMIZATION RESULTS =='
    
    do k1=1,size(errors)
        write(matmodfit_res,*) 'Parameter set '//int2str(k1)//'/'//int2str(size(errors))
        write(matmodfit_res,*) 'Final material parameters:'
        call write_vector(matmodfit_res, ipars(:,k1), dbl_format)
        write(matmodfit_res,*) 'Error after optimization:'
        write(matmodfit_res,dbl_format) errors(k1)
        write(matmodfit_res,*) ''
    enddo
    
    write(matmodfit_res,*) '== END OF OPTIMIZATION RESULTS =='
    
end subroutine

subroutine write_opt_ana_results(error, dfdx, corr_matrix)
implicit none
    double precision    :: error
    double precision    :: dfdx(:)  ! Gradient of objective function wrt. optimized material parameters
    double precision    :: corr_matrix(:,:)
    
    write(matmodfit_res,*) '== OPTIMUM ANALYSIS RESULTS =='
    write(matmodfit_res,*) 'Error at optimum point:'
    write(matmodfit_res,"(ES15.5E3)") error
    write(matmodfit_res,*) 'Error gradient for optimized material parameters:'
    call write_vector(matmodfit_res, dfdx)
    write(matmodfit_res,*) 'Correlation matrix for optimized material parameters:'
    call write_matrix(matmodfit_res, corr_matrix)
    write(matmodfit_res,*) '== END OPTIMUM ANALYSIS RESULTS =='
    
end subroutine

subroutine write_errhist_header()
implicit none
	write(matmodfit_errhist, "(A9, A15, A15)") '% sim nr.', 'error', 'parameters...'
end subroutine

subroutine write_errhist(simnr, error, mpar)
implicit none
	integer, intent(in)				:: simnr
	double precision, intent(in)	:: error
    double precision, intent(in)    :: mpar(:)

    character(len=40)               :: format_string
	
    
    format_string = '(I9, ES15.5E3, '//int2str(size(mpar))//'ES15.5E3)'
    
	write(matmodfit_errhist, format_string) simnr, error, mpar
end subroutine

subroutine write_optstage_output(setnr, error, mpar)
implicit none
	integer, intent(in)				:: setnr
	double precision, intent(in)	:: error
    double precision, intent(in)    :: mpar(:)

    character(len=40)               :: format_string
	
    
    format_string = '(I9, ES15.5E3, '//int2str(size(mpar))//'ES15.5E3)'
    
	write(matmodfit_log, format_string) setnr, error, mpar
    write(*, format_string) setnr, error, mpar
    
end subroutine

    
end module output_mod
