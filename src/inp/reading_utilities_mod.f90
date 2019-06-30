! Utilities helping during reading of input files. 
    
module reading_utilities_mod
use output_mod
use gen_util_mod
    implicit none
    
    private
    
    public                  :: readline                                     ! Base routine for reading a textline
    public                  :: textline                                     ! Variable for saving the string data on each line
    public                  :: linecount, status                            ! Status variables (line number and read/convert status)
    public                  :: setup_input, close_input                     ! Initiate and end input reading
    public                  :: maxsim, maxopt                               ! Maximum number of simulations/optimizations
    
    ! Input types reading routines
    public                  :: read_str                                     ! Input format 1
    public                  :: read_logical                                 ! Input format 2
    public                  :: read_int                                     ! Input format 3
    public                  :: read_int_vector, read_int_vector_allocate    ! Input format 4
    public                  :: read_int_mvector                             ! Input format 5
    public                  :: read_dbl                                     ! Input format 6
    public                  :: read_dbl_vector, read_dbl_vector_allocate    ! Input format 7
    public                  :: read_dbl_mvector                             ! Input format 8
    public                  :: read_dblcomp3                                ! Input format 9
    
    ! Checks
    public                  :: end_of_category, end_of_subcategory          ! Has the end of the current category been reached?
    public                  :: is_category                                  ! Determine if current textline is a category
    
    
    integer, parameter      :: textline_length = 1000
    integer, parameter      :: inputfile_length = 80
    integer, parameter      :: max_vector_length = 150
    character, parameter    :: comment_symbol = '!'
    
    integer, parameter      :: maxsim = 10000 ! Maximum number of simulations
    integer, parameter      :: maxopt = 10000 ! Maximum number of simulations

    integer                 :: inpfilenum
    integer                 :: linecount
    integer                 :: status
    character(len=textline_length)  :: textline
    
    contains

subroutine setup_input(inputfile)
implicit none
character(len=*) :: inputfile

    linecount   = 1
    status      = 0
    open(newunit=inpfilenum, file=inputfile, status='old', action='read', IOSTAT=status)
    if (status/=0) then !Check if file can be opened.
        call write_output('Inputfile "'//trim(inputfile)//'" cannot be found', 'error', 'inp')
    endif
    
end subroutine

subroutine close_input()
implicit none
    close(inpfilenum)
end subroutine

subroutine error_on_line()
implicit none
    call close_input()
    call write_output('Input file format error on line '//int2str(linecount), 'error', 'inp')
    
end subroutine

subroutine readline()
    implicit none
    textline = ""
    do while (textline=="") !Don't use lines that are completely empty (with exception on what is behind the comment symbol)
        read(inpfilenum, "(A)", IOSTAT=status) textline
        if (status>0) then
            call error_on_line()
        elseif (status<0) then
            textline="End of file"    !End of file
        else
            call remove_comment()
            linecount = linecount + 1
        endif
    enddo
    
end subroutine readline

subroutine remove_comment()
    implicit none
    integer                     :: k1
    if (textline(1:1)==comment_symbol) then
        textline = ""
    else
        do k1=2,(textline_length-1)
            if (textline(k1:k1)==comment_symbol) then
                exit
            endif
        enddo
        textline = textline(1:(k1-1))
    endif
end subroutine remove_comment

subroutine find_length_int(arraylen, tmp_int)
implicit none
    integer                         :: arraylen
    integer                         :: check_int
    integer                         :: tmp_int(max_vector_length+1)
    
    ! Check that textline is ok in format
    read(textline, *, iostat=status) check_int  ! Use check_int as temp variable for now
    if (status==0) then    
        status = 0
        check_int = -1234   ! Some unlikely number
        tmp_int = 0
        read(textline, *, iostat=status) tmp_int
        ! Iterate until no input matches the guessed unlikely number
        do while (any(tmp_int==check_int))
            check_int = check_int + 1
        enddo
        ! When ok, assign the number that doesn't match any input and read the input
        tmp_int = check_int
        read(textline, *, iostat=status) tmp_int
        if (any(status==(/-1, 0/))) then    ! End of file (-1) expected if line shorter than max_vector_length+1
            status = 0
            arraylen = max_vector_length+1 - count(tmp_int==check_int)
            if (arraylen==max_vector_length+1) then
                call write_output('Max number of elements in vector input is '//int2str(max_vector_length)//'(line '//int2str(linecount)//')', 'warning', 'inp')
                arraylen=max_vector_length
            endif
        else
            call error_on_line()
        endif
    else
        call error_on_line()
    endif
    
end subroutine

subroutine find_length_dbl(arraylen, tmp_dbl)
implicit none
    integer                         :: arraylen
    double precision                :: check_dbl
    double precision                :: tmp_dbl(max_vector_length+1)
    
    ! Check that textline is ok in format
    read(textline, *, iostat=status) check_dbl  ! Use check_int as temp variable for now
    if (status==0) then
        check_dbl = -1234.5678   ! Some unlikely number
        tmp_dbl = 0
        read(textline, *, iostat=status) tmp_dbl
        ! Iterate until no input matches the guessed unlikely number
        do while (any(dbl_comp_array(tmp_dbl, check_dbl)))
            check_dbl = check_dbl + 1.d0
        enddo
        ! When ok, assign the number that doesn't match any input and read the input
        tmp_dbl = check_dbl
        read(textline, *, iostat=status) tmp_dbl
        if (any(status==(/-1, 0/))) then    ! End of file (-1) expected if line shorter than max_vector_length+1
            status = 0
            arraylen = max_vector_length+1 - count(dbl_comp_array(tmp_dbl, check_dbl))
            if (arraylen==max_vector_length+1) then
                call write_output('Max number of elements in vector input is '//int2str(max_vector_length)//'(line '//int2str(linecount)//')', 'warning', 'inp')
                arraylen=max_vector_length
            endif
        else
            call error_on_line()
        endif
    else
        call error_on_line()
    endif
    
end subroutine

! Input format 1
subroutine read_str(string, strlen)
    implicit none
    integer                         :: strlen
    character(len=strlen)           :: string

    call readline()         ! Read line into textline
    if (status==0) then     ! Check that reading was OK
        string = textline   ! Save textline in the "string" variable
    else
        call error_on_line()
    endif
end subroutine read_str

! Input format 2
subroutine read_logical(logvar)
implicit none
    logical                     :: logvar
    integer                     :: intval
    
    call readline()         ! Read the line in the file
    if (status==0) then     ! Check that reading went OK and line put into textline
        read(textline, *, IOSTAT=status) logvar! Convert textline into integer
        if (status/=0) then                     ! Check that conversion went OK
            read(textline, *, IOSTAT=status) intval
            if (status==0) then
                if (intval==0) then
                    logvar = .false.
                elseif (intval==1) then
                    logvar = .true.
                else
                    call error_on_line()
                endif
            else
                call error_on_line()
            endif
        endif
    else
        call error_on_line()
    endif
    
end subroutine

! Input format 3
subroutine read_int(intval)
    implicit none
    integer                     :: intval
    
    call readline()         ! Read the line in the file
    if (status==0) then     ! Check that reading went OK and line put into textline
        read(textline, *, IOSTAT=status) intval ! Convert textline into integer
        if (status/=0) then                     ! Check that conversion went OK
            call error_on_line()
        endif
    else
        call error_on_line()
    endif
    
end subroutine read_int

! Input format 4 - known length
subroutine read_int_vector(intval)
    implicit none
    integer                         :: intval(:), tmpint
    
    call readline() ! Save the file line into textline
    if (status==0) then
        read(textline, *, IOSTAT=status) intval ! Convert the textline into an integer array
        ! If length of intval is longer than input, status becomes -1. 
        ! To not interpret this as end of file, check for a scalar int:
        ! This should perhaps be used to generate a warning as not all values are set
        if (status==-1) then
            read(textline, *, IOSTAT=status) tmpint
        endif
        if (status/=0) then ! Check that the conversion seem to work
            call error_on_line()
        endif
    else
        call error_on_line()
    endif
end subroutine read_int_vector

! Input format 4 - unknown length
subroutine read_int_vector_allocate(intval)
    implicit none
    integer, allocatable            :: intval(:)
    integer                         :: arraylen
    integer                         :: tmp_int(max_vector_length+1)
    
    call readline()     ! Read in the file line into textline
    if (status==0) then ! Check status of reading the file line
        call find_length_int(arraylen, tmp_int) ! Find the number of integers on the line
        allocate(intval(arraylen))  
        intval = tmp_int(1:arraylen)    ! Save the final integer array
    endif

end subroutine read_int_vector_allocate

! Input format 5 - unknown number of columns, number of rows specified as a first integer
subroutine read_int_mvector(intval)
implicit none
    integer, allocatable            :: intval(:,:)
    integer                         :: nrows, k1
    integer, allocatable            :: tmp_int(:)
    
    call read_int(nrows)
    call read_int_vector_allocate(tmp_int)
    allocate(intval(size(tmp_int),nrows))
    intval(:,1) = tmp_int
    do k1=2,nrows
        call read_int_vector(intval(:,k1))
    enddo
    
end subroutine

! Input format 6 - single double input
subroutine read_dbl(dblval)
    implicit none
    double precision            :: dblval
    
    call readline()     ! Read in the file line into textline
    if (status==0) then ! Check status of reading the file line
        read(textline, *, IOSTAT=status) dblval ! Convert the textline to a double
        if (status/=0) then
            call error_on_line()
        endif
    else
        call error_on_line()
    endif
end subroutine read_dbl

! Input format 7 - vector double input (known size)
subroutine read_dbl_vector(dblval)
    implicit none
    double precision                :: dblval(:), tmpdbl
    call readline()
    if (status==0) then
        read(textline, *, IOSTAT=status) dblval
        ! If length of dblval is longer than input, status becomes -1. 
        ! To not interpret this as end of file, check for a scalar double:
        if (status==-1) then
            read(textline, *, IOSTAT=status) tmpdbl
        endif
        
        if (status/=0) then
            call error_on_line()
        endif
    else
        call error_on_line()
    endif
end subroutine read_dbl_vector

! Input format 7 - vector double input (unknown size - allocate)
subroutine read_dbl_vector_allocate(dblval)
    implicit none
    double precision, allocatable   :: dblval(:)
    integer                         :: arraylen
    double precision                :: tmp_dbl(max_vector_length+1)
    
    call readline()
    if (status==0) then
        call find_length_dbl(arraylen, tmp_dbl)
        allocate(dblval(arraylen))
        dblval = tmp_dbl(1:arraylen)
    else
        call error_on_line()
    endif

end subroutine read_dbl_vector_allocate

! Input format 8 - multiple double vectors (unknown number of cols, specified number of rows)
subroutine read_dbl_mvector(dblval, ncols, defval)
implicit none
    double precision, allocatable   :: dblval(:,:)
    integer                         :: nrows, k1, k0
    double precision, allocatable   :: tmp_dbl(:)
    integer, optional               :: ncols    ! Number of columns if known (to allow different number of columns per row)
    double precision, optional      :: defval   ! Default value, only if ncols is used
    
    call read_int(nrows)
    if (present(ncols)) then
        k0 = 1
        allocate(dblval(ncols, nrows))
        if (present(defval)) then
            dblval=defval
        endif
    else
        call read_dbl_vector_allocate(tmp_dbl)
        allocate(dblval(size(tmp_dbl), nrows))
        dblval(:,1) = tmp_dbl
        k0 = 2
    endif
    
    do k1=k0,nrows
        call read_dbl_vector(dblval(:,k1))
    enddo
    
end subroutine

! Input format 9 - 3 dimensional double array
subroutine read_dblcomp3(dblarray)
    implicit none
    double precision                :: dblarray(:,:,:), temp(size(dblarray, 3))
    integer                         :: k1, k2
    logical                         :: incategory
    dblarray = 0.d0
    
    incategory = .true.
    do while (incategory)
        call readline()
        if ((textline(2:2)/='*').and.(status==0)) then
            temp = 0.d0
            read(textline, *, IOSTAT=status) k1, k2, temp
            if (status==0) then
                incategory = .true.
                dblarray(k1,k2,:) = temp
            else
                call error_on_line()
            endif
        else
            incategory = .false.
        endif
    enddo

end subroutine read_dblcomp3

! Check if end of category has been reached
function end_of_category(keystr1, keystr2, testlength)
! Determine if we are at the end of one category (or transitioning to new number of same by having textline=keystring)
implicit none
character(len=*), optional  :: keystr1
character(len=*), optional  :: keystr2
integer, optional           :: testlength
logical                     :: end_of_category
character(len=len(textline)):: testline

if (present(testlength)) then
    testline = textline(1:testlength)
else
    testline = textline
endif

if (present(keystr2)) then
    end_of_category = (testline==keystr1).or.(testline==keystr2).or.(status/=0)
elseif (present(keystr1)) then
    end_of_category = (testline==keystr1).or.(status/=0)
else
    end_of_category = (status/=0)
endif

end function

function end_of_subcategory()
implicit none
    logical :: end_of_subcategory
    character(len=len(textline)):: tmpline
    tmpline = adjustl(textline)
    end_of_subcategory = (tmpline(1:1)=='<').or.(status/=0)
end function

function is_category()
implicit none
    logical :: is_category
    is_category = textline(1:1)=='<'
end function

end module
