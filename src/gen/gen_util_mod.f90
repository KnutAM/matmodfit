! General utility functions to be used around in the program
    
module gen_util_mod
    implicit none
    
    private
    ! Comparing doubles
    public  dbl_compare     ! Compare two float variables with tolerance
    public  dbl_comp_array  ! Compare a float array with a scalar float with tolerance
    public  dbl_comp_array2 ! Compare a float array with another float array with tolerance
    
    ! String tools
    public  upcase, lowcase ! Convert string to upper case only or lower case only
    
    ! Array tools
    public  get_n_smallest_inds ! Find the n smallest indicies in an array (sorted by size small to large)
    
    double precision, parameter     :: def_tol_comp = 1.d-10
    
    contains
    
function dbl_compare(val1, val2, tol_setting) result(is_equal)
! Compare if abs(val1-val2)<tol (Relative tolerance used for values>1)
implicit none
    double precision            :: val1, val2
    double precision, optional  :: tol_setting
    logical                     :: is_equal
    
    double precision            :: tol
    
    if (present(tol_setting)) then
        tol = tol_setting
    else
        tol = 1.d-10
    endif
    
    is_equal = abs(val1-val2)<max(1.d0, abs(val1+val2)/2)*tol
    
end function

function dbl_comp_array(val1, val2, tol_setting) result(is_equal)
! Compare if abs(val1-val2)<tol (Relative tolerance used for values>1)
implicit none
    double precision            :: val1(:), val2
    double precision, optional  :: tol_setting
    logical, allocatable        :: is_equal(:)
    
    double precision            :: tol
    integer                     :: k1
    
    allocate(is_equal(size(val1)))
    
    if (present(tol_setting)) then
        tol = tol_setting
    else
        tol = 1.d-10
    endif
    
    do k1=1,size(val1)
        is_equal(k1) = abs(val1(k1)-val2)<max(1.d0, abs(val1(k1)+val2)/2)*tol
    enddo
    
end function

function dbl_comp_array2(val1, val2, tol_setting) result(is_equal)
! Compare if abs(val1-val2)<tol (Relative tolerance used for values>1)
implicit none
    double precision            :: val1(:), val2(:)
    double precision, optional  :: tol_setting
    logical, allocatable        :: is_equal(:)
    
    double precision            :: tol
    integer                     :: k1
    
    allocate(is_equal(size(val1)))
    
    if (present(tol_setting)) then
        tol = tol_setting
    else
        tol = 1.d-10
    endif
    
    do k1=1,size(val1)
        is_equal(k1) = abs(val1(k1)-val2(k1))<max(1.d0, abs(val1(k1)+val2(k1))/2)*tol
    enddo
    
end function
      

function upcase(string) result(upper)
    character(len=*), intent(in) :: string
    character(len=len(trim(string))) :: upper
    integer :: j
    do j = 1,len(trim(string))
        if(string(j:j) >= "a" .and. string(j:j) <= "z") then
            upper(j:j) = achar(iachar(string(j:j)) - 32)
        else
            upper(j:j) = string(j:j)
        end if
    end do
end function upcase

function lowcase(string) result(lower)
    character(len=*), intent(in) :: string
    character(len=len(trim(string))) :: lower
    integer :: j
    do j = 1,len(trim(string))
        if(string(j:j) >= "A" .and. string(j:j) <= "Z") then
            lower(j:j) = achar(iachar(string(j:j)) + 32)
        else
            lower(j:j) = string(j:j)
        end if
    end do
end function lowcase

function get_n_smallest_inds(x, n) result(inds)
! Given a double precision array x, find the n indicies of the smallest values
implicit none
    double precision    :: x(:)
    integer             :: n
    integer             :: inds(n)
    integer             :: k1
    logical, allocatable:: unused(:)
    allocate(unused(size(x)))
        
    unused = .true.
        
    do k1=1,n
        inds(k1) = minloc(x, dim=1, mask=unused)
        unused(inds(k1)) = .false.
    enddo
end function

end module gen_util_mod
    
