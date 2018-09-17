! Specifics for linux
! Provide the subroutine load_umat(lib_name, umat_ptr)
module dll_mod
	use, intrinsic :: iso_c_binding
    use gen_util_mod
    use output_mod
    implicit none
    
    private
    public load_user_library
    
    character(len=3), parameter :: lib_ext = '.so'
    integer(c_int), parameter	:: RTLD_LAZY = 1
    integer(c_int), parameter   :: RTLD_NOW  = 2
    CHARACTER(C_CHAR), DIMENSION(1), SAVE, TARGET, PRIVATE :: unknown_errormsg="Unknown error"

   ! interface to linux API (https://rosettacode.org/wiki/Call_a_function_in_a_shared_library#GNU_Fortran_on_Linux)
   interface
      function dlopen(filename,mode) bind(c,name="dlopen")
         ! void *dlopen(const char *filename, int mode);
         use iso_c_binding
         implicit none
         type(c_ptr) :: dlopen
         character(c_char), intent(in) :: filename(*)
         integer(c_int), value :: mode
      end function
 
      function dlsym(handle,name) bind(c,name="dlsym")
         ! void *dlsym(void *handle, const char *name);
         use iso_c_binding
         implicit none
         type(c_funptr) :: dlsym
         type(c_ptr), value :: handle
         character(c_char), intent(in) :: name(*)
      end function
 
      function dlclose(handle) bind(c,name="dlclose")
         ! int dlclose(void *handle);
         use iso_c_binding
         implicit none
         integer(c_int) :: dlclose
         type(c_ptr), value :: handle
      end function
      
      function dlerror() bind(c,name="dlerror")
        use iso_c_binding
        implicit none
        type(c_ptr)   :: dlerror
      end function
      
      
        
   end interface
    
    contains
    
    function c_2_f_string(c_str_ptr) RESULT(f_str_ptr)
    implicit none
    ! Source: https://cims.nyu.edu/~donev/Fortran/DLL/DLL.Forum.txt
    ! Convert a null-terminated C string into a Fortran character array pointer
        type(c_ptr), intent(in)   :: c_str_ptr ! Pointer to c-string
        character(kind=c_char), dimension(:), pointer :: f_str_ptr

        interface ! strlen is a standard C function from <string.h>
            ! int strlen(char *string)
            function strlen(string) bind(c,name="strlen")
                use iso_c_binding
                implicit none
                type(c_ptr), value :: string
                integer            :: strlen
            end function
        end interface   

        if(c_associated(c_str_ptr)) then
            call c_f_pointer(fptr=f_str_ptr, cptr=c_str_ptr, shape=[strlen(c_str_ptr)])
        else
            f_str_ptr=>unknown_errormsg
        endif

    end function
    
    
    subroutine load_user_library(fun_cptr, path_name, fun_name)
    implicit none
        type(c_funptr), allocatable :: fun_cptr(:)  ! C-type pointer to the dll's subroutine
        character(len=*)            :: path_name    ! dynamically linked library name
        character(len=*)            :: fun_name(:)  ! Base function name (without underscores)
        
        type(c_ptr)                 :: dll_handle   ! C-type handle to entire so (dll on win)
        character(len=100)		    :: fun_mangle(4)! Different attempts of name-mangling scheme to find function
        integer					    :: k1           ! Loop through the fun_mangle variable
        integer                     :: nfun, kf     ! Number of functions to load from dll, and looping iterator for the same        
        character(len=500)          :: errormsg
        
        nfun = size(fun_name)
        allocate(fun_cptr(nfun))
        
        dll_handle = c_null_ptr
        fun_cptr = c_null_funptr
        
        dll_handle = dlopen(trim(path_name)//trim(lib_ext)//c_null_char, RTLD_NOW)
        if (.not.c_associated(dll_handle)) then
            write(errormsg,*) c_2_f_string(dlerror())
            call write_output('dlopen error: '//trim(errormsg), 'error', 'inp')
    	else
            do kf=1,nfun
                fun_mangle(1) = upcase(fun_name(kf))
                fun_mangle(2) = trim(fun_mangle(1))//'_'
                fun_mangle(3) = lowcase(fun_name(kf))
                fun_mangle(4) = trim(fun_mangle(3))//'_'
                
    		    do k1=1,size(fun_mangle)
    			    fun_cptr(kf) = dlsym(dll_handle, trim(fun_mangle(k1))//c_null_char)
			        if (c_associated(fun_cptr(kf))) then
				        exit
			        endif
			    enddo
			    if (.not.c_associated(fun_cptr(kf))) then
			        write(errormsg,*) c_2_f_string(dlerror())
                    call write_output('dlsym error: '//trim(errormsg), 'error', 'inp')
    			endif
            enddo
        endif
        
    end subroutine
    
end module dll_mod
    
