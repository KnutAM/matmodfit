! Specifics for windows
    
module dll_mod
    use, intrinsic :: iso_c_binding    
    use kernel32    ! Allows use of LoadLibrary, GetProcAddress and FreeLibrary on Windows
    use gen_util_mod
    use output_mod
    implicit none   

    private
    public load_user_library
    
    character(len=4), parameter     :: lib_ext = '.dll'
    
    contains

    subroutine load_user_library(fun_cptr, path_name, fun_name)
    implicit none
        type(c_funptr), allocatable         :: fun_cptr(:)  ! C-type pointer to the dll's subroutine
        integer(c_intptr_t)                 :: fun_tmp      ! Temporary integer with pointer location
        character(len=*)         			:: path_name    ! dynamically linked library name
        character(len=*)                    :: fun_name(:)  ! Base function name (without underscores)
        
        integer(handle)                     :: dll_handle   ! C-type handle to entire dll
        character(len=100)		            :: fun_mangle(4)! Different attempts of name-mangling scheme to find function
        integer					            :: k1           ! Loop through the fun_mangle variable
        integer                             :: nfun, kf     ! Number of functions to load from dll, and looping iterator for the same
        
        nfun = size(fun_name)
        allocate(fun_cptr(nfun))
        
        dll_handle = LoadLibrary( trim(path_name)//trim(lib_ext)//c_null_char )
        if (dll_handle==0) then
            call write_output('Cannot locate "'//trim(path_name)//trim(lib_ext)//'"', 'error', 'inp', halt=.false.)
            call write_output('Check that you have compiled it as 64 bit', 'error', 'inp', halt=.false.)
            call write_output('Please verify that it is available in path', 'error', 'inp', loc=.false.)
            
        else
            do kf=1,nfun
                fun_mangle(1) = upcase(fun_name(kf))
                fun_mangle(2) = trim(fun_mangle(1))//'_'
                fun_mangle(3) = lowcase(fun_name(kf))
                fun_mangle(4) = trim(fun_mangle(3))//'_'
                
    		    do k1=1,size(fun_mangle)
    			    fun_tmp = GetProcAddress(dll_handle, trim(fun_mangle(k1))//c_null_char)
    			    if (fun_tmp/=0) then
    				    exit
				    endif
			    enddo
			    if (fun_tmp==0) then
                    call write_output('Cannot find subroutine '//fun_name(kf)//' in "'//trim(path_name)//trim(lib_ext)//'"', 'error', 'inp', halt=.false.)
				    call write_output('Please verify that '//fun_name(kf)//' exists using "dumpbin /exports '//trim(path_name)//trim(lib_ext), 'error', 'inp', loc=.false.)
                endif
                fun_cptr(kf) = transfer(fun_tmp, fun_cptr(kf))
            enddo
        endif
        
    end subroutine
        
end module dll_mod
    
