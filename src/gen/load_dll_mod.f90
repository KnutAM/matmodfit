module load_dll_mod
use, intrinsic :: iso_c_binding
use dll_mod
use umat_mod
use usr_interface_mod

    private
    public  :: load_umat
    public  :: load_user_scaling
    public  :: load_user_sim
    public  :: load_user_opt
    public  :: load_user_error
    
    contains
    
    subroutine load_umat(path_name, umat)
    implicit none
        character(len=*)         			:: path_name     ! umat dynamically linked library name
        procedure (umat_template), pointer  :: umat			! Procedure pointer to umat
        type(c_funptr), allocatable         :: umat_cptr(:) ! C-type pointer to the dll's umat subroutine
        character(len=5)		            :: fun_name(1)  ! umat name
        
        fun_name(1) = 'umat'
        call load_user_library(umat_cptr, path_name, fun_name)
        
        call c_f_procpointer(transfer(umat_cptr(1), c_null_funptr), umat)

    end subroutine
    
    subroutine load_user_scaling(path_name, xvar_to_mpar_addr, ipar_to_xvar_addr, xvar_to_ipar_addr)
    implicit none
        character(len=*)         			:: path_name     ! dynamically linked library name
        procedure(user_xvar_to_mpar_template),pointer :: xvar_to_mpar_addr ! Pointer address to user_xvar_to_mpar subroutine
        procedure(user_ipar_to_xvar_template),pointer :: ipar_to_xvar_addr ! Pointer address to user_ipar_to_xvar subroutine
        procedure(user_xvar_to_ipar_template),pointer :: xvar_to_ipar_addr ! Pointer address to user_xvar_to_mpar subroutine
        type(c_funptr), allocatable         :: fun_cptr(:)  ! C-type pointers to the dll's subroutines
        character(len=20)                   :: fun_name(3)  ! Names of functions
        
        fun_name(1) = 'xvar_to_mpar'
        fun_name(2) = 'ipar_to_xvar'
        fun_name(3) = 'xvar_to_ipar'
        
        call load_user_library(fun_cptr, path_name, fun_name)
		
        call c_f_procpointer(transfer(fun_cptr(1), c_null_funptr), xvar_to_mpar_addr)
        call c_f_procpointer(transfer(fun_cptr(2), c_null_funptr), ipar_to_xvar_addr)
        call c_f_procpointer(transfer(fun_cptr(3), c_null_funptr), xvar_to_ipar_addr)
        
    end subroutine
    
    subroutine load_user_sim(path_name, sim_addr)
    implicit none
        character(len=*)                    :: path_name    ! dynamically linked library name
        procedure(user_sim_template),pointer:: sim_addr     ! Pointer address to user_sim subroutine
        type(c_funptr), allocatable         :: fun_cptr(:)  ! C-type pointers to the dll's subroutines
        character(len=20)                   :: fun_name(1)  ! Names of functions
        
        fun_name(1) = 'sim'
        
        call load_user_library(fun_cptr, path_name, fun_name)
		
        call c_f_procpointer(transfer(fun_cptr(1), c_null_funptr), sim_addr)
        
    end subroutine
    
    subroutine load_user_opt(path_name, opt_addr)
    implicit none
        character(len=*)                    :: path_name    ! dynamically linked library name
        procedure(user_opt_template),pointer:: opt_addr     ! Pointer address to user_opt subroutine
        type(c_funptr), allocatable         :: fun_cptr(:)  ! C-type pointers to the dll's subroutines
        character(len=20)                   :: fun_name(1)  ! Names of functions
        
        fun_name(1) = 'opt'
        
        call load_user_library(fun_cptr, path_name, fun_name)
		
        call c_f_procpointer(transfer(fun_cptr(1), c_null_funptr), opt_addr)
        
    end subroutine
    
    subroutine load_user_error(path_name, err_addr)
    implicit none
        character(len=*)                    :: path_name    ! dynamically linked library name
        procedure(user_error_template),pointer:: err_addr     ! Pointer address to user_opt subroutine
        type(c_funptr), allocatable         :: fun_cptr(:)  ! C-type pointers to the dll's subroutines
        character(len=20)                   :: fun_name(1)  ! Names of functions
        
        fun_name(1) = 'uerror'
        
        call load_user_library(fun_cptr, path_name, fun_name)
		
        call c_f_procpointer(transfer(fun_cptr(1), c_null_funptr), err_addr)
        
    end subroutine
    
end module
