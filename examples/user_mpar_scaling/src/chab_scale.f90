! Example for the chaboche model
! mpar: E, nu, Y0, Hiso, 1/kinf, Hk1, 1/ainf1, Hk2, 1/ainf2, ....
! ipar: E, nu, Y0, Htot, Ysat, rk1, ra1, rk2, ra2, ...
! Where
!Htot = Hiso + sum(Hki)
!Ysat = Y0 + kinf + sum(ainfi)
!rki = Hki/Htot
!rai = ainfi/(Ysat-Y0)
! Scale all parameters linearly

module util_funs
implicit none
    private
    public :: lin_scale
    public :: lin_scale_back
    
    contains

    function lin_scale(val, val_max, val_min)
    implicit none
        double precision    :: val, val_max, val_min
        double precision    :: lin_scale
        
        lin_scale = (val-val_min)/(val_max-val_min)
    end function

    function lin_scale_back(val_scaled, val_max, val_min) result(val)
    implicit none
        double precision    :: val_scaled, val_max, val_min
        double precision    :: val
        
        val = (val_max-val_min)*val_scaled + val_min
    end function

end module
    
subroutine xvar_to_mpar(mpar, ipar, xvar, optim, ipar_min, ipar_max)
!DEC$ ATTRIBUTES DLLEXPORT :: xvar_to_mpar
use util_funs
implicit none
    double precision, intent(inout) :: mpar(:)      ! The material parameters (used by umat): To be updated
    double precision, intent(inout) :: ipar(:)      ! Contains the initial input parameters.
    double precision, intent(in)    :: xvar(:)      ! Contains scaled parameters [0,1] that are optimized
    integer, intent(in)             :: optim(:)     ! Contains the position of the optimized parameters in ipar
    double precision, intent(in)    :: ipar_min(:)  ! Min values for input material parameters
    double precision, intent(in)    :: ipar_max(:)  ! Max values for input material parameters
    ! Internal variables
    integer                         :: ind, k1    
    
    ! Convert xvar to ipar first
    do k1=1,size(xvar)
        ind = optim(k1)
        ipar(ind) = lin_scale_back(xvar(k1), ipar_max(k1), ipar_min(k1))
    enddo
    
    ! Then convert mpar to ipar:
    !       1,  2,  3,    4,      5,   6,       7,   8,       9
    ! mpar: E, nu, Y0, Hiso, 1/kinf, Hk1, 1/ainf1, Hk2, 1/ainf2, ....
    ! ipar: E, nu, Y0, Htot,   Ysat, rk1,     ra1, rk2,     ra2, ...
    mpar = ipar
    ! Plastic modulii
    mpar(4) = ipar(4)*(1.d0 - sum(ipar(6::2)))  ! Isotropic
    mpar(6::2) = ipar(4)*ipar(6::2)             ! Kinematic
    ! Saturation stresses
    mpar(5) = 1.d0/((ipar(5)-ipar(3))*(1.d0 - sum(ipar(7::2))))
    mpar(7::2) = 1.d0/((ipar(5)-ipar(3))*ipar(7::2))
    

end subroutine
    
subroutine ipar_to_xvar(ipar, xvar, optim, ipar_min, ipar_max)
!DEC$ ATTRIBUTES DLLEXPORT :: ipar_to_xvar
use util_funs
implicit none
    double precision, intent(in)    :: ipar(:)      ! Contains the current material parameters
    double precision, intent(inout) :: xvar(:)      ! Contains scaled parameters [0,1] that are optimized, to be updated
    integer, intent(in)             :: optim(:)     ! Contains the position of the optimized parameters in mpar
    double precision, intent(in)    :: ipar_min(:)  ! Min values for input material parameters
    double precision, intent(in)    :: ipar_max(:)  ! Max values for input material parameters
    
    !Internal variables
    integer                         :: ind, k1
    
    do k1=1,size(xvar)
        ind = optim(k1)
        xvar(k1) = lin_scale(ipar(ind), ipar_max(k1), ipar_min(k1))
    enddo

end subroutine
    
subroutine xvar_to_ipar(ipar, xvar, optim, ipar_min, ipar_max)
!DEC$ ATTRIBUTES DLLEXPORT :: xvar_to_ipar
use util_funs
implicit none
    ! Convert scaled parameters into input material parameters. 
    double precision, intent(inout) :: ipar(:)      ! Contains the initial input parameters, to be updated
    double precision, intent(in)    :: xvar(:)      ! Contains scaled parameters [0,1] that are optimized
    integer, intent(in)             :: optim(:)     ! Contains the position of the optimized parameters in ipar
    double precision, intent(in)    :: ipar_min(:)  ! Min values for input material parameters
    double precision, intent(in)    :: ipar_max(:)  ! Max values for input material parameters
    
    ! Internal variables
    integer                         :: ind, k1
    
    do k1=1,size(xvar)
        ind = optim(k1)
        ipar(ind) = lin_scale_back(xvar(k1), ipar_max(k1), ipar_min(k1))
    enddo
    
end subroutine    