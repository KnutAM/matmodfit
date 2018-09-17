! Convert and scale between "true" (material) parameters (mpar) used for analysis and the scaled variables (xvar) used for optimization
! glob%ipar_sctype controls how the parameters are scaled:
! 0: No scaling
! 1: Linear scaling so that between 0 and 1
! 2: Log scaling, so that between 0 and 1
! 3: Inverse scaling, then linear scaling so that between 0 and 1
!
! To add further scaling options:
! Add in routines true_to_scaled and scaled_to_true
! If this scaling does not produce a range (0,1), then you should treat this is in the 
! get_lower_bounds and get_upper_bounds functions
    
module convert_mpar_mod
    use types_mod
    use output_mod
    use usr_interface_mod
    implicit none
    
    private
    ! xvar is the array of parameters seen by the optimization routines
    ! mpar is the array of parameters seen by the material routine
    ! ipar is the array of parameters given by the input file
    public :: xvar_to_mpar
    public :: ipar_to_xvar
    public :: xvar_to_ipar 
    public :: get_lower_bounds, get_upper_bounds
  
    contains

 
subroutine xvar_to_ipar(ipar, xvar, glob)
    implicit none
    type(glob_typ)      :: glob
    double precision    ::  ipar(:), xvar(:)
    procedure(user_xvar_to_ipar_template),pointer:: user_xvar_to_ipar ! Addresss to user subroutine

    ipar = glob%ipar_init
    if (glob%user_scale) then
        user_xvar_to_ipar => glob%xvar_to_ipar_addr
        call user_xvar_to_ipar(ipar, xvar, glob%ipar_optim, glob%ipar_min, glob%ipar_max)
    else
        ! "un-scale" the optimized material parameters, and put them on the correct place in mpar
        ipar(glob%ipar_optim) = scaled_to_true(xvar, glob)
    endif
    
end subroutine xvar_to_ipar
    
subroutine xvar_to_mpar(mpar, xvar, glob)
    implicit none
    type(glob_typ)      ::  glob
    double precision    ::  mpar(:), xvar(:)
    procedure(user_xvar_to_mpar_template),pointer:: user_xvar_to_mpar ! Addresss to user subroutine

    mpar = glob%ipar_init
    if (glob%user_scale) then
        user_xvar_to_mpar => glob%xvar_to_mpar_addr
        call user_xvar_to_mpar(mpar, glob%ipar_init, xvar, glob%ipar_optim, glob%ipar_min, glob%ipar_max)
    else
        ! "un-scale" the optimized material parameters, and put them on the correct place in mpar
        mpar(glob%ipar_optim) = scaled_to_true(xvar, glob)
    endif
    
    
end subroutine xvar_to_mpar

subroutine ipar_to_xvar(ipar, xvar, glob)
    implicit none
    type(glob_typ)      ::  glob
    double precision    ::  ipar(:), xvar(:)
    procedure(user_ipar_to_xvar_template),pointer:: user_ipar_to_xvar ! Addresss to user subroutine
    
    if (glob%user_scale) then
        user_ipar_to_xvar => glob%ipar_to_xvar_addr
        call user_ipar_to_xvar(ipar, xvar, glob%ipar_optim, glob%ipar_min, glob%ipar_max)
    else
        xvar = true_to_scaled(ipar(glob%ipar_optim), glob)
    endif
    
end subroutine ipar_to_xvar

function true_to_scaled(true_par, glob) result(scal_par)
implicit none
    double precision, parameter     :: invtol = 1.d-30  !Minimum value to allow normal inverse scaling
    double precision    :: true_par(:)
    type (glob_typ)     :: glob
    double precision    :: scal_par(size(true_par))
    integer             :: k1, sctype
    double precision    :: minval, maxval
    
    do k1=1,size(true_par)
        sctype = glob%ipar_sctype(k1)
        
        minval = glob%ipar_min(k1) ! minimum value
        maxval = glob%ipar_max(k1) ! maximum value
        if      (sctype==0) then    ! No scaling
            scal_par(k1) = true_par(k1)
        elseif  (sctype==1) then    ! Linear scaling
            scal_par(k1) = (true_par(k1)-minval)/(maxval-minval)
        elseif  (sctype==2) then    ! Log scaling
            if ((minval<=0.d0).or.(maxval<=0.d0)) then
                call write_output('Parameter nr. '//int2str(k1)//': Range values must greater than zero for log scaling', 'error', 'opt')
            endif
            scal_par(k1) = (log(true_par(k1)) - log(minval)) / (log(maxval) - log(minval))
        elseif  (sctype==3) then    ! Inverse scaling
            if ((abs(minval)<invtol).or.(abs(maxval)<invtol)) then
                call write_output('Parameter nr. '//int2str(k1)//': Range values cannot be too close to zero for inverse scaling', 'error', 'opt', halt=.false.)
                call write_output('Min abs value is '//dbl2str(invtol)//', current min = '//dbl2str(minval)//' and max = '//dbl2str(maxval), 'error', 'opt', loc=.false.)
            endif
            if (abs(true_par(k1))<invtol) then
                call write_output('Parameter nr. '//int2str(k1)//'= '//dbl2str(true_par(k1))//'. Inverse scaling is dangerous', 'error', 'opt', halt=.false.)
                call write_output('Setting parameter to (almost) infinity', 'error', 'opt', loc=.false.)
                scal_par(k1) = sign(1.d0, true_par(k1))*huge(1.d0)
            else
                scal_par(k1) = 1/true_par(k1)
            endif
            scal_par(k1) = (scal_par(k1)-1/maxval)/(1/minval-1/maxval)
        else
            call write_output('Parameter nr. '//int2str(k1)//': ipar_sctype = '//int2str(sctype)//' is not supported', 'error', 'opt')
        endif
        
    enddo
    
end function true_to_scaled

function scaled_to_true(scal_par, glob) result(true_par)
implicit none
    double precision, parameter     :: invtol = 1.d-30  !Minimum value to allow normal inverse scaling
    double precision    :: scal_par(:)    
    type (glob_typ)     :: glob
    double precision    :: true_par(size(scal_par))
    integer             :: k1, sctype
    double precision    :: minval, maxval
    
    
    do k1=1,size(scal_par)
        sctype = glob%ipar_sctype(k1)
        minval = glob%ipar_min(k1) ! minimum value
        maxval = glob%ipar_max(k1) ! maximum value
        if      (sctype==0) then    ! No scaling
            true_par(k1) = scal_par(k1)
        elseif  (sctype==1) then    ! Linear scaling
            true_par(k1) = scal_par(k1)*(maxval-minval) + minval
        elseif  (sctype==2) then    ! Log scaling
            true_par(k1) = minval*(maxval/minval)**(scal_par(k1))
        elseif  (sctype==3) then    ! Inverse scaling
            true_par(k1) = scal_par(k1)*(1/minval-1/maxval) + 1/maxval
            if (abs(true_par(k1))<invtol) then
                call write_output('Parameter nr. '//int2str(k1)//'= '//dbl2str(true_par(k1))//'. Inverse scaling is dangerous', 'error', 'opt', halt=.false.)
                call write_output('Setting parameter to (almost) infinity', 'error', 'opt', loc=.false.)
                true_par(k1) = sign(1.d0, true_par(k1))*huge(1.d0)
            else
                true_par(k1) = 1/true_par(k1)
            endif
        else
            call write_output('Parameter nr. '//int2str(k1)//': ipar_sctype = '//int2str(sctype)//' is not supported', 'error', 'opt')
        endif
    enddo
    
end function scaled_to_true
    
function get_lower_bounds(glob) result(lb)
implicit none
    type(glob_typ)       :: glob
    double precision    :: lb(size(glob%ipar_sctype))
    integer             :: k1
    
    lb = 0.d0   ! Default
    do k1=1,size(lb)
        if (glob%ipar_sctype(k1)==0) then !parameters not scaled, use initial range
            lb(k1) = glob%ipar_min(k1)
        endif
    enddo
end function

function get_upper_bounds(glob) result(ub)
implicit none
    type(glob_typ)       :: glob
    double precision    :: ub(size(glob%ipar_sctype))
    integer             :: k1
    
    ub = 1.d0   ! Default
    do k1=1,size(ub)
        if (glob%ipar_sctype(k1)==0) then !parameters not scaled, use initial range
            ub(k1) = glob%ipar_max(k1)
        endif
    enddo
end function

end module convert_mpar_mod
