module sim_util_mod
use sim_tensor_mod
implicit none

    private
    public  solve_mateq  ! Solve constrained equation C(fd,fd)*du(fd)+R(fd)=0 for du(fd), fd=free_dofs
    public  GL_stran_calc
    
    !Identity tensor
    double precision, parameter, dimension(9)   :: I2 = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)
    
    contains
    
    !Solve one step in newton iteration
subroutine solve_mateq(C, R, du, lconv, free_dofs, known_dofs)
    ! Solve C*du+R=0 => C(known_dofs,free_dofs)*du(free_dofs) = -R(known_dofs) subjected to bc governed by ctrl
    ! dgesv documentation: https://software.intel.com/en-us/node/468876
    implicit none
        double precision, intent(in)    :: C(:,:)           !Stiffness
        double precision, intent(in)    :: R(:)             !Residual at current F
        double precision, intent(out)   :: du(size(R))      !Suggested step
        logical, intent(out)            :: lconv            !Tells if solution succeeded
        integer                         :: info             !Information about solution
        integer                         :: nvar             !Number of variables
        integer, intent(in)             :: free_dofs(:)
        integer, intent(in)             :: known_dofs(:)
        integer                         :: ipiv(size(free_dofs))    !pivot indices for LU decomposition
        double precision                :: tmp_b(size(free_dofs)), tmp_C(size(free_dofs), size(free_dofs))

        du = 0.d0
        !du(free_dofs) = -R(free_dofs)    !Input b in dgesv is overwritten to answer x
        nvar = size(free_dofs)
        tmp_C = C(known_dofs, free_dofs)
        tmp_b = -R(known_dofs)
        !    dgesv( n,   nrhs,        a,  lda, ipiv,  b,  ldb, info )
        call dgesv(nvar,    1, tmp_C , nvar, ipiv, tmp_b, nvar, info)
        lconv = (info==0)
        du(free_dofs) = tmp_b

end subroutine solve_mateq

subroutine GL_stran_calc(F9_old, F9_new, stran, dstran)
implicit none
    double precision :: F9_old(9), F9_new(9), stran(6), dstran(6)
    double precision :: E9_old(9), dE9(9)
    E9_old = 0.5d0*(F9_old+transpose_v9(F9_old)) - I2
    dE9 = 0.5d0*(F9_new+transpose_v9(F9_new)) - I2 - E9_old
    stran = stran_v6(E9_old)
    dstran = stran_v6(dE9)
end subroutine

end module