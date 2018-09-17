module sim_tensor_mod
use output_mod
implicit none

private ! Private by default
public :: dPdF_calc, P9_2_sigma6
public :: scont22, scont42, dcont44, dcont42
public :: openprod, openprod_a, openprod_b
public :: invert_v9, determinant_v9, transpose_v9
public :: voigt9, voigt4, matrix, matrix4
public :: stran_v9, stress_v9, stran_v6, stress_v6
public :: v9x9, v6x6, minsym9x9, m1_m3_mul

!Tensor modulue convention used for index notation
integer, parameter, dimension(3)  :: norm = (/1,2,3/)         !Normal components
integer, parameter, dimension(3)  :: shr1 = (/4,7,5/)         !Upper shear components
integer, parameter, dimension(3)  :: shr2 = (/8,6,9/)         !Lower shear components
integer, parameter, dimension(6)  :: first= (/1,2,3,4,7,5/)   !Normal + upper shear components
integer, parameter, dimension(9)  :: ind1 = (/1, 2, 3, 1, 2, 3, 1, 2, 3/)
integer, parameter, dimension(9)  :: ind2 = (/1, 2, 3, 2, 3, 1, 3, 1, 2/)

!Standard tensors
double precision, parameter, dimension(9)   :: I2 = (/1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/)

    contains
    
!! ==================== Tensor functions ====================

!Given ddsdde, calculate P9 and dPdF 
subroutine dPdF_calc(P9, dPdF, stress, ddsdde, F9)
implicit none
    double precision, intent(in)    :: stress(6), ddsdde(6,6), F9(9)
    double precision, intent(out)   :: P9(9), dPdF(9,9)
    double precision                :: Finv(9), Finvt(9), detF, tmp9(9), tmp99a(9,9), tmp99b(9,9)
    
    Finv    = invert_v9(F9)
    Finvt   = transpose_v9(Finv)
    detF    = determinant_v9(F9)
    tmp9    = stress_v9(stress) ! 9-component sigma
    P9      = detF*scont22(tmp9, Finvt)
    !P9      = detF*scont22(stress_v9(stress), Finvt)
    tmp99a  = openprod_a(I2,Finv)
    tmp99b  = v9x9(ddsdde)
    tmp99a  = dcont44(tmp99a, tmp99b)
    tmp99b  = openprod_a(I2, Finvt)
    dPdF    = detF*dcont44(tmp99a, tmp99b) - openprod_b(P9, Finv)
    !dPdF    = detF*dcont44(dcont44(openprod_a(I2,Finv), v9x9(ddsdde)), openprod_a(I2, Finvt)) - openprod_b(P9, Finv)
end subroutine dPdF_calc

!Given P9 and F9, calculate sigma6
function P9_2_sigma6(P9, F9) result(sigma6)
implicit none
    double precision    :: P9(9), F9(9), sigma6(6)
    double precision    :: tmp9(9)
    tmp9 = transpose_v9(F9)     ! F transpose
    tmp9 = scont22(P9, tmp9)    ! tau
    sigma6 = stress_v6(tmp9)/determinant_v9(F9)
end function

!Single contraction between 2 2nd order 9-voigt tensors
function scont22(a2, b2)
implicit none
    double precision :: a2(9), b2(9), scont22(9)
    double precision :: am(3,3), bm(3,3)
    am = matrix(a2)
    bm = matrix(b2)
    am = matmul(am,bm)
    scont22 = voigt9(am)
    !scont22 = voigt9(matmul(matrix(a2), matrix(b2)))
end function scont22

!Single contraction between a 4th order tensor [9x9] (last index) and 2nd order tensor [9] (first index)
function scont42(a4, b2)
implicit none
    double precision :: a4(9,9), b2(9), scont42(9,9)
    double precision :: aM(3,3,3,3), bM(3,3), scontM(3,3,3,3)
    integer :: k1,k2,k3,k4
    
    aM = matrix4(a4)
    bM = matrix(b2)
    do k1=1,3
        do k2=1,3
            do k3=1,3
                do k4=1,3
                    scontM(k1,k2,k3,k4) = sum(aM(k1,k2,k3,:)*bM(:,k4))
                enddo
            enddo
        enddo
    enddo
    scont42 = voigt4(scontM)
end function

!Double contraction between two 4th order tensors
function dcont44(a4,b4) result(c4)
implicit none
    double precision :: a4(9,9), b4(9,9), c4(9,9)
    c4 = matmul(a4,b4)
end function

!Double contraction between a 4th and a 2nd order tensor
function dcont42(a4,b2) result(c2)
implicit none
    double precision :: a4(9,9), b2(9), c2(9)
    c2 = matmul(a4,b2)
end function

!Open product between 2 2nd order tensors [9]
function openprod(a2,b2) result(ab)
implicit none
    double precision :: a2(9), b2(9), ab(9,9)
    integer :: k1, k2
    ! This could be replaced by matmul?
    do k1=1,9
        do k2 = 1,9
            ab(k1,k2) = a2(k1)*b2(k2)
        enddo
    enddo
end function openprod

!Open product "above" between 2 2nd order tensors
function openprod_a(b,c) result(a)
implicit none
    double precision :: b(9), c(9), a(9,9)
    ! Copied from tensor module, op_a_V9
    ! Derived from matlab program derive_expr
    a(1,1)=b(1)*c(1); a(1,2)=b(4)*c(4); a(1,3)=b(7)*c(7);
    a(1,4)=b(1)*c(4); a(1,5)=b(4)*c(7); a(1,6)=b(7)*c(1);
    a(1,7)=b(1)*c(7); a(1,8)=b(4)*c(1); a(1,9)=b(7)*c(4);
    a(2,1)=b(8)*c(8); a(2,2)=b(2)*c(2); a(2,3)=b(5)*c(5);
    a(2,4)=b(8)*c(2); a(2,5)=b(2)*c(5); a(2,6)=b(5)*c(8);
    a(2,7)=b(8)*c(5); a(2,8)=b(2)*c(8); a(2,9)=b(5)*c(2);
    a(3,1)=b(6)*c(6); a(3,2)=b(9)*c(9); a(3,3)=b(3)*c(3);
    a(3,4)=b(6)*c(9); a(3,5)=b(9)*c(3); a(3,6)=b(3)*c(6);
    a(3,7)=b(6)*c(3); a(3,8)=b(9)*c(6); a(3,9)=b(3)*c(9);
    a(4,1)=b(1)*c(8); a(4,2)=b(4)*c(2); a(4,3)=b(7)*c(5);
    a(4,4)=b(1)*c(2); a(4,5)=b(4)*c(5); a(4,6)=b(7)*c(8);
    a(4,7)=b(1)*c(5); a(4,8)=b(4)*c(8); a(4,9)=b(7)*c(2);
    a(5,1)=b(8)*c(6); a(5,2)=b(2)*c(9); a(5,3)=b(5)*c(3);
    a(5,4)=b(8)*c(9); a(5,5)=b(2)*c(3); a(5,6)=b(5)*c(6);
    a(5,7)=b(8)*c(3); a(5,8)=b(2)*c(6); a(5,9)=b(5)*c(9);
    a(6,1)=b(6)*c(1); a(6,2)=b(9)*c(4); a(6,3)=b(3)*c(7);
    a(6,4)=b(6)*c(4); a(6,5)=b(9)*c(7); a(6,6)=b(3)*c(1);
    a(6,7)=b(6)*c(7); a(6,8)=b(9)*c(1); a(6,9)=b(3)*c(4);
    a(7,1)=b(1)*c(6); a(7,2)=b(4)*c(9); a(7,3)=b(7)*c(3);
    a(7,4)=b(1)*c(9); a(7,5)=b(4)*c(3); a(7,6)=b(7)*c(6);
    a(7,7)=b(1)*c(3); a(7,8)=b(4)*c(6); a(7,9)=b(7)*c(9);
    a(8,1)=b(8)*c(1); a(8,2)=b(2)*c(4); a(8,3)=b(5)*c(7);
    a(8,4)=b(8)*c(4); a(8,5)=b(2)*c(7); a(8,6)=b(5)*c(1);
    a(8,7)=b(8)*c(7); a(8,8)=b(2)*c(1); a(8,9)=b(5)*c(4);
    a(9,1)=b(6)*c(8); a(9,2)=b(9)*c(2); a(9,3)=b(3)*c(5);
    a(9,4)=b(6)*c(2); a(9,5)=b(9)*c(5); a(9,6)=b(3)*c(8);
    a(9,7)=b(6)*c(5); a(9,8)=b(9)*c(8); a(9,9)=b(3)*c(2);
end function

!Open product "below" between 2 2nd order tensors
function openprod_b(b,c) result(a)
implicit none
    double precision :: b(9), c(9), a(9,9)
    ! Copied from tensor module, op_b_V9
    ! Derived from matlab program derive_expr
    a(1,1)=b(1)*c(1);a(1,2)=b(4)*c(4);a(1,3)=b(7)*c(7);
    a(1,4)=b(4)*c(1);a(1,5)=b(7)*c(4);a(1,6)=b(1)*c(7);
    a(1,7)=b(7)*c(1);a(1,8)=b(1)*c(4);a(1,9)=b(4)*c(7);
    a(2,1)=b(8)*c(8);a(2,2)=b(2)*c(2);a(2,3)=b(5)*c(5);
    a(2,4)=b(2)*c(8);a(2,5)=b(5)*c(2);a(2,6)=b(8)*c(5);
    a(2,7)=b(5)*c(8);a(2,8)=b(8)*c(2);a(2,9)=b(2)*c(5);
    a(3,1)=b(6)*c(6);a(3,2)=b(9)*c(9);a(3,3)=b(3)*c(3);
    a(3,4)=b(9)*c(6);a(3,5)=b(3)*c(9);a(3,6)=b(6)*c(3);
    a(3,7)=b(3)*c(6);a(3,8)=b(6)*c(9);a(3,9)=b(9)*c(3);
    a(4,1)=b(1)*c(8);a(4,2)=b(4)*c(2);a(4,3)=b(7)*c(5);
    a(4,4)=b(4)*c(8);a(4,5)=b(7)*c(2);a(4,6)=b(1)*c(5);
    a(4,7)=b(7)*c(8);a(4,8)=b(1)*c(2);a(4,9)=b(4)*c(5);
    a(5,1)=b(8)*c(6);a(5,2)=b(2)*c(9);a(5,3)=b(5)*c(3);
    a(5,4)=b(2)*c(6);a(5,5)=b(5)*c(9);a(5,6)=b(8)*c(3);
    a(5,7)=b(5)*c(6);a(5,8)=b(8)*c(9);a(5,9)=b(2)*c(3);
    a(6,1)=b(6)*c(1);a(6,2)=b(9)*c(4);a(6,3)=b(3)*c(7);
    a(6,4)=b(9)*c(1);a(6,5)=b(3)*c(4);a(6,6)=b(6)*c(7);
    a(6,7)=b(3)*c(1);a(6,8)=b(6)*c(4);a(6,9)=b(9)*c(7);
    a(7,1)=b(1)*c(6);a(7,2)=b(4)*c(9);a(7,3)=b(7)*c(3);
    a(7,4)=b(4)*c(6);a(7,5)=b(7)*c(9);a(7,6)=b(1)*c(3);
    a(7,7)=b(7)*c(6);a(7,8)=b(1)*c(9);a(7,9)=b(4)*c(3);
    a(8,1)=b(8)*c(1);a(8,2)=b(2)*c(4);a(8,3)=b(5)*c(7);
    a(8,4)=b(2)*c(1);a(8,5)=b(5)*c(4);a(8,6)=b(8)*c(7);
    a(8,7)=b(5)*c(1);a(8,8)=b(8)*c(4);a(8,9)=b(2)*c(7);
    a(9,1)=b(6)*c(8);a(9,2)=b(9)*c(2);a(9,3)=b(3)*c(5);
    a(9,4)=b(9)*c(8);a(9,5)=b(3)*c(2);a(9,6)=b(6)*c(5);
    a(9,7)=b(3)*c(8);a(9,8)=b(6)*c(2);a(9,9)=b(9)*c(5);
end function

!Invert a 9 component voigt 2nd order tensor
function invert_v9(A) result (inv_A)
implicit none
    double precision :: A(9), inv_A(9),det_A
    !Code mostly copied from Tensor module, inv_v9
    det_A = determinant_v9(A)
    
    if (abs(det_A)<1.0d-25) then
        call write_output('Matrix close to being singular in inv_v9', 'warning', 'atp')
        det_A=1.d-15
    endif
    
    inv_A(1)=(-A(5)*A(9)+A(2)*A(3))/det_A
    inv_A(4)=( A(7)*A(9)-A(4)*A(3))/det_A
    inv_A(7)=( A(5)*A(4)-A(2)*A(7))/det_A

    inv_A(8)=( A(5)*A(6)-A(8)*A(3))/det_A
    inv_A(2)=(-A(7)*A(6)+A(1)*A(3))/det_A
    inv_A(5)=(-A(5)*A(1)+A(8)*A(7))/det_A

    inv_A(6)=(-A(2)*A(6)+A(8)*A(9))/det_A
    inv_A(9)=( A(4)*A(6)-A(1)*A(9))/det_A
    inv_A(3)=( A(2)*A(1)-A(8)*A(4))/det_A
    
end function invert_v9

!Calculate determinant of 9-component voigt 2nd order tensor
function determinant_v9(A) result (detA)
implicit none
    double precision :: A(9), detA
    detA=-A(7)*A(2)*A(6)+A(4)*A(5)*A(6)+  &
          A(7)*A(8)*A(9)-A(1)*A(5)*A(9)-  &
          A(4)*A(8)*A(3)+A(1)*A(2)*A(3)
end function determinant_v9

!Take transpose of 9-component voigt 2nd order tensor
function transpose_v9(v9) result (v9t)
implicit none
    double precision :: v9(9), v9t(9)
    v9t(norm) = v9(norm)
    v9t(shr1) = v9(shr2)
    v9t(shr2) = v9(shr1)
end function transpose_v9


!! ==================== Conversion functions ====================
!Convert 3x3 matrix to 9 voigt
function voigt9(matrix)
implicit none
    double precision :: voigt9(9), matrix(3,3)
    integer :: k1
    do k1=1,9
        voigt9(k1) = matrix(ind1(k1), ind2(k1))
    enddo
    
end function voigt9

!Convert 3x3x3x3 matrix to 9x9 voigt (4th order tensor)
function voigt4(matrix)
implicit none
    double precision :: voigt4(9,9), matrix(3,3,3,3)
    integer          :: k1, k2
    
    do k1 = 1,9
        do k2 = 1,9
            voigt4(k1,k2) = matrix(ind1(k1),ind2(k1),ind1(k2),ind2(k2))
        enddo
    enddo
end function voigt4

!Convert 9 voigt to 3x3 matrix
function matrix(voigt9)
implicit none
    double precision :: matrix(3,3), voigt9(9)
    integer          :: k1
    do k1 = 1,9
        matrix(ind1(k1), ind2(k1)) = voigt9(k1)
    enddo
    
end function matrix

!Convert 9x9 voigt to 3x3x3x3 matrix
function matrix4(voigt4)
implicit none
    double precision :: voigt4(9,9), matrix4(3,3,3,3)
    integer          :: k1, k2
    
    do k1 = 1,9
        do k2 = 1,9
            matrix4(ind1(k1),ind2(k1),ind1(k2),ind2(k2)) = voigt4(k1,k2)
        enddo
    enddo

end function matrix4

!Convert strain in 6 voigt (Abaqus convention) to strain in 9 voigt (Tensor module convention)
function stran_v9(v6) result(v9)
implicit none
    double precision :: v9(9), v6(6)
    v9(norm) = v6(1:3)
    v9(shr1) = v6(4:6)/2.d0
    v9(shr2) = v6(4:6)/2.d0
end function

!Convert strain in 9 voigt (Tensor module convention) to strain in 6 voigt (Abaqus convention)
function stran_v6(v9) result(v6)
implicit none
    double precision :: v9(9), v6(6)
    v6(1:3) = v9(norm)
    v6(4:6) = v9(shr1)+v9(shr2)
end function

!Convert stress in 6 voigt (Abaqus convention) to stress in 9 voigt (Tensor module convention) 
function stress_v9(v6) result(v9)
implicit none
    double precision :: v9(9), v6(6)
    v9(first) = v6
    v9(shr2)  = v6(4:6)
end function

!Convert stress in 9 voigt (Tensor module convention) to stress in 6 voigt (Abaqus convention)
function stress_v6(v9) result(v6)
implicit none
    double precision :: v6(6), v9(9)
    v6(1:3) = v9(norm)
    v6(4:6) = (v9(shr1) + v9(shr2))/2.d0
end function

!Convert 6x6 voigt to 9x9 voigt (Abaqus to Tensor module index convention)
!Assume dsde with double shear components for e6
function v9x9(v6) result(v9)
implicit none
    double precision v9(9,9), v6(6,6)
    
    v9(norm,norm) = v6(1:3, 1:3)
    v9(shr1,norm) = v6(4:6, 1:3)
    v9(shr2,norm) = v6(4:6, 1:3)
    
    v9(norm,shr1) = v6(1:3, 4:6)
    v9(shr1,shr1) = v6(4:6, 4:6)
    v9(shr2,shr1) = v6(4:6, 4:6)
    
    v9(norm,shr2) = v6(1:3, 4:6)
    v9(shr1,shr2) = v6(4:6, 4:6)
    v9(shr2,shr2) = v6(4:6, 4:6)
    
end function

!Convert 9x9 voigt to 6x6 voigt (Tensor module to abaqus index convention)
!Assume dsde with double shear components for e6
function v6x6(v9) result(v6)
implicit none
    double precision v6(6,6), v9(9,9)
    v9 = minsym9x9(v9)
    v6 = v9(first,first)
end function v6x6

!Make a 9x9 voigt minor symmetric
function minsym9x9(b) result (a)
implicit none
    double precision, intent(in) :: b(9,9)
    double precision:: a(9,9)
    a = b
    
    a(:, shr1) = (a(:, shr1) + a(:, shr2))/2.d0  !Symmetrize columns
    a(:, shr2) = a(:, shr1)
    
    a(shr1, :) = (a(shr1, :) + a(shr2, :))/2.d0  !Symmetrize rows
    a(shr2, :) = a(shr1,:)

end function

!! ==================== matrix functions ====================
function m1_m3_mul(m1,m3) result(mout)
implicit none
    double precision :: m1(:), m3(:,:,:)
    double precision :: mout(size(m3,dim=2), size(m3,dim=3))
    integer          :: k1, k2
    do k1=1,size(m3,dim=2)
        do k2=1,size(m3,dim=3)
            mout(k1,k2) = sum(m1*m3(:,k1,k2))
        enddo
    enddo
end function m1_m3_mul

end module
