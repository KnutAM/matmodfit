!DEC$ FREEFORM
module StVernant_tensor_module
implicit none

    contains
!   This module uses the following voigt tensor notation:
!   T = [T_11, T_22, T_33, T_12, T_23, T_31, T_13, T_21, T_32]
!   and contains the following routines:
!
!   eye_v9                  Make a 9-voigt identity tensor
!   trans_v9                Take the transpose of 9-voigt tensor
!   tr_v9                   Take the trace of a 9-voigt tensor
!   det_v9                  Take the determinant of a 9-voigt tensor
!   m_2_v9                  Convert a 3x3 matrix to a 9-voigt tensor
!   v9_d_v9                 Single contraction (2nd order): c_ik = a_ij b_jk
!   v9x9_dd_v9x9            double contraction (4th order): c_ijmn = a_ijkl b_klmn
!   op_v9                   Open product (2nd order):       a_ijkl = b_ij c_kl
!   op_a_v9                 Open product above (2nd order): a_ijkl = b_ik c_jl
!   op_b_v9                 Open product below (2nd order): a_ijkl = b_il c_jk
!   minsym_v9x9             Make 4th order minor symmetric
!   dSdE_2_dtaudF           Convert dS/dE to d(tau)/dF
!   dtaudF_2_ddsdde         Convert d(tau)/dF to abaqus ddsdde
!   sigma_v9_2_sigma_abaqus Convert sigma_v9 in voigt format to sigma with 6 components in abaqus format

    !==============================================================================     
    function eye_v9() result(eye)
       implicit none
       double precision             :: eye(9)

       eye = 0.0d0
       eye(1:3) = (/ 1.0d0, 1.0d0, 1.0d0 /)

    end function eye_v9
    !==============================================================================  
    
    !============================================================================== 
    function trans_v9(b) result (a)
    implicit none

        double precision, intent (in) :: b(9)
        double precision :: a(9)

        a(1:3)=b(1:3)
        a((/4,5,6,7,8,9/))=b((/8,9,7,6,4,5/))

    end function
    !==============================================================================

    !============================================================================== 
    function tr_v9(a) result(tra)
    implicit none

          double precision, intent (in) :: a(9)
          double precision :: tra
          tra=sum(a(1:3)); 
      
    end function 
    !==============================================================================
        
    !==============================================================================       
    function det_v9(A) result(detA)
    implicit none

        double precision, intent (in) :: A(9)
        double precision :: detA
        detA=-A(7)*A(2)*A(6)+A(4)*A(5)*A(6)+  &
              A(7)*A(8)*A(9)-A(1)*A(5)*A(9)-  &
              A(4)*A(8)*A(3)+A(1)*A(2)*A(3)  
           
    end function 
    !==============================================================================

    !==============================================================================
    function m_2_v9(b) result (a)
    implicit none
    
        double precision, intent (in) :: b(3,3)
        double precision :: a(9)
      
        a(1)=b(1,1)
        a(2)=b(2,2)
        a(3)=b(3,3)
        a(4)=b(1,2)
        a(9)=b(3,2)
        a(5)=b(2,3)
        a(8)=b(2,1)
        a(7)=b(1,3)
        a(6)=b(3,1)

    end function 
    !==============================================================================
    
    !============================================================================== 
    function v9_d_v9(a,b) result (c)
    implicit none

        integer i,j
        double precision, intent (in) :: a(9),b(9)
        double precision :: c(9)

        c(1)=a(1)*b(1) + a(7)* b(6) + a(4)* b(8)
        c(2)=a(2)*b(2) + a(8)* b(4) + a(5)* b(9)
        c(3)=a(3)*b(3) + a(9)* b(5) + a(6)* b(7)
        c(4)=a(4)*b(2) + a(1)* b(4) + a(7)* b(9)
        c(5)=a(5)*b(3) + a(2)* b(5) + a(8)* b(7)
        c(6)=a(6)*b(1) + a(3)* b(6) + a(9)* b(8)
        c(7)=a(7)*b(3) + a(4)* b(5) + a(1)* b(7)
        c(8)=a(8)*b(1) + a(5)* b(6) + a(2)* b(8)
        c(9)=a(9)*b(2) + a(6)* b(4) + a(3)* b(9)

    end function
    !============================================================================== 

    !==============================================================================
    function v9x9_dd_v9x9(a,b) result (c)
    implicit none
          double precision, intent (in) :: a(9,9),b(9,9)
          double precision :: c(9,9)
          integer I,J,K
          call dgemm('N','N', 9, 9, 9, 1.d0,a, 9, b, 9, 0.d0, c, 9)
          !c = matmul(a,b)  ! Optional if not using mkl
    end function
    !============================================================================== 

    !==============================================================================  
    function op_v9(b,c) result (a)
    implicit none

        double precision, intent (in) :: b(9),c(9)
        double precision :: a(9,9)
        integer I,J
        do I=1,9
           do J=1,9
                 a(I,J)=b(I)*c(J)
           enddo
        enddo

    end function 
    !==============================================================================  

    !==============================================================================   
    function op_a_v9(b,c) result (a) 
    implicit none

        double precision, intent (in) :: b(9),c(9)
        double precision :: a(9,9)
    
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
    !==============================================================================

    !============================================================================== 
    function op_b_v9(b,c) result (a)
    implicit none

        double precision, intent (in) :: b(9),c(9)
        double precision :: a(9,9)
    
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
    !==============================================================================
    
    !==============================================================================
    function minsym_v9x9(b) result (a)
    implicit none

        double precision, intent (in) :: b(9,9)
        double precision :: a(9,9)
        integer:: upper(3), lower(3)
        upper = (/4, 5, 7/) !Indicies on upper diagonal
        lower = (/8, 9, 6/) !Indicies on lower diagonal
        a = b
    
        a(:, upper) = (a(:, upper) + a(:, lower))/2.d0  !Symmetrize columns
        a(:, lower) = a(:, upper)
    
        a(upper, :) = (a(upper, :) + a(lower, :))/2.d0  !Symmetrize rows
        a(lower, :) = a(upper,:)

    end function
    !==============================================================================

    !==============================================================================     
    function dSdE_2_dtaudF(dSdE,F,S) result(dtaudF)
       implicit none
       double precision, intent(in) :: dSdE(9,9), F(9), S(9)
       double precision             :: dtaudF(9,9), I2(9), Ft(9)
        I2 = eye_v9()
        Ft = trans_v9(F)
        dtaudF = op_a_v9(I2, v9_d_v9(F, S)) + op_b_v9(v9_d_v9(F, S), I2) &
            + 0.5*v9x9_dd_v9x9(v9x9_dd_v9x9(op_a_v9(F, F), dSdE), op_b_v9(I2, Ft)+op_a_v9(Ft, I2))

    end function dSdE_2_dtaudF
    !==============================================================================  

    !==============================================================================
    function dtaudF_2_ddsdde(dtaudF, F) result(ddsdde)
        implicit none
        double precision, intent(in) :: dtaudF(9,9), F(9)
        double precision             :: tmp(9,9), ddsdde(6,6)
        integer                      :: posvec(6)
    
        tmp = (1.d0/det_v9(F)) * v9x9_dd_v9x9(dtaudF, op_a_v9(eye_v9(), trans_v9(F)))
        tmp = minsym_v9x9(tmp)
    
        posvec = (/1, 2, 3, 4, 6, 5/)   ! Index conversion
        ddsdde = tmp(posvec, posvec)
    
    end function dtaudF_2_ddsdde
    !==============================================================================

    !==============================================================================
    function sigma_v9_2_sigma_abaqus(sigma_v9) result(sigma_abaqus)
        !Convert sigma in tensor module format to abaqus output format
       implicit none
       double precision, intent(in) :: sigma_v9(9)
       double precision             :: sigma_abaqus(6)
       double precision             :: tmp(9)
        tmp = 0.5*(sigma_v9 + trans_v9(sigma_v9))   !Symmetrize
        sigma_abaqus = tmp((/1, 2, 3, 4, 6, 5/))      !Convert indicies

    end function sigma_v9_2_sigma_abaqus
    !==============================================================================

end module StVernant_tensor_module
!   St. Verntant type of elasticity example umat for mpfit
!   Note the "!DEC$ ATTRIBUTES DLLEXPORT :: UMAT" directly following the subroutine declaration, 
!   ensuring that the subroutine is exported when making a dll using intel ifort.
!   The material routine takes two material parameters: mu (also known as shear modulus G) and lambda
!   These are the Lamé constants for an elastic solid.
subroutine umat(stress,statev,ddsdde,sse,spd,scd, &
    rpl,ddsddt,drplde,drpldt,&
    stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
    ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
    celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!DEC$ ATTRIBUTES DLLEXPORT :: UMAT
	use StVernant_tensor_module
    implicit none
        
    character(len=80) cmname
    integer:: ndi,nshr,ntens,nstatv,nprops,noel,npt,&
    layer, kspt, kstep, kinc
    double precision:: dtime,temp,dtemp,pnewdt,celent,sse,spd,scd,rpl,drpldt
    double precision:: stress(ntens),statev(nstatv),&
    ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),&
    stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),&
    props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

	double precision :: mu, la, S9(9), E9(9), I2(9), F9(9), C9(9), dSdE(9,9), dtau_dF(9,9), stress9(9)
	! Extract properties (Lamé constants)
    mu = props(1)
	la = props(2)
		
    ! Calculate the initial strain variables
	I2 = eye_v9()
	F9 = m_2_v9(dfgrd1)
	C9 = m_2_v9(matmul(transpose(dfgrd1), dfgrd1))
	E9 = (C9-I2)/2.d0
    
    
    ! Calculate 2nd Piola-Kirchoff stress
	S9 = la*tr_v9(E9)*I2 + 2*mu*E9
    
    ! Convert to Cauchy stress
	stress9 = v9_d_v9(v9_d_v9(F9, S9), trans_v9(F9))/det_v9(F9)
	
    ! Convert to abaqus' format for stress
	stress = sigma_v9_2_sigma_abaqus(stress9)
		
    ! Calculate tangent stiffness for 2nd Piola-Kirchoff wrt Green-Lagrange strain
	dSdE = la*op_v9(I2,I2) + mu*(op_a_v9(I2,I2) + op_b_v9(I2,I2))
		
    ! Convert to d(tau)/dF
	dtau_dF = dSdE_2_dtaudF(dSdE,F9,S9)
		
    ! Convert to Abaqus' tangent stiffness ddsdde
	ddsdde = dtaudF_2_ddsdde(dSdE, F9)
		
end subroutine