module atp_element_mod
use sim_tensor_mod
use sim_util_mod
use output_mod
implicit none

private 

public      :: element_setup    ! Setup the element type
public      :: element_nlgeom   ! Element for nlgeom=true (non-linear geometric analysis)
public      :: element_lingeom  ! Element for nlgeom=false (linear geometric analysis)
public      :: get_gpinfo       ! Get information about current gauss point
public      :: shapefun         ! Get shape function variables at given radius
public      :: FBeps_calc       ! Calculate deformation gradient and strain at given radius. 
                                ! Deformation gradient is calculated consistently (i.e. no small strain approximations here!)


!Parameters
double precision, parameter :: pi = acos(-1.d0)  !Define pi

! General element properties
integer         :: ngp      !Number of gauss points
integer         :: nnod     !Number of nodes
integer         :: ndof     !Number of dofs (nnod+2)
logical         :: bbar     !Abaqus b-bar method (Selectively reduced integration)

    contains

! Setup of elements (needs to run before using the other subroutines in the module)
subroutine element_setup(set_ngp, set_nnod, set_bbar)
implicit none
    integer :: set_ngp
    integer :: set_nnod
    logical :: set_bbar
    
    ngp = set_ngp
    nnod = set_nnod
    ndof = nnod + 2
    bbar = set_bbar
        
end subroutine element_setup

! Main element routine (Calculate Ke and Re, as well as gpstress, gpstatev and gpF)
subroutine element_nlgeom( Ke, Re, ue, gpstress, gpstatev, gpF, Rpos, H0, &
kinc, kstep, noel, pnewdt, props, cmname, dtemp, temp, time, dt, umat)
use umat_mod
implicit none
    !other input variables
    double precision:: Ke(ndof,ndof), Re(ndof), ue(ndof), gpstress(6,ngp), gpstatev(:,:), gpF(9,ngp), Rpos(nnod), H0
    double precision:: Ke_gp(ndof,ndof), Re_gp(ndof), time(2)
    !umat variables
    procedure(umat_template),pointer:: umat     ! Addresss to umat subroutine
    character(len=80)           :: cmname
    integer                     :: nstatv, nprops, noel, kstep, kinc
    double precision            :: dt,temp,dtemp,pnewdt,celent,sse,spd,scd,rpl,drpldt
    double precision            :: stress(6), ddsdde(6,6), ddsddt(6), drplde(6)
    double precision            :: stran(6), dstran(6), abatime(2), props(:), coords(3), drot(3,3), dfgrd0(3,3)
    double precision            :: dfgrd1(3,3), predef(1), dpred(1)
    double precision, allocatable :: statev(:)
    !internal variables
    double precision:: dNpdR(nnod), Np(nnod), phi, uz, r, R0, w, F9(9), drdR, ue_r(ndof-2), B_ki(9,ndof), C_kij(9,ndof,ndof)
    integer         :: k1
    double precision:: bbar_F22, bbar_B_ki_row2(ndof), bbar_C_kij_row2(ndof,ndof)
    
    allocate(statev(size(gpstatev, dim=1)))
    
    celent = Rpos(nnod)-Rpos(1)
    nstatv = size(gpstatev, dim=1)
    nprops = size(props)
    
    ! Convert to abaqus time definition, as it requires the time in the beginning of the increment.
    abatime = time - dt
    
    sse = 0.d0; spd = 0.d0; scd = 0.d0; rpl = 0.d0
    predef = 0.d0; dpred  = 0.d0
    ddsddt = 0.d0; drplde = 0.d0; drpldt = 0.d0
    
    drot = 0.d0;        drot(1,1) = 1.d0
    drot(2,2) = 1.d0;   drot(3,3) = 1.d0
    
    uz  = ue(1); phi = ue(2); ue_r= ue(3:ndof)    ! Radial displacements
    
    Ke = 0.d0
    Re = 0.d0
    
    if (bbar) then
        R0 = (Rpos(1)+Rpos(size(Rpos)))/2.d0    ! Middle point
        call shapefun(R0, Rpos, dNpdR, Np, ue_r, r, drdR)
        call FBC_calc(F9, B_ki, C_kij, dNpdR, Np, drdR, r, phi, uz, R0, H0)
        bbar_F22 = F9(2)
        bbar_B_ki_row2 = B_ki(2,:)
        bbar_C_kij_row2 = C_kij(2,:,:)
    endif
    
    do k1=1,ngp
        dfgrd0 = matrix(gpF(:,k1))
        call get_gpinfo(R0, w, k1, Rpos)
        call shapefun(R0, Rpos, dNpdR, Np, ue_r, r, drdR)
        call FBC_calc(F9, B_ki, C_kij, dNpdR, Np, drdR, r, phi, uz, R0, H0)
        !call FBC_calc2(F9, dNpdR, Np, drdR, r, phi, uz, R0, H0)
        
        if (bbar) then
            F9(2) = bbar_F22
            B_ki(2,:) = bbar_B_ki_row2
            C_kij(2,:,:) = bbar_C_kij_row2
        endif
        
        dfgrd1 = matrix(F9)
        stress = gpstress(:,k1)
        statev = gpstatev(:,k1)
        coords = (/R0, 0.d0, 0.d0/)
        call GL_stran_calc(gpF(:,k1), F9, stran, dstran)
        call umat(stress,statev,ddsdde,sse,spd,scd, & !Can call umat_internal instead for a test routine without external umat
                rpl,ddsddt,drplde,drpldt,&
                stran,dstran,abatime,dt,temp,dtemp,predef,dpred,cmname,&
                3,3,6,nstatv,props,nprops,coords,drot,pnewdt,&
                celent,dfgrd0,dfgrd1,noel,k1,1,1,kstep,kinc)
        !call Re_calc2(Ke_gp, Re_gp, F9, stress, ddsdde, Np, dNpdR, r, phi, H0, R0) !Calculate Ke_gp and Re_gp
        call Re_calc(Ke_gp, Re_gp, F9, B_ki, C_kij, stress, ddsdde) !Calculate Ke_gp and Re_gp
        Ke = Ke + 2.d0*pi*R0*H0*w*Ke_gp
        Re = Re + 2.d0*pi*R0*H0*w*Re_gp
        
        ! Add updated values to gp-saved variables (the old values of these are saved in the main program)
        gpF(:,k1)       = F9
        gpstress(:,k1)  = stress
        gpstatev(:,k1)  = statev
        
        if (pnewdt<1.d0) then
            call write_output('material routine (stp='//int2str(kstep)//', incr='//int2str(kinc)//', e='//int2str(noel)//', gp='//int2str(k1)//') requested a smaller timestep: pnewdt='//dbl2str(pnewdt,'F0.4'), 'status', 'sim:atp')
            exit
        endif
        
    enddo
    
end subroutine element_nlgeom

subroutine element_lingeom( Ke, Re, ue, gpstrain, gpstress, gpstatev, gpF, Rpos, H0, &
kinc, kstep, noel, pnewdt, props, cmname, dtemp, temp, time, dt, umat)
use umat_mod
implicit none
    !other input variables
    double precision:: Ke(ndof,ndof), Re(ndof), Ke_gp(ndof,ndof), Re_gp(ndof)
    double precision:: ue(ndof), Rpos(nnod), H0, time(2)
    double precision:: gpstrain(6,ngp), gpstress(6,ngp), gpstatev(:,:), gpF(9,ngp)
    !umat variables
    procedure(umat_template),pointer:: umat     ! Addresss to umat subroutine
    character(len=80)           :: cmname
    integer                     :: nstatv, nprops, noel, kstep, kinc
    double precision            :: dt,temp,dtemp,pnewdt,celent,sse,spd,scd,rpl,drpldt
    double precision            :: stress(6), statev(size(gpstatev, dim=1)), ddsdde(6,6), ddsddt(6), drplde(6)
    double precision            :: stran(6), dstran(6), abatime(2), props(:), coords(3), drot(3,3), dfgrd0(3,3)
    double precision            :: dfgrd1(3,3), predef(1), dpred(1)
    !internal variables
    double precision            :: dNpdR(nnod), Np(nnod), phi, uz, r, R0, w, F9(9), B_ss(6,ndof), drdR, ue_r(ndof-2)
    integer                     :: k1
    double precision            :: bbar_eps22, bbar_deps22, bbar_B_ss_row2(ndof) ! Variables used if bbar=.true.
    
    celent = Rpos(nnod)-Rpos(1)
    nstatv = size(gpstatev, dim=1)
    nprops = size(props)
    
    ! Convert to abaqus time definition, as it requires the time in the beginning of the increment.
    abatime = time - dt
    
    sse = 0.d0; spd = 0.d0; scd = 0.d0; rpl = 0.d0
    predef = 0.d0; dpred  = 0.d0
    ddsddt = 0.d0; drplde = 0.d0; drpldt = 0.d0
    
    drot = 0.d0;        drot(1,1) = 1.d0
    drot(2,2) = 1.d0;   drot(3,3) = 1.d0
    
    uz  = ue(1); phi = ue(2); ue_r= ue(3:ndof)    ! Radial displacements
    
    Ke = 0.d0
    Re = 0.d0
    
    if (bbar) then
        R0 = (Rpos(1)+Rpos(size(Rpos)))/2.d0    ! Middle point
        call shapefun(R0, Rpos, dNpdR, Np, ue_r, r, drdR)
        call FBeps_calc(F9, B_ss, stran, dstran, dNpdR, Np, drdR, r, phi, uz, R0, H0, gpstrain(:,1))
        bbar_eps22 = stran(2)
        bbar_deps22 = dstran(2)
        bbar_B_ss_row2 = B_ss(2,:)
    endif
    
    
    do k1=1,ngp
        dfgrd0 = matrix(gpF(:,k1))
        call get_gpinfo(R0, w, k1, Rpos)
        call shapefun(R0, Rpos, dNpdR, Np, ue_r, r, drdR)
        call FBeps_calc(F9, B_ss, stran, dstran, dNpdR, Np, drdR, r, phi, uz, R0, H0, gpstrain(:,k1))
        if (bbar) then
            stran(2)  = bbar_eps22
            dstran(2) = bbar_deps22
            B_ss(2,:) = bbar_B_ss_row2
        endif
        
        dfgrd1 = matrix(F9)
        stress = gpstress(:,k1)
        statev = gpstatev(:,k1)
        coords = (/R0, 0.d0, 0.d0/)
        ! Add updated values to gp-saved variables
        ! The old values of these are saved in the main program 
        ! and are used to reset the updated values for each increment
        gpF(:,k1)       = F9
        gpstrain(:,k1)  = stran+dstran
        ! Need to set it here for F and stran variables, to avoid being affected by changes within the UMAT subroutine
        
        call umat(stress,statev,ddsdde,sse,spd,scd, &
                rpl,ddsddt,drplde,drpldt,&
                stran,dstran,abatime,dt,temp,dtemp,predef,dpred,cmname,&
                3,3,6,nstatv,props,nprops,coords,drot,pnewdt,&
                celent,dfgrd0,dfgrd1,noel,k1,1,1,kstep,kinc)
        call Re_ss_calc(Ke_gp, Re_gp, B_ss, stress, ddsdde) !Calculate Ke_gp and Re_gp
        Ke = Ke + 2.d0*pi*R0*H0*w*Ke_gp
        Re = Re + 2.d0*pi*R0*H0*w*Re_gp
        
        ! Add updated values to gp-saved variables
        ! The old values of these are saved in the main program 
        ! and are used to reset the updated values for each increment
        gpstress(:,k1)  = stress
        gpstatev(:,k1)  = statev
        
        if (pnewdt<1.d0) then
            call write_output('material routine (stp='//int2str(kstep)//', incr='//int2str(kinc)//', e='//int2str(noel)//', gp='//int2str(k1)//') requested a smaller timestep: pnewdt='//dbl2str(pnewdt,'F0.4'), 'status', 'sim:atp')
            exit
        endif
        
    enddo
    
end subroutine element_lingeom

! Calculate position R0 and weight for current gauss point   
subroutine get_gpinfo(R0, w, gp, Rpos)
implicit none
    double precision, intent(out)   :: R0, w
    double precision, intent(in)    :: Rpos(nnod)
    integer, intent(in)             :: gp
    double precision                :: pos(ngp), wght(ngp)
    double precision                :: tmp1, tmp2, tmp3
    ! Weight and positions taken from Ottosen & Petersson "Introduction to the Finite Element Method" Table 20.1
    if (ngp==1) then
        pos = (/0.d0/)
        wght = (/2.d0/)
    elseif (ngp==2) then
        pos = (/-1.d0, 1.d0/)/sqrt(3.d0)
        wght = (/1.d0, 1.d0/)
    elseif (ngp==3) then
        tmp1 = sqrt(3.d0/5.d0)
        pos = (/-tmp1, 0.d0, tmp1/)
        wght = (/5.d0, 8.d0, 5.d0/)/9.d0
    elseif (ngp==4) then
        tmp1 = 0.861136311594053d0
        tmp2 = 0.339981043584856d0
        pos = (/-tmp1, -tmp2, tmp2, tmp1/)
        tmp1 = 0.347854845137454d0
        tmp2 = 0.652145154862546d0
        wght= (/tmp1, tmp2, tmp2, tmp1/)
    elseif (ngp==5) then
        tmp1 = 0.906179845938664d0
        tmp2 = 0.538469310105683d0
        pos = (/-tmp1, -tmp2, 0.d0, tmp2, tmp1/)
        tmp1 = 0.236926885056189d0
        tmp2 = 0.478628670499366d0
        tmp3 = 0.568888888888889d0
        wght = (/tmp1, tmp2, tmp3, tmp2, tmp1/)
    elseif (ngp==6) then
        tmp1 = 0.932469514203152d0
        tmp2 = 0.661209386466265d0
        tmp3 = 0.238619186083197d0
        pos = (/-tmp1, -tmp2, -tmp3, tmp3, tmp2, tmp1/)
        tmp1 = 0.171324492379170d0
        tmp2 = 0.360761573048139d0
        tmp3 = 0.467913934572691d0
        wght = (/tmp1, tmp2, tmp3, tmp3, tmp2, tmp1/)
    else
        call write_output('Min ngp=1, max ngp=6', 'error', 'sim:atp')
    endif
    
    R0 = (Rpos(nnod)-Rpos(1))*pos(gp)/2.d0 + (Rpos(nnod)+Rpos(1))/2.d0  !Change domain
    w = (Rpos(nnod)-Rpos(1))*wght(gp)/2.d0                              !Change domain

end subroutine

! Calculate Ke and Re contribution from current gauss point
subroutine Re_calc(Ke, Re, F9, B_ki, C_kij, stress, ddsdde)
implicit none
    double precision, intent(out)   :: Ke(ndof,ndof), Re(ndof)
    double precision, intent(in)    :: F9(9), B_ki(9,ndof), C_kij(9,ndof,ndof), stress(6), ddsdde(6,6)
    double precision                :: P9(9), dPdF(9,9)
    
    call dPdF_calc(P9, dPdF, stress, ddsdde, F9)
    
    Re = matmul(P9,B_ki)
    Ke = m1_m3_mul(P9, C_kij) + matmul(transpose(B_ki), matmul(dPdF, B_ki))
    ! Preformance can be enhanced by using dgemm!
end subroutine

subroutine Re_calc2(Ke, Re, F9, stress, ddsdde, Np, dNpdR, r, phi, H0, R0)
implicit none
    integer, parameter, dimension(4):: inds = (/1,2,3,5/)
    double precision, intent(out)   :: Ke(ndof,ndof), Re(ndof)
    double precision, intent(in)    :: F9(9), stress(6), ddsdde(6,6)
    double precision, intent(in)    :: Np(ndof-2), dNpdR(ndof-2), r, phi, H0, R0
    double precision                :: P9(9), dPdF(9,9)
    double precision                :: Bred(4,ndof)
    
    call dPdF_calc(P9, dPdF, stress, ddsdde, F9)
    
    
    Bred = 0.d0
    Bred(1, 3:ndof) =  dNpdR
    Bred(2, 3:ndof) =  Np/R0
    Bred(3, 1)      =  1.d0/H0
    Bred(4, 2)      =  r/H0
    Bred(4, 3:ndof) =  Np*phi/H0
    !Re = matmul(P9(inds),Bred)
    
    
    Re(1) = P9(3)/H0
    Re(2) = P9(5)*r/H0
    Re(3:ndof) = P9(1)*dNpdR + P9(2)*Np/R0 + P9(5)*Np*phi/H0
    
    ! Calculate C matrix (dBdu)
    !C_kij(5, 3:ndof, 2) = Np/H0
    !C_kij(5, 2, 3:ndof) = Np/H0
    !Ke = m1_m3_mul(P9, C_kij) + matmul(transpose(B_ki), matmul(dPdF, B_ki))
    Ke = matmul(transpose(Bred), matmul(dPdF(inds, inds), Bred))
    Ke(3:ndof, 2) = Ke(3:ndof, 2) + P9(5)*Np/H0
    Ke(2, 3:ndof) = Ke(2, 3:ndof) + P9(5)*Np/H0
    
end subroutine


subroutine Re_ss_calc(Ke, Re, B_ss, stress, ddsdde)
implicit none
    double precision, intent(out)   :: Ke(ndof,ndof), Re(ndof)
    double precision, intent(in)    :: B_ss(6,ndof),stress(6), ddsdde(6,6)
!   double precision                :: tmp(6,ndof)
    
    Re = matmul(stress, B_ss)
    Ke = matmul(transpose(B_ss), matmul(ddsdde, B_ss))
    !            1    2     3     4  5     6       7     8     9   10    11   12   13
    !           ta,  tb,    m,    n, k, alph,      a,  lda,    b, ldb, beta,   c, ldc
!   call dgemm('N', 'N',    6, ndof, 6, 1.d0, ddsdde,    6, B_ss,   6, 0.d0, tmp, 6)
!   call dgemm('T', 'N', ndof, ndof, 6, 1.d0,   B_ss,    6,  tmp,   6, 0.d0,  Ke, ndof)
    ! Using dgemm decrease performance when run code optimized. 
    ! In debug mode a large speed increase is seen, but this is of coarse not interesting...
end subroutine

! Calculate the shape function value and its derivative at certain position
subroutine shapefun(R0, Rnods, dNpdR, Np, ur, r, drdR)
implicit none
    double precision, intent(in)    :: R0, Rnods(nnod), ur(nnod)
    double precision, intent(out)   :: dNpdR(nnod), Np(nnod), r, drdR
    integer             :: k1
    
    ! Calculate shape functions' value and their derivatives at R0
    if(Nnod==2) then
        Np = (/Rnods(2)-R0, R0-Rnods(1)/)/(Rnods(2)-Rnods(1))
        dNpdR = (/-1.d0, 1.d0/)/(Rnods(2)-Rnods(1))
    elseif(Nnod==3) then
        Np(1) = (R0-Rnods(2))*(R0-Rnods(3))/( (Rnods(1)-Rnods(2))*(Rnods(1)-Rnods(3)) )
        Np(2) = (R0-Rnods(1))*(R0-Rnods(3))/( (Rnods(2)-Rnods(1))*(Rnods(2)-Rnods(3)) )
        Np(3) = (R0-Rnods(1))*(R0-Rnods(2))/( (Rnods(3)-Rnods(1))*(Rnods(3)-Rnods(2)) )
        
        dNpdR(1) = (2*R0 - Rnods(2) - Rnods(3))/( (Rnods(1)-Rnods(2))*(Rnods(1)-Rnods(3)) )
        dNpdR(2) = (2*R0 - Rnods(1) - Rnods(3))/( (Rnods(2)-Rnods(1))*(Rnods(2)-Rnods(3)) )
        dNpdR(3) = (2*R0 - Rnods(1) - Rnods(2))/( (Rnods(3)-Rnods(1))*(Rnods(3)-Rnods(2)) )
    else
        call write_output('Number of nodes chosen not implemented', 'error', 'sim:atp')
    endif
    
    ! Calculate the updated r coordinate along with the deformation gradient component dr/dR
    r = R0
    drdR = 1.d0
    do k1=1,nnod
        r    =   r  +    Np(k1)*ur(k1)
        drdR = drdR + dNpdR(k1)*ur(k1)
    enddo
    
    
end subroutine

!Calculate the matrices F, B and C
subroutine FBC_calc(F9, B_ki, C_kij, dNpdR, Np, drdR, r, phi, uz, R0, H0)
implicit none
    double precision, intent (out)  :: F9(9), B_ki(9,ndof), C_kij(9,ndof,ndof)
    double precision, intent(in)    :: dNpdR(Nnod), Np(Nnod), drdR, r, phi, uz, R0, H0
    F9    = 0.d0
    B_ki  = 0.d0
    C_kij = 0.d0
    
    ! Calculate deformation gradient
    F9(1) = drdR
    F9(2) = r/R0
    F9(3) = 1 + uz/H0
    F9(5) = r*phi/H0
    
    ! Calculate B-matrix dFdu
    B_ki(1, 3:ndof) =  dNpdR
    B_ki(2, 3:ndof) =  Np/R0
    B_ki(3, 1)      =  1.d0/H0
    B_ki(5, 2)      =  r/H0
    B_ki(5, 3:ndof) =  Np*phi/H0
    
    ! Calculate C matrix (dBdu)
    C_kij(5, 3:ndof, 2) = Np/H0
    C_kij(5, 2, 3:ndof) = Np/H0
    
end subroutine

subroutine FBC_calc2(F9, dNpdR, Np, drdR, r, phi, uz, R0, H0)
implicit none
    double precision, intent (out)  :: F9(9)
    double precision, intent(in)    :: dNpdR(Nnod), Np(Nnod), drdR, r, phi, uz, R0, H0
    F9    = 0.d0
    
    ! Calculate deformation gradient
    F9(1) = drdR
    F9(2) = r/R0
    F9(3) = 1 + uz/H0
    F9(5) = r*phi/H0
    
end subroutine

!Calculate the matrices F, B and C
subroutine FBeps_calc(F9, B_ss, eps, deps, dNpdR, Np, drdR, r, phi, uz, R0, H0, eps_old)
implicit none
    double precision, intent(out)   :: F9(9), B_ss(6,ndof), eps(6), deps(6)
    double precision, intent(in)    :: dNpdR(Nnod), Np(Nnod), r, phi, uz, R0, H0, drdR, eps_old(6)
    
    ! Calculate deformation gradient (actual, not assuming small strains)
    F9    = 0.d0
    F9(1) = drdR
    F9(2) = r/R0
    F9(3) = 1 + uz/H0
    F9(5) = r*phi/H0

    ! Calculate small strain tensors
    ! Note abaqus order: [11, 22, 33, 12, 13, 23]
    ! Abaqus wants old strain and strain increment!
    eps   = 0.d0
    eps(1) = drdR-1.d0
    eps(2) = (r/R0)-1.d0
    eps(3) = uz/H0
    eps(6) = 2.d0*R0*phi/(2*H0) ! Engineering strain
    
    deps = eps - eps_old
    
    eps = eps_old
    
    
    ! Calculate b-matrix deps_du
    B_ss  = 0.d0
    B_ss(1, 3:ndof) =  dNpdR
    B_ss(2, 3:ndof) =  Np/R0
    B_ss(3,      1) =  1/H0
    B_ss(6,      2) =  R0/H0
    
end subroutine

end module
