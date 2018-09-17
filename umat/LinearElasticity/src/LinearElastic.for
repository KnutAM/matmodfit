!DEC$ FREEFORM
! Material parameters props = (E, nu)
! E:    Youngs modulus of elasticity
! nu:   Poissons ratio
!
subroutine umat(stress,statev,ddsdde,sse,spd,scd, &
    rpl,ddsddt,drplde,drpldt,&
    stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
    ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
    celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
!DEC$ ATTRIBUTES DLLEXPORT :: UMAT
    implicit none
        
    character(len=80) cmname
    integer:: ndi,nshr,ntens,nstatv,nprops,noel,npt,&
    layer, kspt, kstep, kinc
    double precision:: dtime,temp,dtemp,pnewdt,celent,sse,spd,scd,rpl,drpldt
    double precision:: stress(ntens),statev(nstatv),&
    ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),&
    stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),&
    props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
    
    double precision :: E, nu, mu, la   ! Material parameters
    integer          :: k1              ! Looping iterator
	
    E =  props(1)
	nu = props(2)
    
    
    la = (E*nu)/( (1.d0 + nu)*(1.d0-2*nu) )
    mu = E/(1.d0 + nu)
    
    ! Calculate stiffness matrix
    ddsdde(1:3, 1:3) = la
    do k1=1,3
        ddsdde(k1,k1) = ddsdde(k1,k1) + 2*mu
    enddo
    do k1=4,6
        ddsdde(k1,k1) = ddsdde(k1,k1) + mu
    enddo
    
    ! Calculate stress
    stress = matmul(stran+dstran, ddsdde)
		
end subroutine