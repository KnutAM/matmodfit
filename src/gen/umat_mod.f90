module umat_mod
    implicit none
    
    private
    public  ::  umat_template
    
    ! Interface for material routine
    abstract interface
        subroutine umat_template(stress,statev,ddsdde,sse,spd,scd, &
        rpl,ddsddt,drplde,drpldt,&
        stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,&
        ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,&
        celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
        
        character(len=80) cmname
        integer:: ndi,nshr,ntens,nstatv,nprops,noel,npt,&
        layer, kspt, kstep, kinc
        double precision:: dtime,temp,dtemp,pnewdt,celent,sse,spd,scd,rpl,drpldt
        double precision:: stress(ntens),statev(nstatv),&
        ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),&
        stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),&
        props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
        end subroutine
    end interface
    
    contains 
        
end module umat_mod
    
