module mps_util_mod
use sim_tensor_mod
use sim_util_mod
use umat_mod
use output_mod
implicit none
    
    private
    public      ::  mps_gen_free_dofs
    public      ::  mps_solve_incr
    
    contains 
    
subroutine mps_gen_free_dofs(free_dofs, ctrl)
implicit none
integer, allocatable        :: free_dofs(:)
integer, intent(in)         :: ctrl(:)

integer                     :: nfree
integer                     :: k1, k2

nfree = count(ctrl<0)
if (nfree==0) then  ! Special indicator if all dofs are prescribed
    allocate(free_dofs(1))
    free_dofs = -1
else
    allocate(free_dofs(nfree))
    k2 = 0
    do k1=1,size(ctrl)
        if (ctrl(k1)<0) then
            k2 = k2 + 1
            free_dofs(k2) = k1
        endif
    enddo
endif

end subroutine


subroutine mps_solve_incr(load_old, disp_old, temp_old, stat_old, load, disp, temp, stat, time, dt, free_dofs, iter, niter, lconv, pnewdt, &
    kinc, kstep, props, cmname, umat, nlgeom)
    use types_mod
    use umat_mod
    implicit none
        double precision, parameter, dimension(9)   :: I9 = [1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0]
        !Input/Output variables
        double precision, intent(in)    :: load_old(:), disp_old(:), temp_old, stat_old(:)
        double precision, intent(inout) :: load(:), disp(:), temp, stat(:)
        double precision, intent(in)    :: time(2), dt
        integer, intent(in)             :: free_dofs(:)
        type(iter_typ), intent(in)      :: iter
        integer, intent(out)            :: niter
        logical, intent(out)            :: lconv
        double precision, intent(out)   :: pnewdt
        integer, intent(in)             :: kinc, kstep  ! Increment and step number
        double precision, intent(in)    :: props(:)     ! Material properties/parameters
        character(len=80), intent(in)   :: cmname       ! Material name sent to umat
        procedure(umat_template),pointer:: umat         ! Addresss to umat subroutine
        logical, intent(in)             :: nlgeom       ! Bool to determine if nonlinear geometry effects should be accounted for
        
        !Internal variables
        double precision :: K(size(disp), size(disp)), R(size(disp)), load_res(size(load))
        double precision :: du(size(disp)), F9(9), F9_old(9)
        double precision :: error, error_old, alpha
        !Additional UMAT variables
        double precision :: stress(6), stran(6), dstran(6), ddsdde(6,6), dfgrd0(3,3), dfgrd1(3,3), drot(3,3)
        double precision :: sse, spd, scd, rpl, ddsddt(6), drplde(6), drpldt, predef(1), dpred(1)
        double precision :: abatime(2), coords(3)
        integer          :: nprops, nstatv
        
        ! Convert to abaqus time definition, as it requires the time in the beginning of the increment.
        abatime = time - dt
        coords = 0.d0
    
        ! Initiate variables various variables needed for umat, but not necessarily used
        sse = 0.d0; spd = 0.d0; scd = 0.d0; rpl = 0.d0
        predef = 0.d0; dpred  = 0.d0
        ddsddt = 0.d0; drplde = 0.d0; drpldt = 0.d0
    
        drot = 0.d0; ! Set drot to identity matrix
        drot(1,1) = 1.d0; drot(2,2) = 1.d0; drot(3,3) = 1.d0
        
        nprops = size(props)
        nstatv = size(stat)
        
        ! Initiation of variables
        error_old = huge(1.d0)  ! Ensure that error is decreasing the first time, to avoid line search. 
        alpha = 1.d0            ! Initial alpha line search factor
        pnewdt = 1.d0
        lconv = .true.
        niter = 0
        du = disp-disp_old
        do while(.true.)
            R       = 0.d0
            K       = 0.d0
            stat    = stat_old  ! State variables should be reset every time

            if (nlgeom) then
                ! Calculate old stress/dfgrd0/stran and dstran
                F9_old = disp_old+I9
                stress = P9_2_sigma6(load_old, F9_old)
                dfgrd0 = matrix(F9_old)
                F9     = disp+I9
                dfgrd1 = matrix(F9)
                call GL_stran_calc(F9_old, F9, stran, dstran)
                
                ! Call umat
                call umat(stress,stat,ddsdde,sse,spd,scd, &
                rpl,ddsddt,drplde,drpldt,&
                stran,dstran,abatime,dt,temp_old,temp-temp_old,predef,dpred,cmname,&
                3,3,6,nstatv,props,nprops,coords,drot,pnewdt,&
                0.d0,dfgrd0,dfgrd1,1,1,1,1,kstep,kinc)
                
                ! Convert stress=Cauchy to load_res=PK1 and stiffness from ddsdde to K=dPdF
                call dPdF_calc(load_res, K, stress, ddsdde, F9)
                
            else
                stran   = disp_old
                dstran  = disp-disp_old
                stress  = load_old
                
                call umat(stress,stat,ddsdde,sse,spd,scd,rpl,ddsddt,drplde,drpldt,&
                stran,dstran,abatime,dt,temp_old,temp-temp_old,&
                predef,dpred,cmname,3,3,6,nstatv,props,nprops,&
                coords,drot,pnewdt,0.d0,dfgrd0,dfgrd1,1,1,1,1,kstep,kinc)
                
                K = ddsdde
                load_res = stress    
            endif
            R = load_res - load
            
            lconv = pnewdt>=1.d0
            if (.not.lconv) then
                exit
            endif

            ! Step 2: Check residual
            if (free_dofs(1)/=-1) then
                error = sqrt(sum(R(free_dofs)**2))
            else
                error = 0.d0
            endif
            
            ! Step 2: Solve increment update, depending on error
            if (error<iter%tol) then
                lconv = .true.
                load = load_res
                exit
            elseif ((error>error_old).and.(niter>iter%ls_start)) then
                disp = disp - du
                alpha = alpha*iter%ls_alpha_factor
                du = alpha*du
            else !tol<error<error_old
                alpha = 1.d0
                call solve_mateq(K, R, du, lconv, free_dofs, free_dofs)
                
                if (.not.lconv) then
                    pnewdt = 0.5d0
                    call write_output('Could not solve matrix equation, attempting a smaller time increment', 'status', 'sim:mps')
                    lconv = .false.
                    exit
                endif
            endif
            
            ! Step 3: Update solution
            disp = disp + du
            niter = niter + 1
            error_old = error
            
            ! Step 4: Check for max number of iterations
            if (niter>iter%max) then
                lconv = .false.
                pnewdt = 0.75d0
                call write_output('Could not find equilibrum in '//int2str(iter%max)//' attempts, using smaller time step', 'status', 'sim:mps')
                exit
            endif
        enddo
end subroutine

end module
    