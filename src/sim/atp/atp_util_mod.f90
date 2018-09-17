    module atp_util_mod
    use sim_tensor_mod
    use sim_util_mod
    use atp_element_mod
    use output_mod
    implicit none

    private
    
    public  ::  solve_incr      ! Attempt to solve one increment
    public  ::  apply_load      ! Apply the external loads
    public  ::  apply_bc        ! Apply the boundary conditions (prescribed displacements)
    public  ::  get_result      ! Get the results from the particular increment
    public  ::  gen_free_dofs   ! Obtain the free degrees of freedom for particular step
    
    ! Parameters
    double precision, parameter :: pi = acos(-1.d0)  !Define pi

    contains
    !==========================================================================================================
    !Iteration subroutines

!Solve one increment
subroutine solve_incr(Rpos, H0, load, disp, temp, dtemp, time, dt, free_dofs, known_dofs, gp_F, gp_s, gp_sv, gp_strain, u, du, &
    iter, niter, lconv, pnewdt, kinc, kstep, props, cmname, umat_address, nlgeom, iter_err_norm)
    use types_mod
    use umat_mod
    implicit none
        !Input/Output variables
        double precision, intent(in)    :: Rpos(:,:), H0   !Element nodal coordinates, dim = (nnod, nel) (nnod = number of nodes in element)
        double precision, intent(inout) :: load(:), disp(:) !External loads and displacements
        double precision, intent(in)    :: temp, dtemp, time(2), dt ! Time and temperature
        double precision, intent(inout) :: gp_F(:,:,:), gp_s(:,:,:), gp_sv(:,:,:), gp_strain(:,:,:)  !Values at gauss points (dim(*),ngp,nel)
        double precision, intent(inout) :: u(:), du(:)  !Displacements (uz, phi, ur1,ur2,ur3,...,urN)
        type(iter_typ)                  :: iter         ! Iteration settings
        double precision, intent(out)   :: pnewdt       ! Variable used to request shorter time step
        double precision, intent(in)    :: props(:)     ! Material parameters
        integer, intent(out)            :: niter        ! Number of iterations used
        integer, intent(in)             :: kinc, kstep  ! Current increment and step
        integer, intent(in)             :: free_dofs(:) ! List over the non-prescribed (i.e. free) degrees of freedom
        integer, intent(in)             :: known_dofs(:)! List over the degrees of freedom with known load
        logical, intent(out)            :: lconv        ! Local convergence bool
        character(len=80), intent(in)   :: cmname       ! Material name sent to umat
        procedure(umat_template),pointer:: umat_address ! Addresss to umat subroutine
        logical, intent(in)             :: nlgeom       ! Bool to determine if nonlinear geometry effects should be accounted for
        double precision, intent(in)    :: iter_err_norm(:)	!Scaling to obtain averaged stresses from the residuals
        
        !Internal variables
        double precision :: K(size(u), size(u)), R(size(u))
        double precision :: error, error_old, alpha
        double precision :: Ke(size(Rpos,1)+2, size(Rpos,1)+2), Re(size(Rpos,1)+2)
        integer          :: nel, nnod, ndof_tot, k1, k2, dof(size(Rpos,1)+2)
        double precision :: gp_F0(9, size(gp_F, dim=2), size(gp_F, dim=3))
        double precision :: gp_s0(6, size(gp_s, dim=2), size(gp_s, dim=3))
        double precision :: gp_strain0(6, size(gp_strain,dim=2), size(gp_strain,dim=3))
        double precision :: gp_sv0(size(gp_sv,dim=1), size(gp_sv,dim=2), size(gp_sv,dim=3))
        double precision :: u_old(size(u)), ue(size(Rpos,1)+2)
        double precision :: Rscale      ! Scale factor for the equilibrium equations
        
        nnod= size(Rpos,1) !Number of nodes per element
        nel = size(Rpos,2) !Number of elements
        ndof_tot = 2 + 1+(nnod-1)*nel

        lconv = .true.

        dof(1:2) = (/1,2/) ! All elements have the axial and torsion dofs

        error_old = huge(1.d0)  ! Ensure that error is decreasing the first time, to avoid line search. 
        alpha = 1.d0            ! Initial alpha line search factor
    
        u_old = u   !Save old u to make it possible to calculate total du as input to the next step's guess
        u = u + du  !du contains correct displacements for prescribed displacement increment

        gp_F0 = gp_F
        gp_s0 = gp_s
        gp_sv0= gp_sv
        gp_strain0 = gp_strain
        pnewdt = 1.d0
        niter = 0
        do while(.true.)
            ! Step 0: Reset values for each iteration attempt
            gp_F  = gp_F0
            gp_s  = gp_s0
            gp_sv = gp_sv0
            gp_strain = gp_strain0 
            
            R       = 0.d0
            Rscale  = 0.d0
            K       = 0.d0

            ! Step 1: Assemble system
            do k1=1,nel
                do k2=1,nnod
                    dof(k2+2) = (k2+2) + (nnod-1)*(k1-1)
                enddo

                ! Get element displacements
                ue = u(dof)
                
                ! Call element routine
                if (nlgeom) then
                    call element_nlgeom( Ke, Re, ue, gp_s(:,:,k1), gp_sv(:,:,k1), gp_F(:,:,k1), Rpos(:,k1), H0, &
                    kinc, kstep, k1, pnewdt, props, cmname, dtemp, temp, time, dt, umat_address)
                else
                    call element_lingeom( Ke, Re, ue, gp_strain(:,:,k1), gp_s(:,:,k1), gp_sv(:,:,k1), gp_F(:,:,k1), Rpos(:,k1), H0, &
                    kinc, kstep, k1, pnewdt, props, cmname, dtemp, temp, time, dt, umat_address)
                endif
                

                if (pnewdt<1.d0) then
                    lconv = .false.
                    exit
                endif

                K(dof,dof)  = K(dof, dof) + Ke
                R(dof)      = R(dof)      + Re
                Rscale      = Rscale + sqrt(sum((Re*iter_err_norm(dof))**2))	! Add the average stress contributions to Rscale
            enddo
            if (.not.lconv) then
                exit
            endif

            ! Apply external loads (include potential effect on stiffness if nlgeom=true)
            call apply_load(R, K, load, H0, Rpos, u, known_dofs, nlgeom)
            
            ! Step 2: Check residual
            if (free_dofs(1)/=-1) then
                error = sqrt(sum((R(known_dofs)*iter_err_norm(known_dofs))**2))/max(Rscale,1.d0)
            !   error = sqrt(sum(R(known_dofs)**2))
            else
                error = 0.d0
            endif
            
            ! Step 2: Solve increment update, depending on error
            if (error<iter%tol) then
                lconv = .true.
                exit
            elseif ((error>error_old).and.(niter>iter%ls_start)) then
                u = u - du
                alpha = alpha*iter%ls_alpha_factor
                du = alpha*du
            else !tol<error<error_old
                alpha = 1.d0
                call solve_mateq(K, R, du, lconv, free_dofs, known_dofs)
                
                if (.not.lconv) then
                    pnewdt = 0.5d0
                    call write_output('Could not solve matrix equation, attempting a smaller time increment', 'status', 'sim:atp')
                    lconv = .false.
                    exit
                endif
            endif
            
            ! Step 3: Update solution
            u = u + du
            niter = niter + 1
            error_old = error
            
            ! Step 4: Check for max number of iterations
            if (niter>iter%max) then
                lconv = .false.
                pnewdt = 0.75d0
                call write_output('Could not find equilibrum in '//int2str(iter%max)//' attempts, using smaller time step', 'status', 'sim:atp')
                exit
            endif
        enddo
        
        du = u - u_old  !Output total displacement increment
        
        ! Insert the previously unknown values into load/disp vectors
        call get_result(R, load, disp, H0, Rpos, u, free_dofs, known_dofs, nlgeom)

end subroutine


!Load and boundary conditions
subroutine apply_load(R, K, load, H0, Rpos, u, known_dofs, nlgeom)
    implicit none
    double precision    :: R(:), K(:,:), load(:), H0, Rpos(:,:), u(:)
    integer             :: known_dofs(:), n
    logical             :: nlgeom
    double precision    :: h, ri, ro
        
    n    = size(u)
        
    h = H0                                  ! Height
    ri = Rpos(1,1)                          ! Inner radius
    ro = Rpos(size(Rpos,1),size(Rpos,2))    ! Outer radius
            
    if (nlgeom) then    ! get current:
        h  = H0 + u(1)  ! height,
        ri = ri + u(3)  ! inner radius and
        ro = ro + u(n)  ! outer radius
    endif
        
    ! Axial force
    if (any(known_dofs==1)) then
        R(1) = R(1) - load(1)
    endif
        
    ! Torque
    if (any(known_dofs==2)) then
        R(2) = R(2) - load(2)
    endif
        
    ! Internal pressure
    if (any(known_dofs==3)) then
        R(3) = R(3) - load(3)*(2.d0*ri*pi)*h
        if (nlgeom) then
            K(3,1)  = K(3,1) - load(3)*(2.d0*ri*pi)
            K(3,3)  = K(3,3) - load(3)*(2.d0*pi)*h
        endif   
    endif
        
    ! External pressure
    if (any(known_dofs==n)) then
        R(n) = R(n) + load(4)*(2.d0*ro*pi)*h    
        if (nlgeom) then
            K(n,1)  = K(n,1) + load(4)*(2.d0*ro*pi)
            K(n,n)  = K(n,n) + load(4)*(2.d0*pi)*h
        endif
    endif

end subroutine
    
subroutine apply_bc(disp, disp_new, du, ctrl, disp_conv)
    implicit none
    double precision, intent(in)    :: disp(:), disp_new(:), disp_conv(:)
    integer, intent(in)             :: ctrl(:)
    double precision                :: du(:)
    integer                         :: ctrl_dofs(size(ctrl)), k1

    ctrl_dofs = (/1, 2, 3, size(du)/)
    
    do k1=1,4
        if (ctrl(k1)>0) then !Displacement control
            du(ctrl_dofs(k1)) = (disp_new(k1)-disp(k1))*disp_conv(k1)
        endif
    enddo

end subroutine

subroutine get_result(R, load, disp, H0, Rpos, u, free_dofs, known_dofs, nlgeom)
    implicit none
    double precision    :: R(:), load(:), disp(:), H0, Rpos(:,:), u(:)
    integer             :: free_dofs(:), known_dofs(:)
    logical             :: nlgeom
    integer             :: n, pos(4), k1
    double precision    :: h, ri, ro, disp_conv(4)
    
    n    = size(u)
    pos = (/1,2,3,n/)
    h = H0                                  ! Height
    ri = Rpos(1,1)                          ! Inner radius
    ro = Rpos(size(Rpos,1),size(Rpos,2))    ! Outer radius
    
    ! Note that this is not the same disp_conv as the globally defined (if ri<1.d-12)
    disp_conv = (/H0, 1.d0, ri, ro/)        ! Must be before nlgeom changes below.
    if (ri<1.d-12) then  ! If inner radius = 0, then strain is infinity. Hence the output is the inner displacement which should be zero!
        disp_conv(3) = huge(1.d0)
    endif
    
            
    if (nlgeom) then    ! get current:
        h  = H0 + u(1)  ! height,
        ri = ri + u(3)  ! inner radius and
        ro = ro + u(n)  ! outer radius
    endif
    
    ! Update outputs
    do k1 = 1,4
        if (any(free_dofs==pos(k1))) then       ! Displacement to be determined
            disp(k1) = u(pos(k1))/disp_conv(k1) ! Obtain resulting strains
        endif
        if (.not.any(known_dofs==pos(k1))) then ! Force to be determined
            if (pos(k1)<=2) then                ! Axial force or torque
                load(k1) = R(pos(k1))           ! Obtain resulting load for axial force or torque
            elseif (pos(k1)==3) then            ! Inner pressure
                if (ri<1.d-12) then             ! If inner radius is equal to zero (i.e. very small)
                    load(3) = 0.d0              ! Set load to zero
                else
                    load(3) = R(pos(k1))/(2*pi*ri*h)
                endif
            else                                ! Outer pressure
                load(4) = -R(pos(k1))/(2*pi*ro*h)
            endif
        endif
    enddo
    
end subroutine

!Generate free_dofs
subroutine gen_free_dofs(free_dofs, known_dofs, ndof_tot, ctrl, intpres_extstrn)
    implicit none
        integer, allocatable :: free_dofs(:)            ! The free degrees of freedom
        integer, allocatable :: known_dofs(:)           ! The degrees of freedom with known load
        integer              :: ctrl(4)                 ! The control parameters (axial, torsion, inside, outside)
        logical              :: intpres_extstrn         ! If internal pressure and external strain used with strain control
        integer              :: ctrl_dofs(4)            ! Degree of freedom for dofs directly controlled by ctrl 
        integer              :: ndof_tot, nfree
        integer              :: k1, k2, k3
    
        nfree = ndof_tot-count(ctrl>0)  ! Get number of free dofs
        
        if (nfree==0) then  ! If no free degrees of freedom, we handle this special case separately!
            allocate(free_dofs(1))
            free_dofs = -1
        else
            allocate(free_dofs(nfree))
            ctrl_dofs = (/1, 2, 3, ndof_tot/)
            ! Find free dofs amongst the first specified in ctrl
            k2 = 1
            do k1 = 1,3
                if (ctrl(k1)<0) then ! Force controlled
                    free_dofs(k2) = ctrl_dofs(k1)
                    k2 = k2 + 1
                endif
            enddo
        
            if (k2<=nfree) then  ! There are additional dofs 
                k3 = 4  ! Three dofs (axial, tors, inner are before the additional dofs)
                do k1=k2,nfree
                    free_dofs(k1) = k3
                    k3 = k3 + 1
                enddo
            endif
        
            ! We want sorted dofs, hence if ctrl(4)<0 then this dof should be put at the end
            ! Then the last is overwritten, but it is actually ndof_tot so the code below is not needed:
            !if (ctrl(4)<0) then
            !    free_dofs(nfree) = ndof_tot
            !endif
            ! Normally, the known_dofs (load known) corresponds to free_dofs (disp unknown)
            allocate(known_dofs, source=free_dofs)
            ! In the special case of known external pressure and external strain, 
            ! with unknown internal disp and internal pressure, they are different:
            if (intpres_extstrn) then
                k1 = count(ctrl(1:2)<0)+1
                if (free_dofs(nfree)==ndof_tot) then
                    ! For special intpres_extstrn mode the external strain should be controlled
                    call write_output('For the special intpres_extstrn mode the external strain must be controlled', 'error', 'sim:atp')
                elseif (free_dofs(k1)/=3) then
                    call write_output('For the special intpres_extstrn mode the internal pressure should be set as controlled', 'error', 'sim:atp')
                endif
                
                known_dofs(k1:) = free_dofs(k1:)+1
            endif
            
                
        endif

end subroutine

end module
