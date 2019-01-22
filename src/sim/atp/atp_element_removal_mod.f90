
module atp_element_removal_mod
use atp_mod
use types_mod
use atp_element_mod
use output_mod
implicit none
    private
    public  :: atp_element_removal ! Axial-torsion-pressure element removal routine
    
    contains
! Element removal procedure
! ============================================================================================================
! 	1. Import data from continued analysis
!	2. Solve load = 0.d0 (Should already have been solved, but good to do. Can also be used as relaxation step)
!	3. Remesh with new mesh nodes from subsequent analysis (can choose arbitrary analysis if several, just to get the input data)
!		a. Automatic mesh outside the body to be (coordinates specified by the deformed geometry)
!		b. Interpolate variables (dofs and state variables)
!		c. Output option decide if old values and interpolated values should be written to file (nice to be able to document e.g. how the state variables vary with radius to say something about the difficulty of interpolation)
!	4. Solve load=0
!	5. Remove elements to only have the ones requested by the user remaining
!	6. Solve load = 0 (5 and 6 could be run in several smaller steps by removing 1 and 1 element…)
!	7. Check tolerance on new nodal positions, and if not ok change the guess for node coordinates and go to 3

    
subroutine atp_element_removal(error, props, f_data, simnr)
use atp_util_mod
use atp_element_mod
use atp_import_mod
use sim_error_mod
use sim_getincr_mod
use sim_writing_util_mod
use umat_mod
use output_mod

implicit none
    ! Input parameters
    double precision, intent(out)   :: error            ! Error
    double precision, intent(in)    :: props(:)         ! Material parameters (true, full) input
    type(fdata_typ), intent(inout)     :: f_data           ! Settings input
    integer, intent(in)             :: simnr            ! Simulation number (=> f_data%sim(simnr) is pertinent to current simulation)

    ! Previous analysis relaxation
    integer                         :: prev_simnr           ! Previous simulation number
    type(fdata_typ)                    :: f_data_initial_relax ! Data belonging to relaxation simulation

    ! Node position iterations
    double precision, allocatable   :: node_pos_guess(:), node_pos_result(:)            ! Node positions 
    double precision, allocatable   :: node_pos_guess_old(:), node_pos_result_old(:)    ! Node positions for better guess
    double precision, allocatable   :: node_pos_update(:)                               ! Node position update
    integer, allocatable            :: node_disp_ind(:)                                 ! Indices in dofs belonging to node_pos variables
    integer                         :: k1                                               ! Iteration index
    logical                         :: geom_conv                                        ! Convergence indicator for geometry iterations
    logical                         :: relx_conv                                        ! Convergence indicator for relaxation
    double precision                :: pos_error                                        ! Position error

    ! Problem size convenience variables
    integer                         :: nel, ngp, nnod, ndof
    logical                         :: new_is_solid                                     ! Bool to decide if the new is a solid specimen and 
                                                                                        ! inner radius should be kept exactly at 0.d0
                                                                                        
    ! Variables for solve_incr
    type(fdata_typ)                    :: f_data_relax
    double precision, allocatable   :: rpos(:,:)
    double precision                :: h0, load(4), disp(4), time(2)
    integer                         :: free_dofs(1)
    double precision, allocatable   :: gp_F(:,:,:), gp_stress(:,:,:), gp_sv(:,:,:), gp_strain(:,:,:), u0(:), du(:)
    double precision                :: iter_tol
    integer                         :: niter, iter_max
    logical                         :: lconv
    double precision                :: pnewdt
    double precision, allocatable   :: iter_err_norm(:)

    ! Simulation timing
    double precision                :: starttime, stoptime
    integer                         :: clock_count, clock_rate

    ! Simulation timing
    call system_clock ( clock_count, clock_rate)
    starttime = real(clock_count)/real(clock_rate)

    ! Convenience size variables
    nel = size(f_data%sim(simnr)%mesh1d%node_pos)-1
    ngp = f_data%sim(simnr)%mesh1d%ngp
    nnod= f_data%sim(simnr)%mesh1d%element_order + 1
    ndof= 2 + nel*(nnod-1)+1

    ! == RELAX PREVIOUS SIMULATION == 
    prev_simnr = abs(f_data%sim(simnr)%init%cont_analysis)
    call atp_relax(relx_conv, props, f_data, simnr, prev_simnr, f_data_initial_relax, '_relax0')
    if (.not.relx_conv) then
        error = huge(1.d0)
        return
    endif
    

    ! == SETUP NODE POSITION ITERATIONS ==
    allocate(node_disp_ind(nel+1))
    ! Index in dof for specified node displacement values
    node_disp_ind = (/(k1, k1=3,ndof,(nnod-1))/)
    allocate(node_pos_guess(nel+1), node_pos_result(nel+1))
    allocate(node_pos_guess_old(nel+1), node_pos_result_old(nel+1))
    allocate(node_pos_update(nel+1))
    
    node_pos_guess = f_data%sim(simnr)%mesh1d%node_pos
    node_pos_result= 0.d0

    ! == SETUP NECESSARY VARIABLES FOR remesh AND solve_incr == 
    allocate(rpos(nnod, nel))       ! Undeformed nodal positions
    ! Gauss point values:
    allocate(gp_stress(6,ngp,nel), gp_sv(f_data%glob%nstatv,ngp,nel), gp_strain(6,ngp,nel), gp_F(9,ngp,nel))
    allocate(u0(ndof), du(ndof), iter_err_norm(ndof))   !Nodal values (and normalization for error (not used))

    free_dofs   = -1    ! All displacements prescribed:
    du          = 0.d0  ! -> prescribed to zero
    load        = 0.d0  ! 
    disp        = 0.d0  ! 
    
    ! These variables are not used, but needed for use with solve_incr
    time        = f_data%sim(1)%atp_er%time_remesh
    iter_tol    = 1.d0  ! No effect, must be >0 as error is set to zero at material convergence
    niter       = 0     ! Dummy output variable (should always become 1)
    iter_max    = 0     ! No effect since solve_incr converges directly if material converges locally (should be >=1)
    pnewdt      = 1.d0  ! This indicator has no effect here, but should be 1.d0 according to abaqus umat standards
    iter_err_norm = 1.d0! Error normalization, but not used here because all displacement prescribed

    new_is_solid = (f_data%sim(simnr)%mesh1d%node_pos(1)<1e-12)
    if (new_is_solid) then ! Force inner node to be zero
        node_pos_guess(1) = 0.d0
    endif
    
    geom_conv = .false.
    
    GEOM_ITER_LOOP: do k1=1,f_data%sim(simnr)%atp_er%geom_iter_max
        call remesh(node_pos_guess, f_data, f_data_initial_relax, simnr, rpos, h0, gp_stress, gp_strain, gp_F, gp_sv, u0, new_is_solid)
        
        if (new_is_solid) then
            u0(3) = 0.d0
        endif
        
        ! Solve zero displacement increment
        call solve_incr(rpos, h0, load, disp, f_data%sim(1)%init%temp_init, 0.d0, time, f_data%sim(1)%atp_er%time_remesh, free_dofs, &
               free_dofs, gp_F, gp_stress, gp_sv, gp_strain, u0, du, f_data%sim(1)%iter, niter, lconv, pnewdt, &
               k1, 1, props, f_data%glob%cmname, f_data%glob%umat_address, f_data%glob%nlgeom, iter_err_norm)
    
        if (.not.lconv) then
           call write_output('No convergence after element removal', 'status', 'sim:atp')
           error = huge(1.d0)
           exit GEOM_ITER_LOOP
        endif
   
        ! Set end values in f_data%sim(simnr)
        f_data%sim(simnr)%mesh1d%h0 = h0
        disp = (/u0(1)/h0, u0(2), u0(3)/rpos(1,1), u0(size(u0))/rpos(size(rpos,1), size(rpos,2))/)
        if (new_is_solid) then
            disp(3) = 0.d0
        endif
        call atp_export_end(f_data, simnr, u0, gp_stress, gp_strain, gp_F, gp_sv, disp, load, h0)
   
        ! Relax 
        call atp_relax(relx_conv, props, f_data, simnr, simnr, f_data_relax, '_relax'//int2str(k1), rpos)
        if (.not.relx_conv) then
            error = huge(1.d0)
            return
        endif
   
        node_pos_result = (/rpos(1,:), rpos(size(rpos,1), size(rpos,2))/) + f_data_relax%sim(1)%end_res%u_end(node_disp_ind)
        pos_error = sqrt(  sum( ( (node_pos_result-f_data%sim(simnr)%mesh1d%node_pos) &
                                  /max(1.d-6,f_data%sim(simnr)%mesh1d%node_pos) )**2 )  )
        if (pos_error<f_data%sim(simnr)%atp_er%node_pos_tol) then
            geom_conv = .true.
            error = 0.d0
            exit GEOM_ITER_LOOP
        elseif (k1==1) then
            ! First iteration - fix point update
            node_pos_guess_old = node_pos_guess
            node_pos_result_old = node_pos_result
            node_pos_guess = node_pos_guess + (f_data%sim(simnr)%mesh1d%node_pos-node_pos_result)
            if (new_is_solid) then ! Force inner node to be zero
                node_pos_guess(1) = 0.d0
            endif
        else
            ! Further iteration, secant update?
            node_pos_update = (f_data%sim(simnr)%mesh1d%node_pos-node_pos_result) &
                                *(node_pos_guess-node_pos_guess_old) &
                                /(node_pos_result-node_pos_result_old)
            node_pos_guess_old = node_pos_guess
            node_pos_result_old = node_pos_result
            node_pos_guess = node_pos_guess + node_pos_update
            if (new_is_solid) then ! Force inner node to be zero
                node_pos_guess(1) = 0.d0
            endif
        endif
   
    enddo GEOM_ITER_LOOP
    
    if (.not.geom_conv) then
        call write_output('No convergence for geometry update due to element removal, error=huge', 'status', 'sim:atp')
        error = huge(1.d0)
    endif
    
    call set_end_values(f_data_relax%sim(1), f_data%sim(simnr))

    call system_clock ( clock_count, clock_rate)
    stoptime = real(clock_count)/real(clock_rate)

end subroutine atp_element_removal

subroutine atp_relax(converged, props, f_data, simnr_new, simnr_old, f_data_relax, append_result_name, rpos)
! Solve relaxing before element removal by setting up another simulation
implicit none
    logical, intent(out)            :: converged
    double precision, intent(in)    :: props(:)
    type(fdata_typ), intent(in)     :: f_data
    integer, intent(in)             :: simnr_new, simnr_old
    type(fdata_typ), intent(out)    :: f_data_relax
    character(len=*), intent(in)    :: append_result_name
    double precision, optional      :: rpos(:,:)
    double precision                :: error
    
    ! Setup the f_data for the relax simulation
    ! 1) Read in the data from the old analysis (prior to element removal)
    allocate(f_data_relax%sim(1))
    !If this allocates the sub-arrays depends on compiler options, hence we should always check if they are allocated...
    f_data_relax%glob = f_data%glob 
    f_data_relax%sim(1) = f_data%sim(simnr_old)
    
    ! 2) Change values required for machining
    ! Result output name (append a _relx to separate from old simulation results)
    write(f_data_relax%glob%outname, "(A)") trim(f_data_relax%glob%outname)//trim(append_result_name)
     
    ! Experiment values - values normally setup during the setup_simulation procedure
    if (.not.allocated(f_data_relax%sim(1)%exp%exp_info)) then
        allocate(f_data_relax%sim(1)%exp%exp_info(11))
    elseif (.not.(size(f_data_relax%sim(1)%exp%exp_info)==11)) then
        deallocate(f_data_relax%sim(1)%exp%exp_info)
        allocate(f_data_relax%sim(1)%exp%exp_info(11))
    endif
    f_data_relax%sim(1)%exp%exp_info = (/0, 1, 2, 0, 3, 0, 4, 0, 5, 0, 0/)
    
    if (allocated(f_data_relax%sim(1)%sim_setup%expdata_array)) then
        deallocate(f_data_relax%sim(1)%sim_setup%expdata_array)
    endif
    allocate(f_data_relax%sim(1)%sim_setup%expdata_array(2,5))
    f_data_relax%sim(1)%sim_setup%expdata_array(1,2:) = 1.d0  ! Initial load (will not be used as ctrl is set to -3)
    f_data_relax%sim(1)%sim_setup%expdata_array(2,2:) = 0.d0  ! Final load
    f_data_relax%sim(1)%sim_setup%expdata_array(:,1) = (/0.d0, f_data%sim(simnr_new)%atp_er%time_relx/)
    
    if (allocated(f_data_relax%sim(1)%sim_setup%stprows)) then
        deallocate(f_data_relax%sim(1)%sim_setup%stprows)
    endif
    allocate(f_data_relax%sim(1)%sim_setup%stprows(2))
    f_data_relax%sim(1)%sim_setup%stprows = (/1,2/)
    
    if (allocated(f_data_relax%sim(1)%sim_setup%steps)) then
        deallocate(f_data_relax%sim(1)%sim_setup%steps)
    endif
    allocate(f_data_relax%sim(1)%sim_setup%steps(1))
    f_data_relax%sim(1)%sim_setup%steps(1) = 1.d0
    
    if (allocated(f_data_relax%sim(1)%sim_setup%disp_error_scale)) then
        deallocate(f_data_relax%sim(1)%sim_setup%disp_error_scale)    
        deallocate(f_data_relax%sim(1)%sim_setup%load_error_scale)    ! Assume if disp isn't allocated neither is load
    endif
    allocate(f_data_relax%sim(1)%sim_setup%disp_error_scale(4,1))
    allocate(f_data_relax%sim(1)%sim_setup%load_error_scale(4,1))
    f_data_relax%sim(1)%sim_setup%disp_error_scale = 1.d0
    f_data_relax%sim(1)%sim_setup%load_error_scale = 1.d0
    
    ! Control values
    if (allocated(f_data_relax%sim(1)%exp%ctrl)) then
        deallocate(f_data_relax%sim(1)%exp%ctrl)
    endif
    allocate(f_data_relax%sim(1)%exp%ctrl(5,1))
    f_data_relax%sim(1)%exp%ctrl(:,1) = (/1, -3, -3, -3, -3/)   ! Load ctrl, same end value as experiment
    
    ! Time incrementation values
    deallocate(f_data_relax%sim(1)%iter%time_incr)
    allocate(f_data_relax%sim(1)%iter%time_incr(5,1))
    f_data_relax%sim(1)%iter%time_incr = f_data%sim(simnr_new)%iter%time_incr
    
    ! Set to be continued from "itself"
    f_data_relax%sim(1)%init%cont_analysis = 1
    
    ! Set new node positions if input
    if (present(rpos)) then
        f_data_relax%sim(1)%mesh1d%node_pos = (/rpos(1,:), rpos(size(rpos,1),size(rpos,2))/)
    endif
    
    ! If sample is zero, don't allow changing the inner radius
    if (f_data_relax%sim(1)%mesh1d%node_pos(1)<1.d-12) then
        f_data_relax%sim(1)%exp%ctrl(4,1) = 1
    endif
    
    
    call atp_simulate(error, props, f_data_relax, 1)
    if (error>(huge(1.d0)/10)) then
        converged = .false.
    else
        converged = .true.
    endif
    

end subroutine

subroutine remesh(node_pos, f_data, f_data_old, simnr, rpos, h0, gp_stress, gp_strain, gp_F, gp_sv, u0, is_solid)
implicit none
    double precision, intent(in)    :: node_pos(:)  ! New (deformed) node positions
    type(fdata_typ), intent(in)    :: f_data, f_data_old
    integer, intent(in)             :: simnr
    double precision, intent(inout) :: rpos(:,:)
    double precision, intent(out)   :: h0
    double precision, intent(inout) :: gp_stress(:,:,:), gp_strain(:,:,:), gp_F(:,:,:), gp_sv(:,:,:)
    double precision, intent(inout) :: u0(:)
    logical, intent(in)             :: is_solid !If the remeshed part should be solid (i.e. ri=0)
    
    integer                         :: nel, ngp, nnod, nnod_tot
    integer                         :: nel_old, ngp_old, nnod_old, nnod_tot_old
    double precision, allocatable   :: all_node_pos0(:)     ! Undeformed locations of all nodes, including internal element nodes (if any)
    double precision, allocatable   :: all_node_pos0_old(:) ! Same as all_node_pos0, but for previous mesh
    double precision, allocatable   :: u_old(:)             ! Old displacement (from f_data_old%end_res%u_end)
    double precision, allocatable   :: gp_pos0(:)           ! Undeformed gauss point positions
    double precision, allocatable   :: gp_pos0_old(:)       ! Undeformed gauss point positions for previous mesh
    
    double precision                :: mult_factor          ! Multiplication factor for interpolation
    
    double precision                :: r, dr_dr0, eps0(6), eps_old0(6)
    double precision, allocatable   :: dNpdR(:), Np(:), ue_r(:)
    double precision, allocatable   :: Bs_mat(:,:)
    integer, allocatable            :: dof(:)
    
    integer                         :: ind, ind_old
    integer                         :: iel, igp, inod       ! Counters for element number and gauss point number
    
    ! Writing of results
    integer                         :: res_fid              ! Result file fid
    character(len=30)               :: format_spec          ! Writing format
    
    h0          = f_data_old%sim(1)%end_res%h0_true  ! Set height to actual from old simulation (no need to adjust rotation or axial displacement values)
    
    nel_old     = size(f_data_old%sim(1)%mesh1d%node_pos) - 1
    nel         = size(f_data%sim(simnr)%mesh1d%node_pos) - 1
    nnod_old    = f_data_old%sim(1)%mesh1d%element_order + 1
    nnod        = f_data%sim(simnr)%mesh1d%element_order + 1
    ngp_old     = f_data_old%sim(1)%mesh1d%ngp
    ngp         = f_data%sim(simnr)%mesh1d%ngp
    nnod_tot    = nel*(nnod-1)+1
    nnod_tot_old= nel_old*(nnod_old-1)+1
    allocate(dof(nnod), ue_r(nnod), dNpdR(nnod), Np(nnod), Bs_mat(6,nnod+2))
    
    
    allocate(u_old(nnod_tot_old+2))
    u_old = f_data_old%sim(1)%end_res%u_end
    ! Get node coordinates
    allocate(all_node_pos0_old(nnod_tot_old))
    allocate(all_node_pos0(nnod_tot))
    all_node_pos0_old(1:(nnod_tot_old):(nnod_old-1)) = f_data_old%sim(1)%mesh1d%node_pos
   
    if (nnod_old==3) then
        ! Need to assign internal element nodal coordinates
        all_node_pos0_old(2:(nnod_tot_old-1):2) = (all_node_pos0_old(1:(nnod_tot_old-2):2) &
                                                +  all_node_pos0_old(3:(nnod_tot_old):2))/2.d0
    endif
    do ind=1,size(node_pos)
        all_node_pos0((ind-1)*(nnod-1)+1) = get_orig_node_pos(node_pos(ind), all_node_pos0_old, u_old)
    enddo
    if (nnod==3) then
        ! Need to assign internal element nodal coordinates
        all_node_pos0(2:(nnod_tot-1):2) = (all_node_pos0(1:(nnod_tot-2):2) &
                                        +  all_node_pos0(3:(nnod_tot):2))/2.d0
    endif
    
    if (is_solid) then
        all_node_pos0(1) = 0.d0
    endif
    
    ! Calculate the new rpos:
    do iel=1,nel
        rpos(:,iel) = all_node_pos0((1+(iel-1)*(nnod-1)):(1+(iel)*(nnod-1)))
    enddo
    
    ! Interpolate nodal values
    ! Check that no nodes are outside the initial geometry (we don't allow additative manufacturing...)
    !if (all_node_pos0(1)<(all_node_pos0_old(1)-1e-10)) then
    !    call write_output('Inside node after element removal cannot be outside the original deformed body', 'error', 'sim:atp_er', halt=.false.)
    !    call write_output('Distance = '//dbl2str(all_node_pos0(1)-all_node_pos0_old(1)), 'error', 'sim:atp_er', loc=.false.)
    !endif
    !
    !if (all_node_pos0(nnod_tot)>(all_node_pos0_old(nnod_tot_old)+1e-10)) then
    !    call write_output('Outside node after element removal cannot be outside the original deformed body', 'error', 'sim:atp_er', halt=.false.)
    !    call write_output('Distance = '//dbl2str(all_node_pos0(nnod_tot)-all_node_pos0_old(nnod_tot_old)), 'error', 'sim:atp_er', loc=.false.)
    !endif
    
    ind_old = 2
    do ind=1,nnod_tot
        do while ((all_node_pos0(ind)>all_node_pos0_old(ind_old)).and.((ind_old+1)<size(all_node_pos0_old)))
            ind_old = ind_old + 1
        enddo
        mult_factor = (all_node_pos0(ind)-all_node_pos0_old(ind_old-1))/(all_node_pos0_old(ind_old)-all_node_pos0_old(ind_old-1))
        u0(2+ind) = u_old((ind_old+2)-1) + (u_old(ind_old+2)-u_old((ind_old+2)-1))*mult_factor
    enddo
    
    ! Set the common dofs
    u0(1:2) = u_old(1:2)
    
    ! Interpolate the gauss point values
    
    ! Get gauss point coordinates for old mesh
    allocate(gp_pos0_old(ngp_old*nel_old), gp_pos0(ngp*nel))
    gp_pos0_old = get_gp_coords(nel_old, ngp_old, nnod_old, all_node_pos0_old)
    gp_pos0     = get_gp_coords(nel, ngp, nnod, all_node_pos0)
    
    ! Loop over all new gauss points and find new values
    ! "State type variables" (State variables and stress are interpolated)
    ! "Kinematic variables" (Strain and deformation gradient are calculated from interpolated displacements)
    ind = 1
    ind_old = 2
    call element_setup(ngp, nnod, .false.)  ! New element 
    do iel = 1,nel
        do inod=1,nnod
            dof(inod) = (inod+2) + (nnod-1)*(iel-1)
        enddo
        ue_r = u0(dof)
        
        do igp = 1,ngp
            ! Interpolation
            do while ((gp_pos0(ind)>gp_pos0_old(ind_old)).and.((ind_old+1)<size(gp_pos0_old)))
                ind_old = ind_old + 1
            enddo
            mult_factor = (gp_pos0(ind)-gp_pos0_old(ind_old-1))/(gp_pos0_old(ind_old)-gp_pos0_old(ind_old-1))
            gp_stress(:,igp,iel) = f_data_old%sim(1)%end_res%stress_end(ind_old-1, :) &
                + (f_data_old%sim(1)%end_res%stress_end(ind_old,:)-f_data_old%sim(1)%end_res%stress_end(ind_old-1,:))*mult_factor
            gp_sv(:,igp,iel) = f_data_old%sim(1)%end_res%statev_end(ind_old-1, :) &
                + (f_data_old%sim(1)%end_res%statev_end(ind_old,:)-f_data_old%sim(1)%end_res%statev_end(ind_old-1,:))*mult_factor
            
            !Old code for interpolating even kinematic values
            !gp_strain(:,igp,iel) = f_data_old%sim(1)%end_res%strain_end(ind_old-1, :) &
            !   + (f_data_old%sim(1)%end_res%strain_end(ind_old,:)-f_data_old%sim(1)%end_res%strain_end(ind_old-1,:))*mult_factor
            !gp_F(:,igp,iel) = f_data_old%sim(1)%end_res%dfgrd_end(ind_old-1, :) &
            !   + (f_data_old%sim(1)%end_res%dfgrd_end(ind_old,:)-f_data_old%sim(1)%end_res%dfgrd_end(ind_old-1,:))*mult_factor
            
            ! Kinematic calculations
            eps_old0 = 0.d0
            call shapefun(gp_pos0(ind), rpos(:,iel), dNpdR, Np, ue_r, r, dr_dr0)
            call FBeps_calc(gp_F(:,igp,iel), Bs_mat, eps0, gp_strain(:,igp,iel), dNpdR, Np, dr_dr0, r, u0(2), u0(1), gp_pos0(ind), h0, eps_old0)
            
            ind = ind + 1
        enddo
        
    enddo

    
    ! If requested, write output in forms of old and new values as functions of the undeformed coordinates
    if (f_data%glob%resnr>0) then
        open(newunit=res_fid, file=trim(f_data%glob%outname)//'_sim'//int2str(simnr)//'_'//int2str(f_data%glob%resnr)//'.txt', status='replace')
        write(res_fid, "(A)") 'Result file for interpolation during remeshing'
        write(res_fid, *) ''    ! Extra line break
        
        ! Print old nodal values
        format_spec = '('//int2str(size(all_node_pos0_old))//'ES15.5E3)'
        write(res_fid, "(A)") 'Old nodal values (rows = node positions, displacements)'
        write(res_fid, format_spec) all_node_pos0_old
        write(res_fid, format_spec) u_old(3:)
        write(res_fid, *) ''    ! Line break
        
        ! Print new nodal values
        format_spec = '('//int2str(size(all_node_pos0))//'ES15.5E3)'
        write(res_fid, "(A)") 'New nodal values (rows = node positions, displacements)'
        write(res_fid, format_spec) all_node_pos0
        write(res_fid, format_spec) u0(3:)
        write(res_fid, *) ''    ! Line break
        
        ! Print the old gauss point values
        format_spec = '('//int2str(size(gp_pos0_old))//'ES15.5E3)'
        write(res_fid, "(A)") 'Old gauss point values (rows = node positions, stress(6), '
        write(res_fid, "(A)") 'strain(6), deformation gradient(9) and state variables ('//int2str(f_data%glob%nstatv)//'))'
        write(res_fid, format_spec) gp_pos0_old
        write(res_fid, *) ''    ! Line break
        write(res_fid, format_spec) f_data_old%sim(1)%end_res%stress_end
        write(res_fid, *) ''    ! Line break
        write(res_fid, format_spec) f_data_old%sim(1)%end_res%strain_end
        write(res_fid, *) ''    ! Line break
        write(res_fid, format_spec) f_data_old%sim(1)%end_res%dfgrd_end
        write(res_fid, *) ''    ! Line break
        write(res_fid, format_spec) f_data_old%sim(1)%end_res%statev_end
        write(res_fid, *) ''    ! Line break
        
        ! Print the new gauss point values
        format_spec = '('//int2str(size(gp_pos0))//'ES15.5E3)'
        write(res_fid, "(A)") 'New gauss point values (rows = node positions, stress(6), '
        write(res_fid, "(A)") 'strain(6), deformation gradient(9) and state variables ('//int2str(f_data%glob%nstatv)//'))'
        write(res_fid, format_spec) gp_pos0
        write(res_fid, *) ''    ! Line break
        write(res_fid, format_spec) transpose(reshape(gp_stress, (/6, ngp*nel/)))
        write(res_fid, *) ''    ! Line break
        write(res_fid, format_spec) transpose(reshape(gp_strain, (/6, ngp*nel/)))
        write(res_fid, *) ''    ! Line break
        write(res_fid, format_spec) transpose(reshape(gp_F, (/9, ngp*nel/)))
        write(res_fid, *) ''    ! Line break
        write(res_fid, format_spec) transpose(reshape(gp_sv, (/f_data%glob%nstatv, ngp*nel/)))
        write(res_fid, *) ''    ! Line break
        
        close(res_fid)
    endif
    
    
end subroutine
        
function get_gp_coords(nel, ngp, nnod, all_node_pos0) result(gp_pos0)
implicit none
    integer, intent(in)             :: nel, ngp, nnod
    double precision, intent(in)    :: all_node_pos0(:)
    double precision, allocatable   :: gp_pos0(:)
    integer                         :: ind, iel, igp
    double precision                :: rpos_el(nnod)
    double precision                :: gp_weight    ! Unused variable, just needed for passing to get_gpinfo
    
    call element_setup(ngp, nnod, .false.)
    allocate(gp_pos0(ngp*nel))
    ind = 1
    do iel = 1,nel
        do igp = 1,ngp
            rpos_el = all_node_pos0( (1+(iel-1)*(nnod-1)):(1+(iel)*(nnod-1)))
            call get_gpinfo(gp_pos0(ind), gp_weight, igp, rpos_el)
            ind = ind + 1
        enddo
    enddo
    
end function

function get_orig_node_pos(rhat, rold, uold) result(rnew)
! Find initial deformed node position of the deformed location rhat
implicit none
    double precision, intent(in)        :: rhat     !Deformed node location
    double precision, intent(in)        :: rold(:)  !The old nodal positions (old mesh)
    double precision, intent(in)        :: uold(:)  !The old displacements (old mesh)
    double precision                    :: rnew     !The new (undeformed/initial) node position of rhat
    integer                             :: pos
    double precision                    :: unew
    double precision                    :: du_dr
    integer                             :: k1, maxiter
    double precision                    :: tol
    double precision                    :: residual
    
    maxiter = 10
    tol = 1.d-10
    
    pos = 2
    rnew = rhat
    residual = 2*tol
    k1 = 1
    do while (residual>tol)
        if (k1>maxiter) then
            call write_output('Cannot find node position', 'error', 'sim:atp_element_removal')
        endif
        
        do while(pos<size(rold))
            pos = pos + 1
            if (rold(pos)>rnew) then
                exit
            endif
        enddo
        
        do while(pos>2)
            pos = pos -1
            if (rold(pos-1)<rnew) then
                exit
            endif
        enddo

        du_dr = (uold(pos)-uold(pos-1))/(rold(pos)-rold(pos-1))
        unew = uold(pos) + du_dr*(rnew-rold(pos-1))
        residual = rnew+unew-rhat
        rnew = rnew - (residual)/(1.d0+du_dr)
        k1 = k1 + 1
    enddo
    
end function

subroutine set_end_values(sim_in, sim_out)
implicit none
    type(sim_typ), intent(in)       :: sim_in
    type(sim_typ), intent(inout)    :: sim_out
    integer                         :: ngp_total, nstatv
    
    ngp_total = size(sim_in%end_res%stress_end, 1)
    nstatv = size(sim_in%end_res%statev_end, 1)
    
        ! Allocate if needed
    if(.not.allocated(sim_out%end_res%u_end)) then
        allocate(sim_out%end_res%u_end(size(sim_in%end_res%u_end)))
    endif
    
    if(.not.allocated(sim_out%end_res%stress_end)) then
        allocate(sim_out%end_res%stress_end(ngp_total, 6))
    endif
    
    if(.not.allocated(sim_out%end_res%strain_end)) then
        allocate(sim_out%end_res%strain_end(ngp_total, 6))
    endif
    
    if(.not.allocated(sim_out%end_res%dfgrd_end)) then
        allocate(sim_out%end_res%dfgrd_end(ngp_total, 9))
    endif
    
    if(.not.allocated(sim_out%end_res%statev_end)) then
        allocate(sim_out%end_res%statev_end(ngp_total, nstatv))
    endif
    
    sim_out%end_res%u_end       = sim_in%end_res%u_end
    sim_out%end_res%stress_end  = sim_in%end_res%stress_end
    sim_out%end_res%strain_end  = sim_in%end_res%strain_end
    sim_out%end_res%dfgrd_end   = sim_in%end_res%dfgrd_end
    sim_out%end_res%statev_end  = sim_in%end_res%statev_end
    
    sim_out%end_res%disp_end    = sim_in%end_res%disp_end
    sim_out%end_res%load_end    = sim_in%end_res%load_end
    sim_out%end_res%h0_true     = sim_in%end_res%h0_true
end subroutine

end module 

