
module atp_element_removal_mod
use atp_mod
use types_mod
use atp_element_mod
use output_mod
implicit none
    private
    public  :: atp_element_removal ! Axial-torsion-pressure element removal routine
    
    logical :: interpolate_failed = .false.
    
    contains
! Element removal procedure
! ============================================================================================================
! 	1. Import data from continued analysis
!	2. Solve load = 0.d0 (I.e. relax)
!	3. Remesh with new mesh nodes as specified for current simulation
!		- Interpolate radial displacements, stress and state variables to new mesh
!       - Based on interpolated displacements, calculate new strain and deformation gradient
!       - Output if requested old values and interpolated values
!	4. Solve load=0 (I.e. relax)
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
    type(fdata_typ), intent(inout)  :: f_data           ! Settings input
    integer, intent(in)             :: simnr            ! Simulation number (=> f_data%sim(simnr) is pertinent to current simulation)

    ! Previous analysis relaxation
    integer                         :: prev_simnr           ! Previous simulation number
    type(fdata_typ)                 :: f_data_initial_relax ! Data belonging to relaxation simulation

    ! Node position iterations
    double precision, allocatable   :: node_disp_old_mesh(:), node_pos_old_mesh(:)
    integer                         :: ntnod_old, nenod_old, nel_old
    double precision, allocatable   :: node_pos_guess(:), node_pos_result(:)            ! Node positions 
    double precision, allocatable   :: node_pos_guess_old(:), node_pos_result_old(:)    ! Node positions for better guess
    double precision, allocatable   :: node_pos_update(:)                               ! Node position update
    integer, allocatable            :: node_disp_ind(:)                                 ! Indices in dofs belonging to node_pos variables
    integer                         :: k1                                               ! Iteration index
    logical                         :: geom_conv                                        ! Convergence indicator for geometry iterations
    logical                         :: relx_conv                                        ! Convergence indicator for relaxation
    double precision                :: pos_error                                        ! Position error

    ! Problem size convenience variables
    integer                         :: nel, ngp, nenod, ntnod, ndof
    logical                         :: new_is_solid                                     ! Bool to decide if the new is a solid specimen and 
                                                                                        ! inner radius should be kept exactly at 0.d0
                                                                                        
    ! Variables for solve_incr
    type(fdata_typ)                 :: f_data_relax
    double precision, allocatable   :: rpos(:,:)
    double precision                :: h0, load(4), disp(4), time(2)
    integer, allocatable            :: free_dofs(:), known_dofs(:)
    double precision, allocatable   :: gp_F(:,:,:), gp_stress(:,:,:), gp_sv(:,:,:), gp_strain(:,:,:), u0(:), du(:)
    double precision                :: iter_tol
    integer                         :: niter, iter_max
    logical                         :: lconv
    double precision                :: pnewdt
    double precision, allocatable   :: iter_err_norm(:)

    ! Simulation timing
    double precision                :: starttime, stoptime
    integer                         :: clock_count, clock_rate
    
    ! Error tracking
    double precision, allocatable   :: geom_error(:)
    integer                         :: res_fid
    
    allocate(geom_error(f_data%sim(simnr)%atp_er%geom_iter_max))
    
    ! Simulation timing
    call system_clock ( clock_count, clock_rate)
    starttime = real(clock_count)/real(clock_rate)

    ! Convenience size variables
    nel = size(f_data%sim(simnr)%mesh1d%node_pos)-1
    ngp = f_data%sim(simnr)%mesh1d%ngp
    nenod= f_data%sim(simnr)%mesh1d%element_order + 1
    ntnod = nel*(nenod-1)+1
    ndof= 2 + nel*(nenod-1)+1

    ! == RELAX PREVIOUS SIMULATION == 
    prev_simnr = abs(f_data%sim(simnr)%init%cont_analysis)
    call atp_relax(relx_conv, props, f_data, simnr, prev_simnr, f_data_initial_relax, '_relax'//int2str(simnr)//'_0')
    if (.not.relx_conv) then
        error = huge(1.d0)
        return
    endif
    

    ! == SETUP NODE POSITION ITERATIONS ==
    allocate(node_disp_ind(nel+1))
    ! Index in dof for specified node displacement values
    node_disp_ind = (/(k1, k1=3,ndof,(nenod-1))/)
    allocate(node_pos_guess(nel+1), node_pos_result(nel+1))
    allocate(node_pos_guess_old(nel+1), node_pos_result_old(nel+1))
    allocate(node_pos_update(nel+1))
    
    ! Guess for undeformed nodal positions
    node_pos_guess = f_data%sim(simnr)%mesh1d%node_pos  ! Guess for deformed positions, needs to be changed to undeformed:
     ! Get dimensions for old mesh
    nenod_old = f_data_initial_relax%sim(1)%mesh1d%element_order + 1
    nel_old = (size(f_data_initial_relax%sim(1)%mesh1d%node_pos) - 1)
    ntnod_old = nel_old*(nenod_old-1)+1
    allocate(node_disp_old_mesh(nel_old+1), node_pos_old_mesh(nel_old+1))
     ! Find interpolated displacements by interpolating on the deformed configuration as an initial guess
    node_disp_old_mesh = f_data_initial_relax%sim(1)%end_res%u_end(3:(ntnod_old+2):(nenod_old-1))
    node_pos_old_mesh = f_data_initial_relax%sim(1)%mesh1d%node_pos_undef
    node_pos_update = interpolate(node_pos_old_mesh + node_disp_old_mesh, node_disp_old_mesh, node_pos_guess)
    
    if (interpolate_failed) then
        call write_output('Interpolation failed on initial guess', 'status', 'atp:element_removal')
        error = huge(1.d0)
        return
    endif
    
     ! Update the guess to be of the undeformed positions
    node_pos_guess = node_pos_guess - node_pos_update ! Guess for undeformed positions
    node_pos_result= 0.d0

    ! == SETUP NECESSARY VARIABLES FOR remesh AND solve_incr == 
    allocate(rpos(nenod, nel))       ! Undeformed nodal positions
    ! Gauss point values:
    allocate(gp_stress(6,ngp,nel), gp_sv(f_data%glob%nstatv,ngp,nel), gp_strain(6,ngp,nel), gp_F(9,ngp,nel))
    allocate(u0(ndof), du(ndof), iter_err_norm(ndof))   !Nodal values (and normalization for error (not used))
    
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
        
        call remesh(node_pos_guess, f_data, f_data_initial_relax, simnr, rpos, h0, gp_stress, gp_strain, gp_F, gp_sv, u0)
        if (interpolate_failed) then
            write(*,*) node_pos_guess
            write(*,*) node_pos_old_mesh
            error = huge(1.d0)
            return
        endif
        
        if (new_is_solid) then
            u0(3) = 0.d0
        endif
        
        ! Solve zero displacement increment
        call element_setup(ngp, nenod, .false., -simnr)
        allocate(free_dofs(1))
        free_dofs   = -1    ! All displacements prescribed:
        call solve_incr(rpos, h0, load, disp, f_data%sim(simnr)%init%temp_init, 0.d0, time, f_data%sim(simnr)%atp_er%time_remesh, free_dofs, &
               free_dofs, gp_F, gp_stress, gp_sv, gp_strain, u0, du, f_data%sim(1)%iter, niter, lconv, pnewdt, &
               k1, -simnr, props, f_data%glob%cmname, f_data%glob%umat_address, f_data%glob%nlgeom, iter_err_norm)
    
        if (.not.lconv) then
           call write_output('No convergence for local material problem after element removal', 'status', 'sim:atp')
           exit GEOM_ITER_LOOP
        endif
        
        ! Solve zero external displacement increment
        deallocate(free_dofs)
        call gen_free_dofs(free_dofs, known_dofs, ndof, [1,1,1,1], .false.)
        call solve_incr(rpos, h0, load, disp, f_data%sim(simnr)%init%temp_init, 0.d0, time, f_data%sim(simnr)%atp_er%time_remesh, free_dofs, &
               known_dofs, gp_F, gp_stress, gp_sv, gp_strain, u0, du, f_data%sim(1)%iter, niter, lconv, pnewdt, &
               k1, -simnr, props, f_data%glob%cmname, f_data%glob%umat_address, f_data%glob%nlgeom, iter_err_norm)
        
        deallocate(free_dofs, known_dofs)
        if (.not.lconv) then
           call write_output('No convergence for for internal relaxation after element removal', 'status', 'sim:atp')
           exit GEOM_ITER_LOOP
        endif
   
        ! Set end values in f_data%sim(simnr)        
        if (new_is_solid) then  ! Set disp(3) = 0
            disp = (/u0(1)/h0, u0(2), u0(3)*0.0, u0(size(u0))/rpos(size(rpos,1), size(rpos,2))/)
        else
            disp = (/u0(1)/h0, u0(2), u0(3)/rpos(1,1), u0(size(u0))/rpos(size(rpos,1), size(rpos,2))/)
        endif
        call atp_export_end(f_data, simnr, u0, gp_stress, gp_strain, gp_F, gp_sv, disp, load, h0)
   
        ! Relax 
        call atp_relax(relx_conv, props, f_data, simnr, simnr, f_data_relax, '_relax'//int2str(simnr)//'_'//int2str(k1), rpos, h0)
        if (.not.relx_conv) then
            error = huge(1.d0)
            return
        endif
   
        node_pos_result = (/rpos(1,:), rpos(size(rpos,1), size(rpos,2))/) + f_data_relax%sim(1)%end_res%u_end(node_disp_ind)
        pos_error = sqrt(  sum( (node_pos_result-f_data%sim(simnr)%mesh1d%node_pos)**2 &
                                  /max(1.d-6,f_data%sim(simnr)%mesh1d%node_pos)**2 )  )
        geom_error(k1) = pos_error
        
        if (pos_error<f_data%sim(simnr)%atp_er%node_pos_tol) then
            geom_conv = .true.
            error = 0.d0
            exit GEOM_ITER_LOOP
        elseif (k1==1) then
            ! First iteration - fix point update
            node_pos_guess_old = node_pos_guess
            node_pos_result_old = node_pos_result
            node_pos_guess = node_pos_guess + (f_data%sim(simnr)%mesh1d%node_pos-node_pos_result)
        else
            ! Further iteration, secant update?
            node_pos_update = (f_data%sim(simnr)%mesh1d%node_pos-node_pos_result) &
                                *(node_pos_guess-node_pos_guess_old) &
                                /sign(max(abs(node_pos_result-node_pos_result_old),f_data%sim(simnr)%mesh1d%node_pos(nel+1)*1.d-20), node_pos_result-node_pos_result_old)
            
            node_pos_guess_old = node_pos_guess
            node_pos_result_old = node_pos_result
            node_pos_guess = node_pos_guess + node_pos_update
        endif
        if (new_is_solid) then ! Force inner node to be zero
                node_pos_guess(1) = 0.d0
        endif
        
    enddo GEOM_ITER_LOOP
    
    if (f_data%glob%resnr>0) then
        open(newunit=res_fid, file=trim(f_data%glob%outname)//'_sim'//int2str(simnr)//'_'//int2str(f_data%glob%resnr)//'.txt', position='append')
        write(res_fid, '(A)') 'Nodal positioning relative error evolution:'
        write(res_fid, '('//int2str(k1)//'ES15.5)') geom_error(1:k1)
        close(res_fid)
    endif
    
    if (.not.geom_conv) then
        if (lconv) then
            call write_output('No convergence for geometry update due to element removal, error=huge', 'status', 'sim:atp')
        endif
        
        error = huge(1.d0)
    else
        call set_end_values(f_data_relax%sim(1), f_data%sim(simnr))
    endif
    

    call system_clock ( clock_count, clock_rate)
    stoptime = real(clock_count)/real(clock_rate)

end subroutine atp_element_removal

subroutine atp_relax(converged, props, f_data, simnr_new, simnr_old, f_data_relax, append_result_name, rpos, h0)
! Solve relaxing before element removal by setting up another simulation
implicit none
    logical, intent(out)            :: converged
    double precision, intent(in)    :: props(:)
    type(fdata_typ), intent(in)     :: f_data
    integer, intent(in)             :: simnr_new, simnr_old
    type(fdata_typ), intent(out)    :: f_data_relax
    character(len=*), intent(in)    :: append_result_name
    double precision, optional      :: rpos(:,:)
    double precision, optional      :: h0
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
        if (present(h0)) then
            f_data_relax%sim(1)%mesh1d%h0 = h0
        else
            call write_output('h0 must be given if rpos is given', 'status', 'atp:element_removal:relax')
        endif       
    elseif (f_data%sim(simnr_old)%init%cont_analysis<0) then
        ! if rpos (and h0) is not present, this is the first relax. 
        ! If we relax from a continued analysis we must account for that the mesh is given for the deformed configuration
        ! Otherwise, if rpos is given we have the correct (initial) mesh (including correct (initial) h0)
        f_data_relax%sim(1)%init%cont_analysis = -1 
    endif
    
    ! If sample is zero, don't allow changing the inner radius
    if (f_data_relax%sim(1)%mesh1d%node_pos(1)<1.d-12) then
        f_data_relax%sim(1)%exp%ctrl(4,1) = 1
    endif
    
    ! Reset error count to avoid wrong dimensioning of error vector sizes
    f_data_relax%sim(1)%err%hist_rows = -1
    
    call atp_simulate(error, props, f_data_relax, 1)
    if (error>(huge(1.d0)/10)) then
        converged = .false.
    else
        converged = .true.
    endif
    

end subroutine

function get_gp_coords(nel, ngp, nnod, r0) result(gp_pos0)
implicit none
    integer, intent(in)             :: nel, ngp, nnod
    double precision, intent(in)    :: r0(:)
    double precision, allocatable   :: gp_pos0(:)
    integer                         :: ind, iel, igp
    double precision                :: rpos_el(nnod)
    double precision                :: gp_weight    ! Unused variable, just needed for passing to get_gpinfo
    
    call element_setup(ngp, nnod, .false., 0)
    allocate(gp_pos0(ngp*nel))
    ind = 1
    do iel = 1,nel
        do igp = 1,ngp
            rpos_el = r0( (1+(iel-1)*(nnod-1)):(1+(iel)*(nnod-1)))
            call get_gpinfo(gp_pos0(ind), gp_weight, igp, rpos_el)
            ind = ind + 1
        enddo
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

subroutine remesh(r0outer, f_data, f_data_old, simnr, rpos, h0, gp_stress, gp_strain, gp_F, gp_sv, u0)
! Given a guess of initial new nodal positions r0, determine
! u0         Nodal displacements (including axial and rotation on positions 1 and 2), interpolated 
! gp_stress  Stress at gauss points (interpolated)
! gp_sv      State variables at gauss points (interpolated)
! gp_strain  Strain at gauss points (calculated based on interpolated nodal displacements u0)
! gp_F       Deformation gradient at gauss points (calculated based on interpolated nodal displacements u0)
implicit none
    double precision, intent(in)    :: r0outer(:)  ! New (undeformed) node positions. Does not contain potential middle nodes for 2nd order elements
    type(fdata_typ), intent(in)     :: f_data, f_data_old
    integer, intent(in)             :: simnr
    double precision, intent(out)   :: rpos(:,:)
    double precision, intent(out)   :: h0
    double precision, intent(out)   :: gp_stress(:,:,:), gp_strain(:,:,:), gp_F(:,:,:), gp_sv(:,:,:)
    double precision, intent(out)   :: u0(:)
    
    integer                         :: iel, nel, nel_old     ! Element counter and number of elements
    integer                         :: inod, nenod, ntnod    ! Nodal counter, nodes per element and number of total nodes (length of r0)
    integer                         :: nenod_old, ntnod_old  ! Nodes per element (old mesh) and total number of nodes in old mesh
    integer                         :: igp, ngp, ngp_old     ! Gauss point counter, number of gauss point (new mesh) and number of gauss points (old mesh) (per element)
    integer                         :: ndof, ndof_old        ! Number of degrees of freedom (new and old mesh)
    integer                         :: ind, ind_old          ! Counter for values att all gauss points (not only inside the element)
    
    double precision, allocatable   :: r0(:)                ! New (undeformed) node position. Contains also the middle nodes (calculated from r0outer)
    double precision                :: mult_factor          ! Factor used in interpolation
    double precision, allocatable   :: r0_old(:)            ! Old nodal positions
    double precision, allocatable   :: u0_old(:)            ! Displacements on old mesh
    double precision, allocatable   :: gp_pos0(:)           ! Undeformed gauss point positions
    double precision, allocatable   :: gp_pos0_old(:)       ! Undeformed gauss point positions for old mesh
    
    double precision                :: rgp, dr_dr0, eps0(6), eps_old0(6)
    double precision, allocatable   :: dNpdR(:), Np(:), ue_r(:)
    double precision, allocatable   :: Bs_mat(:,:)
    integer, allocatable            :: dof(:)
    double precision                :: hold                 ! Height (gauge length) for old simulation
    
    integer                         :: res_fid              ! Result file fid
    character(len=30)               :: format_spec          ! Writing format
    
    ! Get sizes of old and new mesh
    nel_old     = size(f_data_old%sim(1)%mesh1d%node_pos) - 1
    nel         = size(f_data%sim(simnr)%mesh1d%node_pos) - 1
    nenod_old   = f_data_old%sim(1)%mesh1d%element_order + 1
    nenod       = f_data%sim(simnr)%mesh1d%element_order + 1
    ngp_old     = f_data_old%sim(1)%mesh1d%ngp
    ngp         = f_data%sim(simnr)%mesh1d%ngp
    ntnod       = nel*(nenod-1)+1
    ntnod_old   = nel_old*(nenod_old-1)+1
    ndof        = ntnod + 2
    ndof_old    = ntnod_old + 2
    
    
    ! Calculate the full r0 from r0outer
    allocate(r0(ntnod))
    r0(1:ntnod:(nenod-1)) = r0outer
    if (nenod==3) then
        r0(2:(ntnod-1):2) = (r0(1:(ntnod-2):2) + r0(3:ntnod:2))/2.d0
    endif
        
    ! Calculate the new rpos:
    do iel=1,nel
        rpos(:,iel) = r0((1+(iel-1)*(nenod-1)):(1+(iel)*(nenod-1)))
    enddo
    
    ! Get old node coordinates
    allocate(r0_old(ntnod_old))
    if (f_data_old%sim(1)%init%cont_analysis<0) then    
        ! Need to treat if the old simulation is continued<0 such that node_pos are not the initial positions!
        r0_old(1:(ntnod_old):(nenod_old-1)) = f_data_old%sim(1)%mesh1d%node_pos - f_data_old%sim(1)%end_res%u_end(3:ndof_old:(nenod_old-1))
        if (nenod_old==3) then   ! Need to assign internal element nodal coordinates
            r0_old(2:(ntnod_old-1):2) = (r0_old(1:(ntnod_old-2):2) + r0_old(3:(ntnod_old):2))/2.d0 - f_data_old%sim(1)%end_res%u_end(4:(ndof_old-1):2)
        endif
    else
        r0_old(1:(ntnod_old):(nenod_old-1)) = f_data_old%sim(1)%mesh1d%node_pos
        if (nenod_old==3) then   ! Need to assign internal element nodal coordinates
            r0_old(2:(ntnod_old-1):2) = (r0_old(1:(ntnod_old-2):2) + r0_old(3:(ntnod_old):2))/2.d0
        endif
    endif
    
    ! Interpolate the displacements and get height
    allocate(u0_old(ntnod_old+2))
    u0_old = f_data_old%sim(1)%end_res%u_end
    u0(3:ndof) = interpolate(r0_old, u0_old(3:ndof_old), r0)
    hold       = f_data_old%sim(1)%end_res%h0_true ! Height for which old displacements given
    h0         = f_data%sim(simnr)%mesh1d%h0/(1 + u0_old(1)/hold)
    u0(1:2)    = u0_old(1:2)*h0/hold
    
    ! Calculate the gauss point coordinates
    allocate(gp_pos0_old(ngp_old*nel_old), gp_pos0(ngp*nel))
    gp_pos0_old = get_gp_coords(nel_old, ngp_old, nenod_old, r0_old)
    gp_pos0     = get_gp_coords(nel, ngp, nenod, r0)
    
    ! Loop over all new gauss points and find new values
    ! "State type variables" (State variables and stress are interpolated)
    ! "Kinematic variables" (Strain and deformation gradient are calculated from interpolated displacements)
    allocate(dof(nenod), ue_r(nenod), dNpdR(nenod), Np(nenod), Bs_mat(6,nenod+2))
    ind = 1
    ind_old = 2
    call element_setup(ngp, nenod, .false., simnr)  ! New element 
    do iel = 1,nel
        do inod=1,nenod
            dof(inod) = (inod+2) + (nenod-1)*(iel-1)
        enddo
        ue_r = u0(dof)
        
        do igp = 1,ngp
            ! Interpolation (not in function for efficiency, same multiplication factor)
            do while ((gp_pos0(ind)>gp_pos0_old(ind_old)).and.(ind_old<size(gp_pos0_old)))
                ind_old = ind_old + 1
            enddo
            ! Common multiplication factor for interpolation: ( r0 - (r0old-) )/( (r0old+) - (r0old-) )
            ! + and - indicate above and below r0
            mult_factor = (gp_pos0(ind)-gp_pos0_old(ind_old-1))/(gp_pos0_old(ind_old)-gp_pos0_old(ind_old-1))
            
            ! Interpolation: v = v(r0old-) + (v(rold+)-v(rold-))*mult_factor
            gp_stress(:,igp,iel) = f_data_old%sim(1)%end_res%stress_end(ind_old-1, :) &
                + (f_data_old%sim(1)%end_res%stress_end(ind_old,:)-f_data_old%sim(1)%end_res%stress_end(ind_old-1,:))*mult_factor
            
            gp_sv(:,igp,iel) = f_data_old%sim(1)%end_res%statev_end(ind_old-1, :) &
                + (f_data_old%sim(1)%end_res%statev_end(ind_old,:)-f_data_old%sim(1)%end_res%statev_end(ind_old-1,:))*mult_factor
            
            ! Kinematic calculations
            eps_old0 = 0.d0
            call shapefun(gp_pos0(ind), rpos(:,iel), dNpdR, Np, ue_r, rgp, dr_dr0)
            call FBeps_calc(gp_F(:,igp,iel), Bs_mat, eps0, gp_strain(:,igp,iel), dNpdR, Np, dr_dr0, rgp, u0(2), u0(1), gp_pos0(ind), h0, eps_old0)
            
            ind = ind + 1
        enddo
        
    enddo
    
    ! If requested, write output in forms of old and new values as functions of the undeformed coordinates
    if (f_data%glob%resnr>0) then
        open(newunit=res_fid, file=trim(f_data%glob%outname)//'_sim'//int2str(simnr)//'_'//int2str(f_data%glob%resnr)//'.txt', status='replace')
        write(res_fid, "(A)") 'Result file for interpolation during remeshing'
        write(res_fid, *) ''    ! Extra line break
        
        ! Print old nodal values
        format_spec = '('//int2str(size(r0_old))//'ES15.5E3)'
        format_spec = f_data%sim(simnr)%outp%dbl_format
        format_spec = '('//int2str(size(r0_old))//trim(format_spec(2:len(format_spec)))
        write(res_fid, "(A)") 'Old nodal values (rows = node positions, displacements)'
        write(res_fid, format_spec) r0_old
        write(res_fid, format_spec) u0_old(3:)
        write(res_fid, *) ''    ! Line break
        
        ! Print new nodal values
        format_spec = f_data%sim(simnr)%outp%dbl_format
        format_spec = '('//int2str(size(r0))//trim(format_spec(2:len(format_spec)))
        write(res_fid, "(A)") 'New nodal values (rows = node positions, displacements)'
        write(res_fid, format_spec) r0
        write(res_fid, format_spec) u0(3:)
        write(res_fid, *) ''    ! Line break
        
        ! Print the old gauss point values
        format_spec = f_data%sim(simnr)%outp%dbl_format
        format_spec = '('//int2str(size(gp_pos0_old))//trim(format_spec(2:len(format_spec)))
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
        format_spec = f_data%sim(simnr)%outp%dbl_format
        format_spec = '('//int2str(size(gp_pos0))//trim(format_spec(2:len(format_spec)))
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
    
end subroutine remesh

function interpolate(x, y, xv) result(yv)
! Given the values y at x-values x, obtain the values yv at positions xv by interpolation
! x and xv must be sorted in ascending order
implicit none
    double precision    :: y(:), x(:), xv(:)
    double precision, allocatable   :: yv(:)
    integer             :: k1, num_base_points, num_out_points, first_above
    double precision    :: numtol
    
    numtol = minval(x(2:size(x))-x(1:(size(x)-1)))*1.d-1 ! Tolerance for extrapolation, max 10 % of smallest interval
    
    num_base_points = size(y)
    num_out_points = size(xv)
    
    allocate(yv(num_out_points))
    
    if (xv(1)<x(1)-numtol) then
        call write_output('Interpolation failed, x(1) must be <= xv(1)', 'status', 'sim:atp_element_removal')
        call write_output('x(1)='//dbl2str(x(1))//', xv(1)='//dbl2str(xv(1)), 'status', 'sim:atp_element_removal', loc=.false.)
        interpolate_failed = .true.
        return
    endif
    
    if (xv(num_out_points)>x(num_base_points)+numtol) then
        call write_output('Interpolation failed, x(end) must be >= xv(end)', 'status', 'sim:atp_element_removal')
        call write_output('x(end)='//dbl2str(x(num_base_points))//', xv(end)='//dbl2str(xv(num_out_points)), 'status', 'sim:atp_element_removal', loc=.false.)
        interpolate_failed = .true.
        return
    endif
    
    first_above = 2
    do k1 = 1,num_out_points
        do while ((x(first_above)<xv(k1)).and.(first_above<num_base_points))
            first_above = first_above + 1
        enddo
        yv(k1) = y(first_above-1) + (y(first_above)-y(first_above-1))*(xv(k1)-x(first_above-1))/(x(first_above)-x(first_above-1))
    enddo
    
end function

    
end module 

