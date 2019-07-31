module atp_import_mod
    use output_mod
    use types_mod
    implicit none
    private
    
    public  :: atp_import_mesh      ! Import mesh/geometry
    public  :: atp_import_init      ! Import the initial conditions
    public  :: atp_export_end       ! Export the final conditions
    
    contains
! == import routines ==
!Import mesh
subroutine atp_import_mesh(mesh1d, len_adj, u0, h0, ngp, bbar, nnod, nel, ndof_tot, Rpos, disp_conv, iter_err_norm)
    implicit none
    type(mesh1d_typ), intent(inout) :: mesh1d
    double precision, intent(in)    :: len_adj, u0(:)
    double precision                :: h0, disp_conv(4)
    integer                         :: ngp, nnod, nel, ndof_tot, n_nodes, k1, k2, k3
    logical                         :: bbar
    double precision, allocatable   :: Rpos(:,:), iter_err_norm(:)
    double precision                :: ri, ro

    h0      = mesh1d%h0*len_adj             ! Initial element height (z-dimension) (Adjusted if cont_analysis < 0)
    ngp     = mesh1d%ngp                    ! Number of gauss integration points per element
    bbar    = mesh1d%abaqus_bbar            ! Use Abaqus' bbar method or not
    nel     = size(mesh1d%node_pos_undef)-1 ! Number of elements
    nnod    = 1 + mesh1d%element_order      ! Number of nodes per element
    n_nodes = nel*(nnod-1)+1                ! Total number of nodes
    ndof_tot= 2 + n_nodes                   ! Total number of degrees of freedom
    
    ! Define the Rpos and iter_err_norm variables
    allocate(Rpos(nnod, nel), iter_err_norm(ndof_tot))
    if (mesh1d%node_pos_undef(1)<1.d-20) then
        iter_err_norm(3) = 0.d0
    else
        iter_err_norm(3) = 1.d0/(mesh1d%node_pos_undef(1)*h0*2.d0)
    endif
    k3 = 4

    do k1=1,nel
        Rpos(   1,k1) = mesh1d%node_pos_undef(k1)
        Rpos(nnod,k1) = mesh1d%node_pos_undef(k1+1)
        do k2=2,(nnod-1)
            ! Linear interpolate positions between first and last node
            Rpos(k2,k1) = Rpos(1,k1) + (k2-1)*(Rpos(nnod,k1)-Rpos(1,k1))/(nnod-1)
            iter_err_norm(k3) = 1.d0/(Rpos(k2,k1)*h0*2.d0)
            k3 = k3 + 1
        enddo
        iter_err_norm(k3) = 1.d0/(Rpos(nnod,k1)*h0*2.d0)
        k3 = k3 + 1
    enddo

    ri = Rpos(1,1)
    ro = Rpos(nnod, nel)
    
    disp_conv = (/h0, 1.d0, ri, ro/)
    if (ri<1.d-12*ro) then
        disp_conv(3) = 0.d0
    endif
    
    iter_err_norm((/1,2/)) = 1.d0/(/(ro**2-ri**2), (ro**3-ri**3)/3.d0/)

end subroutine atp_import_mesh

!Import initial conditions
subroutine atp_import_init(f_data, simnr, exp_info, nstatv, gp_s0, gp_strain0, gp_F0, gp_sv0, u0, temp, disp, load, disp_exp, load_exp, time, expdata, len_adj)
use sim_getincr_mod
! import the initial conditions
implicit none
    logical, parameter, dimension(4)::  strain_inds = (/.true., .false., .true., .true./)
    integer, parameter, dimension(4)::  disp_inds = (/2,4,6,8/)
    type(fdata_typ), intent(inout)  ::  f_data
    integer, intent(in)             ::  simnr
    integer, intent(in)             ::  exp_info(:)
    integer, intent(out)            ::  nstatv
    double precision, allocatable   ::  gp_s0(:,:,:), gp_strain0(:,:,:), gp_F0(:,:,:), gp_sv0(:,:,:), u0(:)    !intent=out
    double precision, intent(out)   ::  temp, disp(4), load(4), disp_exp(4), load_exp(4), time(2), len_adj
    double precision, intent(inout) ::  expdata(:,:)
        
    ! internal variables
    logical                         ::  adj_geom    ! Should the geometry be adjusted (i.e. if cont_analysis<0)
    double precision                ::  temp_new
    integer                         ::  cont_analysis
    integer                         ::  fpos, col, k1
    integer                         ::  ngp, nel, nnod, n_nodes, ndof_tot

    ! Get some size info (copy of code in import mesh)
    ngp     = f_data%sim(simnr)%mesh1d%ngp                 ! Number of gauss integration points per element
    nel     = size(f_data%sim(simnr)%mesh1d%node_pos)-1    ! Number of elements
    nnod    = 1 + f_data%sim(simnr)%mesh1d%element_order   ! Number of nodes per element
    n_nodes = nel*(nnod-1)+1                        ! Total number of nodes
    ndof_tot= 2 + n_nodes                           ! Total number of degrees of freedom
    
    temp    = f_data%sim(simnr)%init%temp_init   ! Initial temperature
    nstatv  = f_data%glob%nstatv
    time    = expdata(1,exp_info(2))        ! Initial time same as beginning of experiment data
    len_adj = 1.d0                          ! Default length adjustment if cont_analysis>=0
    adj_geom= .false.
        
    allocate(u0(ndof_tot), gp_s0(6, ngp, nel), gp_strain0(6, ngp, nel), gp_F0(9, ngp, nel), gp_sv0(nstatv, ngp, nel))
        
    cont_analysis= abs(f_data%sim(simnr)%init%cont_analysis)
    adj_geom = f_data%sim(simnr)%init%cont_analysis<0
        
    if (cont_analysis>0) then   ! assign values from u_end and statev_end at cont_analysis to simnr's u0 and gp_sv0
        ! Check that u_end has been defined in simulation (cont_analysis) and if so assign to u0
        if (allocated(f_data%sim(cont_analysis)%end_res%u_end)) then
            u0 = f_data%sim(cont_analysis)%end_res%u_end
        else
            call write_output('u0 from analysis '//int2str(cont_analysis)//' to be continued from is not defined', 'error', 'atp')
        endif
        
        ! Check that stress_end has been defined in simulation (cont_analysis) and if so assign to gp_sv0
        if (allocated(f_data%sim(cont_analysis)%end_res%stress_end)) then
            call get_gp_values(gp_s0, f_data%sim(cont_analysis)%end_res%stress_end)
        else
            call write_output('Stress from analysis '//int2str(cont_analysis)//' to be continued from is not defined', 'error', 'atp')
        endif
        
        ! Check that strain_end has been defined in simulation (cont_analysis) and if so assign to gp_strain0
        if (allocated(f_data%sim(cont_analysis)%end_res%strain_end)) then
            call get_gp_values(gp_strain0, f_data%sim(cont_analysis)%end_res%strain_end)
        else
            call write_output('Strain from analysis '//int2str(cont_analysis)//' to be continued from is not defined', 'error', 'atp')
        endif
        
        ! Check that dfgrd_end has been defined in simulation (cont_analysis) and if so assign to gp_F0
        if (allocated(f_data%sim(cont_analysis)%end_res%dfgrd_end)) then
            call get_gp_values(gp_F0, f_data%sim(cont_analysis)%end_res%dfgrd_end)
        else
            call write_output('Deformation gradient from analysis '//int2str(cont_analysis)//' to be continued from is not defined', 'error', 'atp')
        endif
        
        ! Check if statev_init is defined by user, if so warn about potential error
        if (allocated(f_data%sim(simnr)%init%statev_init)) then
            call write_output('statev_init in continued analysis '//int2str(simnr)//' already specified: Neglecting those values', 'warning', 'atp')
        endif
        
        ! Check that statev_end has been defined in simulation (cont_analysis) and if so assign to gp_sv0
        if (allocated(f_data%sim(cont_analysis)%end_res%statev_end)) then
            call get_gp_values(gp_sv0, f_data%sim(cont_analysis)%end_res%statev_end)
        else
            call write_output('statev from analysis '//int2str(cont_analysis)//' to be continued from is not defined', 'error', 'atp')
        endif
        
        disp = f_data%sim(cont_analysis)%end_res%disp_end
        load = f_data%sim(cont_analysis)%end_res%load_end
        
        ! Check that node_pos_undef has been defined in previous simulation and allocate in current if necessary
        if (allocated(f_data%sim(cont_analysis)%mesh1d%node_pos_undef)) then
            if (.not.allocated(f_data%sim(simnr)%mesh1d%node_pos_undef)) then
                allocate(f_data%sim(simnr)%mesh1d%node_pos_undef(size(f_data%sim(cont_analysis)%mesh1d%node_pos_undef)))    
            endif
            f_data%sim(simnr)%mesh1d%node_pos_undef = f_data%sim(cont_analysis)%mesh1d%node_pos_undef
        else
            call write_output('node_pos_undef from analysis '//int2str(cont_analysis)//' to be continued from is not defined', 'error', 'atp')
        endif
        
        
        if (adj_geom) then ! New experiment (extensometer new references) [cont_analysis<0]
            ! Adjustment factor for the gauge length (adjusted in import_mesh):
            len_adj = 1.d0/(1.d0 + disp(1))
            
            ! Adjust the initial axial displacement and rotation
            u0(1:2) = u0(1:2)*(f_data%sim(simnr)%mesh1d%h0/(1.d0+disp(1)))/f_data%sim(cont_analysis)%end_res%h0_true
            
            ! Adjust the old rotation values to match the updated gauge length:
            disp(2) = disp(2)*(f_data%sim(simnr)%mesh1d%h0/(1.d0+disp(1)))/f_data%sim(cont_analysis)%end_res%h0_true
            
            ! Adjust normal strains and rotation:
            do k1=1,4
                col = exp_info(2 + disp_inds(k1))
                if (col/=0) then 
                    if (strain_inds(k1)) then
                        expdata(:,col) = disp(k1) + (1.d0+disp(k1))*expdata(:,col)
                    else
                        expdata(:,col) = u0(2) + expdata(:,col) ! Rotation, addition sufficies
                    endif
                endif
            enddo
        endif
    else
        ! Assign initial displacements and loads 
        u0 = 0.d0
        disp = 0.d0
        load = 0.d0
        
        ! Assign initial stress and strain states
        gp_s0 = 0.d0
        gp_strain0 = 0.d0
        gp_F0 = 0.d0
        gp_F0(1:3, :, :) = 1.d0
                
        
        ! Assign initial state variables
        if (allocated(f_data%sim(simnr)%init%statev_init)) then
            call get_gp_values(gp_sv0, f_data%sim(simnr)%init%statev_init)
        else
            gp_sv0 = 0.d0
        endif
        
        ! Define undeformed mesh
        if (allocated(f_data%sim(simnr)%mesh1d%node_pos_undef)) then
            f_data%sim(simnr)%mesh1d%node_pos_undef = f_data%sim(simnr)%mesh1d%node_pos
        else
            allocate(f_data%sim(simnr)%mesh1d%node_pos_undef, source=f_data%sim(simnr)%mesh1d%node_pos)
        endif
        
    endif
        
    ! Assign initial experiment values
    fpos = 2 ! Interpolate between row 1 and 2
    call get_incr(temp, load_exp, disp_exp, temp_new, time(2), expdata, exp_info, fpos)

end subroutine

subroutine get_gp_values(gp_values, array_values)
implicit none
    double precision    :: array_values(:,:)
    double precision    :: gp_values(:,:,:)
    integer             :: nel, ngp, nva
    integer             :: iel, igp, iva
    
    nel = size(gp_values,3)
    ngp = size(gp_values,2)
    nva = size(gp_values,1) ! Number of array_value columns
    
    if (size(array_values, 1)==1) then !Only one row, all gp to have same values
        do iva=1,nva    ! For each value (column) in array_values
            gp_values(iva,:,:) = array_values(1,iva)
        enddo
    else                                !All rows, each gp has different values
        iva = 1 ! Here it means row number in array_values
        do iel = 1,nel
            do igp = 1,ngp
                gp_values(:,igp,iel) = array_values(iva,:)
                iva = iva + 1
            enddo
        enddo
    endif
    
end subroutine
        
subroutine atp_export_end(f_data, simnr, u_end, gp_s_end, gp_strain_end, gp_F_end, gp_sv_end, disp_end, load_end, h0)
! Take the final simulation results and export to f_data%sim(simnr)%[u_end and statev_end]
implicit none
    type(fdata_typ), intent(inout)  ::  f_data
    integer, intent(in)             ::  simnr
    double precision, intent(in)    ::  u_end(:)
    double precision, intent(in)    ::  gp_s_end(:,:,:) ! Stress at end
    double precision, intent(in)    ::  gp_strain_end(:,:,:)
    double precision, intent(in)    ::  gp_F_end(:,:,:)
    double precision, intent(in)    ::  gp_sv_end(:,:,:)
    double precision, intent(in)    ::  disp_end(4), load_end(4)
    double precision, intent(in)    ::  h0
    integer                         ::  nel, ngp, nstatv
    
    nel = size(gp_sv_end, 3)
    ngp = size(gp_sv_end, 2)
    nstatv=size(gp_sv_end,1)
    
    ! Allocate if needed
    if(.not.allocated(f_data%sim(simnr)%end_res%u_end)) then
        allocate(f_data%sim(simnr)%end_res%u_end(size(u_end)))
    endif
    
    if(.not.allocated(f_data%sim(simnr)%end_res%stress_end)) then
        allocate(f_data%sim(simnr)%end_res%stress_end(nel*ngp, 6))
    endif
    
    if(.not.allocated(f_data%sim(simnr)%end_res%strain_end)) then
        allocate(f_data%sim(simnr)%end_res%strain_end(nel*ngp, 6))
    endif
    
    if(.not.allocated(f_data%sim(simnr)%end_res%dfgrd_end)) then
        allocate(f_data%sim(simnr)%end_res%dfgrd_end(nel*ngp, 9))
    endif
    
    if(.not.allocated(f_data%sim(simnr)%end_res%statev_end)) then
        allocate(f_data%sim(simnr)%end_res%statev_end(nel*ngp, nstatv))
    endif
    
    
    ! Set output
    f_data%sim(simnr)%end_res%u_end = u_end
    call get_array_values(f_data%sim(simnr)%end_res%stress_end, gp_s_end)
    call get_array_values(f_data%sim(simnr)%end_res%strain_end, gp_strain_end)
    call get_array_values(f_data%sim(simnr)%end_res%dfgrd_end, gp_F_end)
    call get_array_values(f_data%sim(simnr)%end_res%statev_end, gp_sv_end)
    f_data%sim(simnr)%end_res%disp_end = disp_end
    f_data%sim(simnr)%end_res%load_end = load_end
    f_data%sim(simnr)%end_res%h0_true  = h0
    
end subroutine

subroutine get_array_values(array_values, gp_values)
! Convert gauss point values into array values
implicit none
double precision    :: array_values(:,:)
double precision    :: gp_values(:,:,:)
integer             :: igp, iel, iva

iva = 1 ! Row counter for array_values.
do iel=1,size(gp_values,3)
    do igp=1,size(gp_values,2)
        array_values(iva, :) = gp_values(:,igp,iel)
        iva = iva + 1
    enddo
enddo


end subroutine
end module
