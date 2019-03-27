module sim_writing_util_mod
use types_mod
implicit none
integer  :: fid_res

private

public  :: setup_result_output, close_result_output
public  :: write_result

character(len=20)   :: dbl_format
    contains
    
subroutine close_result_output(save_results, simtime, error)
implicit none
    logical, intent(in)         :: save_results
    double precision, optional  :: simtime
    double precision, optional  :: error
    
    if (save_results) then
        if (present(simtime)) then
            write(fid_res, "(A, F0.3, A)") '% Simulation time:  ', simtime, 's'
        endif
        if (present(error)) then
            write(fid_res, "(A, E15.5)") '% Simulation error: ', error
        endif
        close(fid_res)
    endif
    
end subroutine

subroutine write_result(step, incr, niter, time, temp, load, load_exp, disp, disp_exp, outp, gp_stress, gp_strain, gp_dfgrd, gp_statev, ur)
    implicit none
!   character(len=20), parameter :: format_spec = '(p1,E15.5E3)'
    double precision, intent(in) :: step
    integer, intent(in)          :: incr, niter
    double precision, intent(in) :: time, temp, disp(:), disp_exp(:), load(:), load_exp(:)
    type(outp_typ), intent(in)   :: outp
    ! Additional output. For mps simulation use (e.g. for stress) reshape(stress, [size(stress), 1, 1])
    double precision, intent(in) :: gp_stress(:,:,:)
    double precision, intent(in) :: gp_strain(:,:,:)
    double precision, intent(in) :: gp_dfgrd(:,:,:)
    double precision, intent(in) :: gp_statev(:,:,:)
    double precision, optional, intent(in) :: ur(:)
    ! Internal variables
    integer                      :: k1, iel, igp
    
    
    !Write step number
    write(fid_res, "(F9.2)", advance="no") step
    
    !Write increment number
    write(fid_res, "(I9)", advance="no") incr
    
    !Write increment number
    write(fid_res, "(I9)", advance="no") niter
        
    !Write time
    write(fid_res, dbl_format, advance="no") time
    
    
    !Write channels (axial, torsional, inner, outer)
    !Load first, then "displacement". 
    do k1=1,size(load)
        ! Load
        write(fid_res, dbl_format, advance="no") load(k1)
        if (outp%result_inclexp) then
            write(fid_res, dbl_format, advance="no") load_exp(k1)
        endif
        
        ! "Displacement"
        write(fid_res, dbl_format, advance="no") disp(k1)
        if (outp%result_inclexp) then
            write(fid_res, dbl_format, advance="no") disp_exp(k1)
        endif
    enddo
    
    !Write temperature
    write(fid_res, dbl_format, advance="no") temp
    
    !Write gauss point output
    do k1=1,size(outp%output_elems)
        iel = outp%output_elems(k1)
        do igp=1,size(gp_stress,2)
            ! Stress output
            if (outp%stress) then
                call write_gp_output(igp, iel, gp_stress, outp%stress_comp)
            endif
            ! Strain output
            if (outp%strain) then
                call write_gp_output(igp, iel, gp_strain, outp%strain_comp)
            endif
            ! Deformation gradient output
            if (outp%dfgrd) then
                call write_gp_output(igp, iel, gp_dfgrd, outp%dfgrd_comp)
            endif
            ! State variable output
            if (outp%statev) then
                call write_gp_output(igp, iel, gp_statev, outp%statev_comp)
            endif
        enddo
    enddo
    
    ! Write nodal values
    if (outp%ur.and.present(ur)) then
        do k1=1,size(outp%output_nodes)
            write(fid_res, dbl_format, advance="no") ur(outp%output_nodes(k1))
        enddo
    endif
    
    !Linebreak
    write(fid_res, "(A)") ''
    
end subroutine
                
subroutine write_gp_output(igp, iel, gp_val, gp_comp)
! Given integration point and element number, write the values of gp_val specified by gp_comp to fid_res
implicit none
integer, intent(in)          :: igp, iel
double precision, intent(in) :: gp_val(:,:,:)
integer, intent(in)          :: gp_comp(:)
integer                      :: k1

    do k1=1,size(gp_comp)
        write(fid_res, dbl_format, advance="no") gp_val(gp_comp(k1), igp, iel)
    enddo

end subroutine


subroutine setup_result_output(filename, stype, nlgeom, outp, ngp)
use output_mod
implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in)          :: stype, ngp
    logical, intent(in)          :: nlgeom
    type(outp_typ), intent(in)   :: outp
    ! Internal variables
    character*4             :: stp_str, inc_str, nit_str, tim_str, temp_str
    character*6,allocatable :: load_str(:), disp_str(:), stress_str(:), statev_str(:)
    character*6             :: strain_str(6), dfgrd_str(9)
    character*40            :: add_str
    integer                 :: k1, col, outp_len, iel, igp
    character(len=100)      :: tmp_str
    character(len=10)       :: dblhead_str_format
    
    open(newunit=fid_res,file=filename,status='replace')
    
    dbl_format = outp%dbl_format
    write(tmp_str,dbl_format) 0.0
    outp_len = len(trim(tmp_str))
    dblhead_str_format = '(A'//int2str(outp_len)//')'
        
    stp_str  = "Step"
    inc_str  = "Incr"
    nit_str  = "Iter"
    tim_str  = "Time"
    temp_str = "Temp"
    
    ! Gauss point output headers
    dfgrd_str = ['F11', 'F22', 'F33', 'F12', 'F23', 'F31', 'F13', 'F21', 'F31']
    if (nlgeom) then
        allocate(stress_str(9))
        stress_str = ['P11', 'P22', 'P33', 'P12', 'P23', 'P31', 'P13', 'P21', 'P31']
        strain_str = ['E11', 'E22', 'E33', 'E12', 'E23', 'E31']
    else
        allocate(stress_str(6))
        stress_str = ['sig11', 'sig22', 'sig33', 'sig12', 'sig13', 'sig23']
        strain_str = ['eps11', 'eps22', 'eps33', 'gam12', 'gam13', 'gam23']
    endif
    
    ! State variable header
    allocate(statev_str(size(outp%statev_comp)))
    do k1=1,size(outp%statev_comp)
        statev_str(k1) = 'sv'//int2str(outp%statev_comp(k1))
    enddo
    
    if (any(stype==[1,11])) then    ! ATP type simulation
        allocate(disp_str(4), load_str(4))
        disp_str = (/"eps_z ", "phi   ", "eps_ci", "eps_co"/)
        load_str = (/"F_z   ", "Torque", "p_i   ", "p_o   "/)
    elseif (stype==2) then          ! MPS
        if (nlgeom)      allocate(disp_str, source=dfgrd_str)
        if (.not.nlgeom) allocate(disp_str, source=strain_str)
        allocate(load_str, source=stress_str)
    endif
    
        
    write(fid_res, "(A1)", advance="no") "%"    !Write comment sign to facilitate reading result into matlab if desired
    
    col = 1
    !Write step number
    write(fid_res, "(A8)", advance="no") stp_str//'('//int2str(col)//')'
    col = col + 1
    
    !Write increment number
    write(fid_res, "(A9)", advance="no") inc_str//'('//int2str(col)//')'
    col = col + 1
    
    !Write number of iterations
    write(fid_res, "(A9)", advance="no") nit_str//'('//int2str(col)//')'
    col = col + 1
    
    !Write time
    write(fid_res, dblhead_str_format, advance="no") tim_str//'('//int2str(col)//')'
    col = col + 1
    
    
    !Write channels (Axial, Torsional, Inner, Outer)
    do k1=1,size(load_str)
        ! Load
        write(fid_res, dblhead_str_format, advance="no") trim(load_str(k1))//'('//int2str(col)//')'
        col = col + 1
        if (outp%result_inclexp) then
            write(fid_res, dblhead_str_format, advance="no") trim(load_str(k1))//'_exp'//'('//int2str(col)//')'
            col = col + 1
        endif
        ! "Displacement"
        write(fid_res, dblhead_str_format, advance="no") trim(disp_str(k1))//'('//int2str(col)//')'
        col = col + 1
        if (outp%result_inclexp) then
            write(fid_res, dblhead_str_format, advance="no") trim(disp_str(k1))//'_exp'//'('//int2str(col)//')'
            col = col + 1
        endif
    enddo
    
    !Write temperature
    write(fid_res, dblhead_str_format, advance="no") temp_str//'('//int2str(col)//')'
    col = col + 1
    
    ! Write additional variables' headers (if requested)
    !Write gauss point output
    do k1=1,size(outp%output_elems)
        iel = outp%output_elems(k1)
        do igp=1,ngp
            ! Stress output
            if (outp%stress) then
                call write_gp_output_header(igp, iel, stress_str, outp%stress_comp, col, stype, dblhead_str_format)
            endif
            ! Strain output
            if (outp%strain) then
                call write_gp_output_header(igp, iel, strain_str, outp%strain_comp, col, stype, dblhead_str_format)
            endif
            ! Deformation gradient output
            if (outp%dfgrd) then
                call write_gp_output_header(igp, iel, dfgrd_str, outp%dfgrd_comp, col, stype, dblhead_str_format)
            endif
            ! State variable output
            if (outp%statev) then
                call write_gp_output_header(igp, iel, statev_str, outp%statev_comp, col, stype, dblhead_str_format)
            endif
        enddo
    enddo
    
    ! Write nodal values
    if (outp%ur) then
        do k1=1,size(outp%output_nodes)
            write(fid_res, dblhead_str_format, advance="no") 'ur'//int2str(outp%output_nodes(k1))
        enddo
    endif
    
    !Write linebreak
    write(fid_res, "(A)") ''
    
end subroutine

subroutine write_gp_output_header(igp, iel, gp_str, gp_comp, col, stype, dblhead_str_format)
use output_mod
implicit none
    integer, intent(in) :: igp, iel, gp_comp(:), stype
    character(len=*)    :: gp_str(:), dblhead_str_format
    integer, intent(inout) :: col

    integer             :: k1
    character(len=10)   :: prefix = ''

    if (stype.ne.2) prefix = 'el'//int2str(iel)//':ip'//int2str(igp)//':'
    

    do k1=1,size(gp_comp)
        write(fid_res, dblhead_str_format, advance="no") trim(prefix)//trim(gp_str(k1))//'('//int2str(col)//')'
        col = col + 1
    enddo

end subroutine


end module
