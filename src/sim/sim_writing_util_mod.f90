module sim_writing_util_mod
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

subroutine write_result(step, incr, niter, time, temp, load, load_exp, disp, disp_exp, additional_output, result_inclexp)
    implicit none
!   character(len=20), parameter :: format_spec = '(p1,E15.5E3)'
    double precision, intent(in) :: step
    integer, intent(in)          :: incr, niter
    double precision, intent(in) :: time, temp, disp(:), disp_exp(:), load(:), load_exp(:), additional_output(:)
    logical, intent(in)          :: result_inclexp
    integer                      :: k1
    
    
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
        if (result_inclexp) then
            write(fid_res, dbl_format, advance="no") load_exp(k1)
        endif
        
        ! "Displacement"
        write(fid_res, dbl_format, advance="no") disp(k1)
        if (result_inclexp) then
            write(fid_res, dbl_format, advance="no") disp_exp(k1)
        endif
    enddo
    
    !Write temperature
    write(fid_res, dbl_format, advance="no") temp
    
    !Write additional output
    do k1=1,size(additional_output)
        write(fid_res, dbl_format, advance="no") additional_output(k1)
    enddo
    
    !Linebreak
    write(fid_res, "(A)") ''
    
end subroutine

subroutine setup_result_output(filename, result_inclexp, stype, nlgeom, outp_format)
use output_mod
implicit none
    character(len=*)        :: filename
    logical                 :: result_inclexp, nlgeom
    integer                 :: stype
    character*4             :: stp_str, inc_str, nit_str, tim_str, temp_str
    character*6,allocatable :: load_str(:), disp_str(:)
    character*40            :: add_str
    integer                 :: k1, col, outp_len
    character(len=*)        :: outp_format
    character(len=100)      :: tmp_str
    character(len=10)       :: dblhead_str_format
    
    open(newunit=fid_res,file=filename,status='replace')
    
    dbl_format = outp_format
    write(tmp_str,dbl_format) 0.0
    outp_len = len(trim(tmp_str))
    dblhead_str_format = '(A'//int2str(outp_len)//')'
        
    stp_str  = "Step"
    inc_str  = "Incr"
    nit_str  = "Iter"
    tim_str  = "Time"
    temp_str = "Temp"
    
    if (any(stype==[1,11])) then    ! ATP type simulation
        allocate(disp_str(4), load_str(4))
        disp_str = (/"eps_z ", "phi   ", "eps_ci", "eps_co"/)
        load_str = (/"F_z   ", "Torque", "p_i   ", "p_o   "/)
    elseif (stype==2) then          ! MPS
        if (nlgeom) then
            allocate(disp_str(9), load_str(9))
            disp_str = ['F11', 'F22', 'F33', 'F12', 'F23', 'F31', 'F13', 'F21', 'F31']
            load_str = ['P11', 'P22', 'P33', 'P12', 'P23', 'P31', 'P13', 'P21', 'P31']
        else
            allocate(disp_str(6), load_str(6))
            disp_str = ['eps11', 'eps22', 'eps33', 'gam12', 'gam13', 'gam23']
            load_str = ['sig11', 'sig22', 'sig33', 'sig12', 'sig13', 'sig23']
        endif
    endif
    
        
    add_str  = "User defined output variables ..."
    
    
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
        if (result_inclexp) then
            write(fid_res, dblhead_str_format, advance="no") trim(load_str(k1))//'_exp'//'('//int2str(col)//')'
            col = col + 1
        endif
        ! "Displacement"
        write(fid_res, dblhead_str_format, advance="no") trim(disp_str(k1))//'('//int2str(col)//')'
        col = col + 1
        if (result_inclexp) then
            write(fid_res, dblhead_str_format, advance="no") trim(disp_str(k1))//'_exp'//'('//int2str(col)//')'
            col = col + 1
        endif
    enddo
    
    !Write temperature
    write(fid_res, dblhead_str_format, advance="no") temp_str//'('//int2str(col)//')'
    col = col + 1
    
    !Write "User defined output variables"
    write(fid_res, '(A'//int2str(3+len(trim(add_str)//'('//int2str(col)//')'))//')', advance="no") trim(add_str)//'('//int2str(col)//')'
    col = col + 1
    
    !Write linebreak
    write(fid_res, "(A)") ''
    
end subroutine


end module
