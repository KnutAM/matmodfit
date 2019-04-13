module sim_inp_mod
use types_mod
use output_mod
use reading_utilities_mod
use check_input_mod
use gen_util_mod
use usr_interface_mod
use load_dll_mod
use iso_c_binding
implicit none
    
    private
    public      read_simulation_settings
    
    contains
    
    subroutine read_simulation_settings(sim)
    implicit none
    
    type(sim_typ), allocatable      :: sim_tmp(:), sim(:)
    integer                         :: stype, k1
    
    k1 = 0
    
    allocate(sim_tmp(maxsim))
    
    do while(.not.end_of_category('<opt>'))
        if (adjustl(textline)=='<sim>') then
            k1 = k1 + 1
            call read_int(stype)
            sim_tmp(k1)%stype = stype
            if      (stype==0) then !External script input
                call get_sim_ext(sim_tmp(k1))
            elseif (stype==-1) then !User simulation
                call get_sim_usr(sim_tmp(k1)%usr_sim)
            elseif  (stype==1) then !ATP-simulation
                call get_sim_atp(sim_tmp(k1))
            elseif  (stype==2) then !MP-simulation (MPS)
                call get_sim_mps(sim_tmp(k1))
            elseif  (stype==11) then !ATP element removal simulation
                call get_sim_atp_er(sim_tmp(k1))
            else                    !stype not defined
                call close_input()
                call write_output('stype = '//int2str(stype)//' not defined', 'error', 'inp')
            endif    
        else
            call readline()
        endif
    enddo
    
    allocate(sim(k1))
    sim = sim_tmp(1:k1)
    
end subroutine read_simulation_settings

subroutine get_sim_ext(sim)
    implicit none
    
    type(sim_typ)   :: sim
    
    do while(status==0)
        call readline()
        if (end_of_category('<sim>', '<opt>')) then
            exit
        elseif (adjustl(textline)=='*script') then
            call read_str(sim%ext_cmd%script, strl)
        elseif (adjustl(textline)=='<<ext_cmd>>') then
            ! Do nothing, as we currently only read one category there is no need to add another subprocedure...
            ! This just allows the default structure with the category
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
        endif
    enddo
    
    ! Check input for external simulation
    call check_sim_ext(sim)
    
end subroutine get_sim_ext

subroutine get_sim_usr(usr_sim)
    implicit none
    type(usr_sim_typ), intent(inout)        :: usr_sim
    
    do while(status==0)
        call readline()
        if (end_of_category('<sim>', '<opt>')) then
            exit
        elseif (adjustl(textline)=='*lib') then
            call read_str(usr_sim%lib, strl)
            call load_user_sim(usr_sim%lib, usr_sim%sim_addr)
        elseif (adjustl(textline)=='*user_data') then
            call read_dbl_vector_allocate(usr_sim%user_data)
        elseif (adjustl(textline)=='<<usr_sim>>') then
            ! Do nothing, as we currently only read in this category there is no need to add another subprocedure...
            ! This just allows the default structure with the category
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
        endif
    enddo
    
    ! Check input for external simulation
    call check_sim_usr(usr_sim)
    
end subroutine get_sim_usr

subroutine get_sim_atp(sim)
    implicit none
    
    type(sim_typ)                   :: sim
    logical                         :: read_next_line

    read_next_line = .true.
    do while(status==0)
        if (read_next_line) then
            call readline()
        endif
        read_next_line = .false.
        if (end_of_category('<sim>', '<opt>')) then
            exit
        elseif (adjustl(textline)=='<<mesh1d>>') then
            call read_sim_mesh1d(sim%mesh1d)
        elseif (adjustl(textline)=='<<exp>>') then
            call read_sim_exp(sim%exp)
        elseif (adjustl(textline)=='<<iter>>') then
            call read_sim_iter(sim%iter)
        elseif (adjustl(textline)=='<<err>>') then
            call read_sim_err(sim%err)
        elseif (adjustl(textline)=='<<init>>') then
            call read_sim_init(sim%init)
        elseif (adjustl(textline)=='<<outp>>') then
            call read_sim_outp(sim%outp)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
            read_next_line = .true.
        endif
    enddo
    
    ! Check inputs for atp simulation
    call check_sim_atp(sim) 
    
end subroutine get_sim_atp
    
subroutine get_sim_atp_er(sim)
    implicit none
    
    type(sim_typ)                   :: sim
    logical                         :: read_next_line

    read_next_line = .true.
    do while(status==0)
        if (read_next_line) then
            call readline()
        endif
        read_next_line = .false.
        if (end_of_category('<sim>', '<opt>')) then
            exit
        elseif (adjustl(textline)=='<<mesh1d>>') then
            call read_sim_mesh1d(sim%mesh1d)
        elseif (adjustl(textline)=='<<exp>>') then
            call read_sim_exp(sim%exp)
        elseif (adjustl(textline)=='<<iter>>') then
            call read_sim_iter(sim%iter)
        elseif (adjustl(textline)=='<<err>>') then
            call read_sim_err(sim%err)
        elseif (adjustl(textline)=='<<init>>') then
            call read_sim_init(sim%init)
        elseif (adjustl(textline)=='<<outp>>') then
            call read_sim_outp(sim%outp)
        elseif (adjustl(textline)=='<<atp_er>>') then
            call read_sim_atp_er(sim%atp_er)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
            read_next_line = .true.
        endif
    enddo
    
    ! Check inputs for atp element removal simulation
    call check_sim_atp_elem_rem(sim)
    
    
end subroutine get_sim_atp_er
  
subroutine get_sim_mps(sim)
    implicit none
    
    type(sim_typ)   :: sim
    logical         :: read_next_line
    
    read_next_line = .true.
    do while(status==0)
        if (read_next_line) then
            call readline()
        endif
        read_next_line = .false.
        if (end_of_category('<sim>', '<opt>')) then
            exit
        elseif (adjustl(textline)=='<<exp>>') then
            call read_sim_exp(sim%exp)
        elseif (adjustl(textline)=='<<iter>>') then
            call read_sim_iter(sim%iter)
        elseif (adjustl(textline)=='<<err>>') then
            call read_sim_err(sim%err)
        elseif (adjustl(textline)=='<<init>>') then
            call read_sim_init(sim%init)
        elseif (adjustl(textline)=='<<outp>>') then
            call read_sim_outp(sim%outp)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
            read_next_line = .true.
        endif
    enddo
    
    ! Check inputs for mp-simulation
    call check_sim_mps(sim)
    
end subroutine get_sim_mps
  
subroutine read_sim_mesh1d(mesh1d)
implicit none
    type(mesh1d_typ) :: mesh1d
    
    do while(status==0)
        call readline()
        if (end_of_subcategory()) then
            exit
        elseif (adjustl(textline)=='*node_pos') then
            call read_dbl_vector_allocate(mesh1d%node_pos)
        elseif (adjustl(textline)=='*h0') then
            call read_dbl(mesh1d%h0)
        elseif (adjustl(textline)=='*ngp') then
            call read_int(mesh1d%ngp)
        elseif (adjustl(textline)=='*abaqus_bbar') then
            call read_logical(mesh1d%abaqus_bbar)
        elseif (adjustl(textline)=='*element_order') then
            call read_int(mesh1d%element_order)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
        endif
    enddo
end subroutine

subroutine read_sim_exp(exp)
implicit none
    type(exp_typ)   :: exp
    
    do while(status==0)
        call readline()
        if (end_of_subcategory()) then
            exit
        elseif (adjustl(textline)=='*ctrl') then
            call read_dbl_mvector(exp%ctrl)
        elseif (adjustl(textline)=='*exp_data') then
            call read_str(exp%exp_file, strl)
        elseif (adjustl(textline)=='*exp_info') then
            call read_int_vector_allocate(exp%exp_info)
        elseif (adjustl(textline)=='*exp_scale') then
            call read_dbl_vector_allocate(exp%exp_scale)
        elseif (adjustl(textline)=='*intpres_extstrn') then
            call read_logical(exp%intpres_extstrn)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
        endif
    enddo
end subroutine

subroutine read_sim_iter(iter)
implicit none
    type(iter_typ)  :: iter
    integer         :: k1, k2
    
    do while(status==0)
        call readline()
        if (end_of_subcategory()) then
            exit
        elseif (adjustl(textline)=='*time_incr') then
            call read_dbl_mvector(iter%time_incr, 5, 0.d0)
            ! Insert default values
            do k1 = 1,size(iter%time_incr, 2)
                do k2=3,5
                    if (dbl_compare(iter%time_incr(k2, k1), 0.d0)) then
                        iter%time_incr(k2, k1) = iter%time_incr(2,k1)  ! Default to dtmain
                    endif
                enddo
            enddo
        elseif (adjustl(textline)=='*iter_tol') then
            call read_dbl(iter%tol)
        elseif (adjustl(textline)=='*iter_max') then
            call read_int(iter%max)
        elseif (adjustl(textline)=='*nconv_incr') then
            call read_int(iter%nconv_incr)
        elseif (adjustl(textline)=='*dt_incr') then
            call read_dbl(iter%dt_incr)
        elseif (adjustl(textline)=='*ls_start') then
            call read_int(iter%ls_start)
        elseif (adjustl(textline)=='*ls_alpha_factor') then
            call read_dbl(iter%ls_alpha_factor)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
        endif
    enddo

end subroutine

subroutine read_sim_err(err)
implicit none
    type(err_typ)   :: err
    
    do while(status==0)
        call readline()
        if (end_of_subcategory()) then
            exit
        elseif (adjustl(textline)=='*error_type') then
            call read_int(err%error_type)
        elseif (adjustl(textline)=='*error_steps') then
            call read_dbl_vector_allocate(err%err_steps)
        elseif (adjustl(textline)=='*err_norm_met') then
            call read_int(err%err_norm_met)
        elseif (adjustl(textline)=='*err_scale') then
            call read_dbl_mvector(err%err_scale)
        elseif (adjustl(textline)=='*error_lib') then
            call read_str(err%error_lib, strl)
            call load_user_error(err%error_lib, err%error_address)
        elseif (adjustl(textline)=='*user_settings') then
            call read_dbl_vector_allocate(err%user_settings)
        elseif (adjustl(textline)=='*nstep_cyc_err_calc') then
            call read_dbl(err%nstep_cyc_err_calc)
        elseif (adjustl(textline)=='*nstep_cyc_initial') then
            call read_dbl(err%nstep_cyc_initial)
        elseif (adjustl(textline)=='*err_type_scale') then
            call read_dbl(err%err_type_scale)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
        endif
    enddo
end subroutine

subroutine read_sim_outp(outp)
implicit none
    type(outp_typ)   :: outp
    
    do while(status==0)
        call readline()
        if (end_of_subcategory()) then
            exit
        elseif (adjustl(textline)=='*result_onlymain') then
            call read_logical(outp%result_onlymain)
        elseif (adjustl(textline)=='*result_inclexp') then
            call read_logical(outp%result_inclexp)
        elseif (adjustl(textline)=='*result_steps') then
            call read_dbl_vector_allocate(outp%result_steps)
        elseif (adjustl(textline)=='*dbl_format') then
            call read_str(outp%dbl_format, strl)
        elseif (adjustl(textline)=='*output_nodes') then
            call read_int_vector_allocate(outp%output_nodes)
        elseif (adjustl(textline)=='*ur') then
            call read_logical(outp%ur)
        elseif (adjustl(textline)=='*output_elems') then
            call read_int_vector_allocate(outp%output_elems)
        elseif (adjustl(textline)=='*stress') then
            call read_logical(outp%stress)
        elseif (adjustl(textline)=='*strain') then
            call read_logical(outp%strain)
        elseif (adjustl(textline)=='*dfgrd') then
            call read_logical(outp%dfgrd)
        elseif (adjustl(textline)=='*statev') then
            call read_logical(outp%statev)
        elseif (adjustl(textline)=='*stress_comp') then
            call read_int_vector_allocate(outp%stress_comp)
        elseif (adjustl(textline)=='*strain_comp') then
            call read_int_vector_allocate(outp%strain_comp)
        elseif (adjustl(textline)=='*dfgrd_comp') then
            call read_int_vector_allocate(outp%dfgrd_comp)
        elseif (adjustl(textline)=='*statev_comp') then
            call read_int_vector_allocate(outp%statev_comp)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
        endif
    enddo
end subroutine

subroutine read_sim_init(init)
implicit none
    type(init_typ)   :: init
    
    do while(status==0)
        call readline()
        if (end_of_subcategory()) then
            exit
        elseif (adjustl(textline)=='*statev_init') then
            call read_dbl_mvector(init%statev_init)
        elseif (adjustl(textline)=='*temp_init') then
            call read_dbl(init%temp_init)
        elseif (adjustl(textline)=='*cont_analysis') then
            call read_int(init%cont_analysis)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
        endif
    enddo
end subroutine

subroutine read_sim_atp_er(atp_er)
implicit none
    type(atp_er_typ)   :: atp_er
    
    do while(status==0)
        call readline()
        if (end_of_subcategory()) then
            exit
        elseif (adjustl(textline)=='*time_relx') then
            call read_dbl(atp_er%time_relx)
        elseif (adjustl(textline)=='*time_remesh') then
            call read_dbl(atp_er%time_remesh)
        elseif (adjustl(textline)=='*geom_iter_max') then
            call read_int(atp_er%geom_iter_max)
        elseif (adjustl(textline)=='*node_pos_tol') then
            call read_dbl(atp_er%node_pos_tol)
        else
            call write_output('Unknown category "'//trim(textline)//'".', 'warning', 'inp')
        endif
    enddo
end subroutine

function same_proc (a,b)
    use, intrinsic :: iso_c_binding
    logical same_proc
    type(c_funptr), intent(in) :: a,b

    same_proc = transfer(a,0_C_INTPTR_T) == transfer(b,0_C_INTPTR_T)
end function same_proc

end module
