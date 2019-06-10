
module simulation_mod
    use sim_setup_mod
    use atp_mod
    use mps_mod
    use atp_element_removal_mod
    use types_mod
    use convert_mpar_mod
    use output_mod
    use usr_interface_mod
    implicit none
    
    private ! by default
    public  :: setup_simulations
    public  :: simulate, simulate_gradient
    public  :: optanalyzer
    
    integer :: mpar_fid
    integer :: eext_fid
    integer :: evec_fid
    
    contains
! Main simulation routine
subroutine simulate(error, xvar, f_data, evec)
    implicit none
    double precision, intent(out)   :: error
    double precision, intent(in)    :: xvar(:)
    type(fdata_typ), intent(inout)  :: f_data
    double precision, optional, allocatable, intent(out) :: evec(:)
    type(evectmp_typ), allocatable  :: ev_tmp(:)
    integer                         :: k1, stype, nerr1, nerr2
    double precision, allocatable   :: mpar(:)
    double precision                :: error_tmp
    logical                         :: failed
    procedure(user_sim_template),pointer:: usr_sim ! Address to user subroutine
    
    failed = .false.
    
    ! Convert xvar to mpar, where mpar includes all material parameters
    allocate(mpar(size(f_data%glob%ipar_init)))
    call xvar_to_mpar(mpar, xvar, f_data%glob)
    
    ! Allocate ev_tmp if needed
    if (present(evec)) then
        allocate(ev_tmp(size(f_data%sim)))
        nerr1 = 0
    endif
    
    error = 0.d0
    do k1=1,size(f_data%sim)
        stype = f_data%sim(k1)%stype
        if      (stype == 0) then
            ! External simulation
            if (present(evec)) then
                call simulate_external(error_tmp, mpar, f_data, k1, ev_tmp(k1)%evec_tmp)
            else
                call simulate_external(error_tmp, mpar, f_data, k1)
            endif
        elseif  (stype == -1) then
            ! User simulation
            usr_sim => f_data%sim(k1)%usr_sim%sim_addr
            
            if (present(evec)) then
                call usr_sim(error_tmp, mpar, f_data%sim(k1)%usr_sim%user_data, ev_tmp(k1)%evec_tmp)
            else
                call usr_sim(error_tmp, mpar, f_data%sim(k1)%usr_sim%user_data)
            endif
        elseif  (stype == 1) then
            ! ATP simulation
            if (present(evec)) then
                call atp_simulate(error_tmp, mpar, f_data, k1, ev_tmp(k1)%evec_tmp)
            else
                call atp_simulate(error_tmp, mpar, f_data, k1)
            endif
        elseif  (stype == 2) then
            ! MP simulation (MPS)
            if (present(evec)) then 
                call mp_simulate(error_tmp, mpar, f_data, k1,  ev_tmp(k1)%evec_tmp)
            else
                call mp_simulate(error_tmp, mpar, f_data, k1)
            endif
        elseif (stype == 11) then
            ! ATP machine ("element removal") (no error output for this stage)
            call atp_element_removal(error_tmp, mpar, f_data, k1)
            if (present(evec)) then 
                allocate(ev_tmp(k1)%evec_tmp(1))
                ev_tmp(k1)%evec_tmp = 0.d0
            endif
        else
            call write_output('Simulation '//int2str(k1)//': stype = '//int2str(stype)//' is not supported', 'error', 'sim')
        endif
        ! Add current error to total error
        error = error + error_tmp
        
        
        if (present(evec)) then
            nerr1 = nerr1 + size(ev_tmp(k1)%evec_tmp)
        endif
        
        if (error>(huge(1.d0)/10)) then
            failed = .true.
            ! No point in keep simulating if error is already infinity for current parameters
            exit
        endif
        
    enddo
    
    if (present(evec)) then
        allocate(evec(nerr1))
        nerr1 = 1
        do k1=1,size(f_data%sim)
            nerr2 = (nerr1 - 1) + size(ev_tmp(k1)%evec_tmp)
            evec(nerr1:nerr2) = ev_tmp(k1)%evec_tmp
            nerr1 = nerr2 + 1
        enddo
    endif
    
end subroutine simulate

! Simulation routine if numerical gradient should be obtained for optimization
subroutine simulate_gradient(error, gradient, xvar, f_data, evec, evec_gradient, resnr)
    implicit none
    double precision, intent(out)   :: error, gradient(:)
    double precision, intent(in)    :: xvar(:)
    type(fdata_typ), intent(inout)  :: f_data
    double precision, optional, allocatable, intent(out) :: evec(:), evec_gradient(:,:)
    integer, optional               :: resnr  ! Save results for the base simulation (without pertubation?) if >0
                                        ! This options overrides whatever is in f_data%glob%resnr for the base sim
    double precision                :: error_pert, pertubation, xvar_pert(size(xvar))
    integer                         :: k1
    double precision, allocatable   :: evec_pert(:)
    integer                         :: tmp_resnr
    
    if (present(resnr)) then
        tmp_resnr = f_data%glob%resnr
        f_data%glob%resnr = resnr
    endif
    
    ! Get error at current xvar
    if (present(evec)) then
        call simulate(error, xvar, f_data, evec)
        if (present(evec_gradient)) then
            allocate(evec_gradient(size(evec), size(xvar)))
            allocate(evec_pert(size(evec)))
        endif
    else
        call simulate(error, xvar, f_data)
    endif
    
    if (present(resnr)) then
        f_data%glob%resnr = tmp_resnr
    endif
        
    xvar_pert = xvar
    pertubation = f_data%glob%num_grad_pert
    do k1=1,size(xvar)
        xvar_pert = xvar
        xvar_pert(k1) = xvar_pert(k1) + pertubation
        if (present(evec)) then
            call simulate(error_pert, xvar_pert, f_data, evec_pert)
            if (present(evec_gradient)) then
                evec_gradient(:,k1) = (evec_pert-evec)/pertubation
            endif
        else
            call simulate(error_pert, xvar_pert, f_data)
        endif
        gradient(k1) = (error_pert-error)/pertubation   !Numerical derivative
    enddo
    
end subroutine simulate_gradient

! Subroutine for setting up simulation (reading in expdata etc)
subroutine setup_simulations(f_data)
implicit none
    type(fdata_typ), intent(inout)  :: f_data
    integer                         :: k1, stype
    
    do k1=1,size(f_data%sim)
        stype = f_data%sim(k1)%stype
        if (stype==0) then      ! External simulation
            ! No setup
        elseif (stype==1) then  ! ATP simulation
            call setup_simulation(f_data, k1) 
        elseif (stype==2) then  ! MPS 
            call setup_simulation(f_data, k1)
        elseif (stype==11) then ! ATP element removal
            ! No setup
        endif
    enddo
    
end subroutine

! Subroutine for analyzing a found optimum
subroutine optanalyzer(xopt, f_data, error, dfdx, evec, corr)
    implicit none
    double precision                :: xopt(:), error
    type(fdata_typ)					:: f_data
    double precision, allocatable   :: dfdx(:), corr(:,:), evec(:), evec_gradient(:,:), evec_norm(:)
    integer                         :: k1, k2, nvar
    integer                         :: tmp_resnr
    
    nvar = size(xopt)
    
    ! allocate output depending on sizes of input
    allocate(dfdx(nvar), corr(nvar, nvar), evec_norm(nvar))
    
    ! Ensure that results are not written for each pertubation
    tmp_resnr = f_data%glob%resnr ! Save old setting
    f_data%glob%resnr = f_data%glob%opt_resnr
    
    ! Get numerical derivatives at optimum, asking for result output at the base simulation
	call simulate_gradient(error, dfdx, xopt, f_data, evec, evec_gradient, tmp_resnr)
    
    f_data%glob%resnr = tmp_resnr ! Reset old setting (incremented if a new result has been written)
    
    ! Calculate correlation
    do k1=1,nvar
        evec_norm(k1) = sqrt(sum(evec_gradient(:,k1)*evec_gradient(:,k1)))
    enddo
    
    corr = 1.d0 ! For diagonal
    do k1 = 1,nvar
        do k2 = (k1+1),nvar
            corr(k1,k2) = sum(evec_gradient(:,k1)*evec_gradient(:,k2))/(evec_norm(k1)*evec_norm(k2))
            corr(k2,k1) = corr(k1,k2)
        enddo
    enddo
    
    
end subroutine optanalyzer

! Subroutine for calling an external simulation (stype=0)
subroutine simulate_external(error, mpar, f_data, simnr, evec)
use compiler_mod
implicit none
    double precision                :: error, mpar(:), tmp
    type(fdata_typ)                    :: f_data    
    integer                         :: simnr, status, nerr, k1
    character(len=(strl+10))        :: command  !Make slightly larger than the script field of f_data%sim(simnr)
    character(len=1)                :: eout_str !String for writing logical if all errors should be written
    double precision, allocatable, optional :: evec(:)
    character(len=strl)             :: mpar_file, efile, evecfile
    
    write(mpar_file, "(A)") 'mpar_'//int2str(simnr)//'.txt'
    write(efile, "(A)") 'eext_'//int2str(simnr)//'.txt'
    write(evecfile, "(A)") 'evec_'//int2str(simnr)//'.txt'
    
    
    if (present(evec)) then
        eout_str = '1'
    else
        eout_str = '0'
    endif
    
    ! Generate command line call with additional arguments:
    ! 1) eout_str (1 if evec is requested, 0 otherwise)
    ! 2) resnr (which result file number it currently is)
    ! 3) mpar_file (mpar_simnr.txt) gives the material parameters
    ! 4) efile (eext_simnr.txt) gives the file to which the scalar error should be written
    ! 5) evecfile (evec_simnr.txt) gives the file to which the error vector should be written (if requested)
    command = trim(f_data%sim(simnr)%ext_cmd%script) // ' ' // eout_str // ' ' // int2str(f_data%glob%resnr) // ' ' &
                                             // trim(mpar_file) // ' ' // trim(efile) // ' ' // trim(evecfile) // ' ' // '>> extsim.log 2>&1 '
    
    ! Write material parameters to default file 'mpar.txt'
    open(newunit=mpar_fid, file=trim(mpar_file), status='replace', action='write', IOSTAT=status)
    write(mpar_fid, "(ES20.10E3)") mpar
    close(mpar_fid)
    
    ! Make command line call
    status = system(command)   ! Option: call execute_command_line() or call system()
    
    if (status/=0) then
        call write_output('Error in external call "'//trim(command)//'"', 'error', 'sim', halt=.false.)
        call write_output('Output is '//int2str(status), 'error', 'sim', loc=.false.)
    endif
    
    ! Aquire error from external simulation
    open(newunit=eext_fid, file=efile, status='old', action='read', IOSTAT=status)
    if (status/=0) then
        close(eext_fid)
        call write_output('Cannot open the error file '//trim(efile)//' from the external script', 'error', 'sim', halt=.false.)
        call write_output('Please verify that this file has been created', 'error', 'sim', loc=.false.)
    endif
    
    rewind(eext_fid)    ! Ensure that file is at beginning before reading
    read(eext_fid, *, IOSTAT=status) error
    if (status/=0) then
        close(eext_fid)
        call write_output('Invalid format in the error file "'//trim(efile)//'"', 'error', 'sim', halt=.false.)
    endif
    
    close(eext_fid)
    
    ! Potentially also get error vector from external
    if (present(evec)) then
        open(newunit=evec_fid, file=evecfile, status='old', action='read', IOSTAT=status)
        if (status/=0) then
            close(evec_fid)
            call write_output('Cannot open the error file '//trim(evecfile)//' from the external script', 'error', 'sim', halt=.false.)
            call write_output('Please verify that this file has been created', 'error', 'sim', loc=.false.)
        endif
        
        rewind(evec_fid)    ! Ensure that file is at beginning before reading
        
        ! Get size of evec if first time:
        if (f_data%sim(simnr)%n_errors==-1) then
            nerr = -1 ! One will also be added after status=-1, hence start with -1
            do while(status==0)
                read(evec_fid, *, IOSTAT=status) tmp
                nerr=nerr+1
            enddo
            rewind(evec_fid)
        else
            nerr = f_data%sim(simnr)%n_errors
        endif
        allocate(evec(nerr))
        
        ! Read in error vector
        do k1=1,nerr
            read(evec_fid, *, IOSTAT=status) tmp
            if (status==0) then
                evec(k1) = tmp
            else
                call write_output('Problem reading error vector input file from external script on line '//int2str(k1), 'error', 'sim', halt=.false.)
                call write_output('Note that the same number of lines should exists every time', 'error', 'sim', loc=.false.)
            endif
        enddo
        
        ! Check that the entire file has been read:
        read(evec_fid, *, IOSTAT=status) tmp
        if (status==0) then 
            call write_output('Not all lines in error vector input file read', 'warning', 'sim', halt=.false.)
            call write_output('Please ensure that the same number of lines exists every time', 'warning', 'sim', loc=.false.)
        endif
        
        close(evec_fid)
    endif
    
        
end subroutine simulate_external

end module simulation_mod
