module sim_setup_mod
use types_mod
use umat_mod
use output_mod
use gen_util_mod
implicit none
    integer, parameter              :: linelength = 300 
    
    contains 
    
subroutine setup_simulation(f_data, simnr)
implicit none
    type(fdata_typ), intent(inout)  :: f_data
    integer, intent(in)             :: simnr
    
    ! Get experiment data and define the nstep and stprows fields of sim
    call get_expdata(f_data, simnr)
    ! Obtain error scaling
    call get_error_scaling(f_data, simnr)
    
end subroutine

! Read experiment data from file
subroutine get_expdata(f_data, simnr)
! import experiment into f_data and assign nstep and stprows
    implicit none
    
    type(fdata_typ), intent(inout)  :: f_data
    integer, intent(in)             :: simnr
    
    ! internal variables
    integer                         :: expdata_fid
    integer                         :: status
    integer                         :: k1, kstep, ncols, nrows
    character(len=linelength)       :: textline
    double precision, allocatable   :: exprow(:)
    integer, allocatable            :: stprows_tmp(:)
    double precision, allocatable   :: steps(:)
    double precision, allocatable   :: time(:)
    
    
    open(newunit=expdata_fid,file=f_data%sim(simnr)%exp%exp_file,status='old',iostat=status)
    if (status/=0) then
        call write_output('Problems opening the file "'//trim(f_data%sim(simnr)%exp%exp_file)//'". Check that the file exists and that you have read rights', 'error', 'sim')
    endif
    
    rewind(expdata_fid)
    
    
    ! Get number of lines (which doesn't start with comment symbol or are completely empty) in exp_data 
    k1 = 0
    read(expdata_fid, "(A)", iostat=status) textline
    call fix_first_line(textline)
    
    do while (status==0)
        if ((textline(1:1)/='!').and.(textline/="")) then
            k1 = k1 + 1
            if (k1 == 1) then
                ncols = get_ncols(textline, f_data%sim(simnr)%exp%exp_file)
            endif
            
        endif
        read(expdata_fid, "(A)", iostat=status) textline
    enddo
    
    nrows = k1
    allocate(f_data%sim(simnr)%sim_setup%expdata_array(nrows, ncols), exprow(ncols), stprows_tmp(nrows))
    rewind(expdata_fid)
    
    if (ncols<maxval(f_data%sim(simnr)%exp%exp_info)) then
        call write_output('There are not '//int2str(maxval(f_data%sim(simnr)%exp%exp_info))//' columns in "'//trim(f_data%sim(simnr)%exp%exp_file)//'". Check the *exp_info input', 'error', 'sim')
    endif
    
    
    if (size(f_data%sim(simnr)%exp%exp_scale)/=1) then
        if (size(f_data%sim(simnr)%exp%exp_scale)/=ncols) then
            call write_output('Number of elements in exp_scale must match the number of columns in the experiment data file', 'error', 'sim')
        endif
    else
        f_data%sim(simnr)%exp%exp_scale = 1.d0
    endif
    
    allocate(steps(nrows))
    status = 0
    k1 = 0
    kstep = 1
    read(expdata_fid, "(A)", iostat=status) textline
    call fix_first_line(textline)
    do while (status==0)
        if ((textline(1:1)/='!').and.(textline/="")) then
            k1 = k1 + 1
            read(textline, *, iostat=status) exprow
            if (status==0) then
                ! Add step information
                if (f_data%sim(simnr)%exp%exp_info(1)/=0) then
                    if (kstep==1) then
                        stprows_tmp(kstep) = 1
                        kstep = kstep + 1
                        steps(kstep-1) = exprow(f_data%sim(simnr)%exp%exp_info(1))
                    elseif (exprow(f_data%sim(simnr)%exp%exp_info(1))>steps(kstep-1)) then
                        stprows_tmp(kstep) = k1
                        kstep = kstep + 1
                        steps(kstep-1) = exprow(f_data%sim(simnr)%exp%exp_info(1))
                    endif
                endif
                if (ncols==size(f_data%sim(simnr)%exp%exp_scale)) then
                    f_data%sim(simnr)%sim_setup%expdata_array(k1, :) = exprow*f_data%sim(simnr)%exp%exp_scale
                else
                    f_data%sim(simnr)%sim_setup%expdata_array(k1, :) = exprow
                endif
            else
                close(expdata_fid)
                call write_output('Bad input in experimental file '//trim(f_data%sim(simnr)%exp%exp_file)//' close to line '//int2str(k1), 'error', 'sim')
            endif
        endif
        read(expdata_fid, "(A)", iostat=status) textline
    enddo
    
    if (f_data%sim(simnr)%exp%exp_info(1)==0) then    ! If no steps/cycles specified in file, use only one for the entire file
        allocate(f_data%sim(simnr)%sim_setup%stprows(2))
        f_data%sim(simnr)%sim_setup%stprows = (/1, k1/)
        allocate(f_data%sim(simnr)%sim_setup%steps(1))
        f_data%sim(simnr)%sim_setup%steps = 1.d0
    else                        
        stprows_tmp(kstep) = k1
        if (stprows_tmp(kstep-1)==k1) then
            call write_output('The last step cannot have only one line in '//trim(f_data%sim(simnr)%exp%exp_file), 'warning', 'sim')
            call write_output('Setting last step to the same as on the previous line', 'warning', 'sim', loc=.false.)
            kstep = kstep - 1
        endif
        
        allocate(f_data%sim(simnr)%sim_setup%stprows(kstep))
        f_data%sim(simnr)%sim_setup%stprows = stprows_tmp(1:kstep)
        allocate(f_data%sim(simnr)%sim_setup%steps(kstep-1))
        f_data%sim(simnr)%sim_setup%steps = steps(1:(kstep-1))
    endif
    
    ! Check that time is increasing
    allocate(time(nrows))
    time = f_data%sim(simnr)%sim_setup%expdata_array(:, f_data%sim(simnr)%exp%exp_info(2))
    if (any(time(1:(nrows-1)) > time(2:nrows)).or.(time(1)==time(nrows))) then
        call write_output('The time in expdata is not (weakly) montonically increasing, check exp_info input. Writing out the first time points', 'error', 'sim', halt=.false.)
        k1 = 0
        do while (k1 < min(nrows, 10))
            k1= k1 + 1
            call write_output(dbl2str(time(k1)), 'error', 'sim', loc=.false., halt=.false.)
        enddo
        call write_output('Exiting', 'error', 'sim', loc=.false.)
    endif
    
    
end subroutine

function get_ncols(string, exp_data) result(ncols)
implicit none
    character(len=*)    :: string
    character(len=*)    :: exp_data
    integer             :: ncols, status
    double precision    :: tmp(100), val0
    logical             :: unique_val
    
    val0 = 142234.d0    ! Some random number
    unique_val = .false.
    do while (.not.unique_val) 
        tmp = 0.d0
        read(string, *, iostat=status) tmp
        if (all(status/=(/0, -1/))) then    ! Can allow -1 (end of file, this is what we expect)
            call write_output('Bad input in experimental file '//trim(exp_data), 'error', 'sim')
        endif
        
        if (all(.not.dbl_comp_array(tmp,val0))) then
            tmp = val0
            read(string, *, iostat=status) tmp
            ncols = count(.not.dbl_comp_array(tmp,val0))
            unique_val = .true.
        else
            val0 = val0 + 1.d0  ! Increment random number if a unique val0 hasn't been determined
        endif
    enddo
    
end function

subroutine fix_first_line(textline)
implicit none
    character(len=linelength)       ::  textline
    logical                         ::  line_ok
    
    line_ok = check_line(textline)
    do while (.not.line_ok)
        textline = textline(2:) ! Remove one character in the beginning
        line_ok = check_line(textline)
    enddo
    
end subroutine
    
function check_line(textline) result(line_ok)
implicit none
    character(len=linelength)       :: textline
    double precision                :: tmp(100)
    logical                         :: line_ok
    integer                         :: status
    
    
    read(textline, *, iostat=status) tmp
    if (any(status==(/0, -1/))) then    ! Can allow -1 (end of "file", this is what we expect)
        line_ok = .true.
    elseif (textline(1:1)=='!') then
        line_ok = .true.
    elseif (textline=='') then
        line_ok = .true.
    else
        line_ok = .false.
    endif
    
end function
    
! Get error scaling
subroutine get_error_scaling(f_data, simnr)
    implicit none
    type(fdata_typ), intent(inout)  :: f_data
    integer, intent(in)             :: simnr
    
    integer, allocatable            :: exp_info(:)
    integer, allocatable            :: stprows(:)
    integer                         :: nstep, k1, kstep, nmaxcol, nchannels
    
    nmaxcol = size(f_data%sim(simnr)%exp%exp_info)
    nchannels = (nmaxcol-3)/2
    
    nstep = size(f_data%sim(simnr)%sim_setup%stprows)-1
    
    allocate(stprows(nstep+1))
    stprows = f_data%sim(simnr)%sim_setup%stprows
    
    allocate(exp_info(nmaxcol))
    exp_info = f_data%sim(simnr)%exp%exp_info
    
    allocate(f_data%sim(simnr)%sim_setup%load_error_scale(nchannels, nstep))
    allocate(f_data%sim(simnr)%sim_setup%disp_error_scale(nchannels, nstep))
    ! Predefine for those not defined (exp_info=0)
    f_data%sim(simnr)%sim_setup%load_error_scale = 0.d0
    f_data%sim(simnr)%sim_setup%disp_error_scale = 0.d0
    
    
    
    if (f_data%sim(simnr)%err%err_norm_met==0) then       ! No scaling
        do k1=1,nchannels
            if (exp_info(2+2*k1-1)/=0) then   ! Load
                f_data%sim(simnr)%sim_setup%load_error_scale(k1,:) = 1.d0
            endif
            if (exp_info(2+2*k1)/=0) then
                f_data%sim(simnr)%sim_setup%disp_error_scale(k1,:) = 1.d0
            endif
        enddo
    elseif (f_data%sim(simnr)%err%err_norm_met==1) then   ! Scale for each step
        do kstep = 1,nstep
            do k1=1,nchannels
                if (exp_info(2+2*k1-1)/=0) then
                    f_data%sim(simnr)%sim_setup%load_error_scale(k1,kstep) = get_scale(f_data%sim(simnr)%sim_setup%expdata_array(stprows(kstep):stprows(kstep+1),exp_info(2+2*k1-1)))
                endif
                if (exp_info(2+2*k1)/=0) then
                    f_data%sim(simnr)%sim_setup%disp_error_scale(k1,kstep) = get_scale(f_data%sim(simnr)%sim_setup%expdata_array(stprows(kstep):stprows(kstep+1),exp_info(2+2*k1)))
                endif
            enddo
        enddo
        
    elseif (f_data%sim(simnr)%err%err_norm_met==2) then   ! Scale with full simulation
        do k1=1,nchannels
            if (exp_info(2+2*k1-1)/=0) then
                f_data%sim(simnr)%sim_setup%load_error_scale(k1,:) = get_scale(f_data%sim(simnr)%sim_setup%expdata_array(:,exp_info(2+2*k1-1)))
            endif
            if (exp_info(2+2*k1)/=0) then
                f_data%sim(simnr)%sim_setup%disp_error_scale(k1,:) = get_scale(f_data%sim(simnr)%sim_setup%expdata_array(:,exp_info(2+2*k1)))
            endif
        enddo
    else
        call write_output('Unknown err_norm_met='//int2str(f_data%sim(simnr)%err%err_norm_met), 'error', 'sim')
    endif
    
end subroutine

function get_scale(values) result(sfac)
implicit none
	double precision  :: sfac
	double precision  :: values(:)
	double precision  :: maxv, minv
	maxv = maxval(values)
	minv = minval(values)
	if ((maxv-minv)>1d-20) then
		sfac = 1.d0/(maxv-minv)
	else
		sfac = 0.d0
	endif
end function

end module sim_setup_mod    
