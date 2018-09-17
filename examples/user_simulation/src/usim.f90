subroutine sim(error, mpar, user_data, evec)
!DEC$ ATTRIBUTES DLLEXPORT :: sim
implicit none
    ! Main user simulation subroutine
    double precision, intent(out)               :: error            ! Error to be minimized
    double precision, intent(in)                :: mpar(:)          ! Material parameters 
    double precision, allocatable, intent(in)   :: user_data(:)     ! Array of user double data
    double precision, allocatable, optional     :: evec(:)          ! Optional array of errors for each time instance, calculate if present
    
    
    error = user_data(1)*mpar(1) + user_data(2)*mpar(2)*mpar(4) - user_data(3)*mpar(2)**2
    !error = mpar(1) + mpar(2)*mpar(4) - mpar(2)**2
    if (present(evec)) then
        if (.not.allocated(evec)) allocate(evec(3))
        evec(1) = mpar(1)
        evec(2) = mpar(2)*mpar(4)
        evec(3) = mpar(2)**2
    endif
    
end subroutine

subroutine sim_input(inpfile, user_data)
!DEC$ ATTRIBUTES DLLEXPORT :: sim_input
implicit none
    character(len=*), intent(in)                :: inpfile      ! Filename of settings file to read from
    double precision, allocatable, intent(out)  :: user_data(:) ! Pointer to user settings to be used by the user_sim subroutine
    
    allocate(user_data(3))
    user_data = [1.d0, -3.d0, -4.d0]

end subroutine
