! Note, in general it is better to use the interface_mod.f90 located in matmodfits src/gen folder to ensure compatability
module interface_mod
    implicit none
    
    private
    public :: sim_template  ! Simulation template
    abstract interface
        
        subroutine sim_template(error, xvar, f_data, evec)
        use iso_c_binding
            double precision, intent(out)               :: error                        ! Error from simulation
            double precision, intent(in)                :: xvar(:)                      ! Variables to be simulated 
            type(c_ptr)                                 :: f_data                       ! Data to be used in simulation
            double precision, intent(out), allocatable, optional     :: evec(:)         ! If requested, supplies vector with errors e_i
        end subroutine
        
    end interface
    
    contains 
        
end module interface_mod
    
    
subroutine opt(simfun, xvar, lb, ub, user_data, error, f_data)
!DEC$ ATTRIBUTES DLLEXPORT :: opt
use iso_c_binding
use interface_mod
implicit none
    ! Input/output variables
    procedure(sim_template)                     :: simfun           ! Simulation procedure
    double precision, intent(inout)             :: xvar(:)          ! Initial parameters, to be updated to final optimized parameters
    double precision, intent(in)                :: lb(:), ub(:)     ! Lower and upper bounds for xvar
    double precision, allocatable, intent(in)   :: user_data(:)     ! User data to be used as optimization settings
    double precision, intent(out)               :: error            ! Error at optimized parameters
    type(c_ptr)                                 :: f_data           ! Data to be used in simulation, to be passed as pointer in user routines
    ! Internal variables
    double precision                            :: xbest(size(xvar))
    double precision                            :: ebest
    double precision                            :: alpha
    integer                                     :: n, nmax
    
    
    ! Make a very simple and stupid optimization algorithm
    ! Make a number of random guesses and see which one produce the lowest error
    nmax = int(user_data(1))
    call simfun(error, xvar, f_data)
    ebest = error
    xbest = xvar
    do n=2,nmax
        call random_number(xvar)    ! Random between 0 and 1
        xvar = xvar*(ub-lb) + lb    ! Random between lb and ub
        call simfun(error, xvar, f_data)
        if (error<ebest) then
            ebest = error
            xbest = xvar
            write(*,*) 'ebest = ', ebest
            write(*,*) xbest
            write(*,*) 
        endif
    enddo
    error = ebest
    xvar = xbest

end subroutine
