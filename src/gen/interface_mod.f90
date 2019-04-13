! Interfaces used internally (e.g. when setting up the optimization user interface, requiring the simulation interface)
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
    
