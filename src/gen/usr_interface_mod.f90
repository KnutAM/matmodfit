module usr_interface_mod

    implicit none
    
    private
    ! Scaling routines
    public :: user_xvar_to_mpar_template
    public :: user_ipar_to_xvar_template
    public :: user_xvar_to_ipar_template
    
    ! Simulation routine
    public :: user_sim_template
    
    !Optimization routine
    public :: user_opt_template
    
    
    abstract interface
    
        ! Interface for scaling routines
        subroutine user_xvar_to_mpar_template(mpar, ipar, xvar, optim, ipar_min, ipar_max)
            ! Convert scaled parameters into real material parameters. 
            ! xvar and optim will have same length, but it is OK to change other parts of mpar as long 
            ! as it is consistent with user_ipar_to_xvar
            double precision, intent(inout) :: mpar(:)      ! The material parameters (used by umat): To be updated
            double precision, intent(in)    :: ipar(:)      ! Contains the initial input parameters.
            double precision, intent(in)    :: xvar(:)      ! Contains scaled parameters [0,1] that are optimized
            integer, intent(in)             :: optim(:)     ! Contains the position of the optimized parameters in ipar
            double precision, intent(in)    :: ipar_min(:)  ! Min values for input material parameters
            double precision, intent(in)    :: ipar_max(:)  ! Max values for input material parameters

        end subroutine
        
        subroutine user_ipar_to_xvar_template(ipar, xvar, optim, ipar_min, ipar_max)
            ! Convert real material parameters into scaled material parameters
            ! xvar and optim will have the same length, but optim is just a suggestion for which parts 
            ! of mpar that will be scaled. Scaling just needs to be consistent with xvar_to_mpar and result
            ! in scaled parameters between 0 and 1
            double precision, intent(in)    :: ipar(:)      ! Contains the current material parameters
            double precision, intent(inout) :: xvar(:)      ! Contains scaled parameters [0,1] that are optimized, to be updated
            integer, intent(in)             :: optim(:)     ! Contains the position of the optimized parameters in mpar
            double precision, intent(in)    :: ipar_min(:)  ! Min values for input material parameters
            double precision, intent(in)    :: ipar_max(:)  ! Max values for input material parameters

        end subroutine
        
        subroutine user_xvar_to_ipar_template(ipar, xvar, optim, ipar_min, ipar_max)
            ! Convert scaled parameters into input material parameters. 
            double precision, intent(inout) :: ipar(:)      ! Contains the initial input parameters, to be updated
            double precision, intent(in)    :: xvar(:)      ! Contains scaled parameters [0,1] that are optimized
            integer, intent(in)             :: optim(:)     ! Contains the position of the optimized parameters in ipar
            double precision, intent(in)    :: ipar_min(:)  ! Min values for input material parameters
            double precision, intent(in)    :: ipar_max(:)  ! Max values for input material parameters
        end subroutine
        
        ! Interface for simulation routine
        subroutine user_sim_template(error, mpar, user_data, evec)
            ! Main user simulation subroutine
            double precision, intent(out)               :: error        ! Error to be minimized
            double precision, intent(in)                :: mpar(:)      ! Material parameters 
            double precision, allocatable, intent(in)   :: user_data(:) ! User settings as read by user_sim_inp
            double precision, allocatable, optional     :: evec(:)      ! Optional array of errors for each time instance, calculate if present
        end subroutine
        
        ! Interface for optimization routine
        subroutine user_opt_template(simfun, xvar, lb, ub, user_data, error, f_data)
        use iso_c_binding
        use interface_mod
            procedure(sim_template)                     :: simfun           ! Simulation procedure
            double precision, intent(inout)             :: xvar(:)          ! Initial parameters, to be updated to final optimized parameters
            double precision, intent(in)                :: lb(:), ub(:)     ! Lower and upper bounds for xvar
            double precision, allocatable, intent(in)   :: user_data(:)     ! User data to be used as optimization settings
            double precision, intent(out)               :: error            ! Error at optimized parameters
            type(c_ptr)                                 :: f_data           ! Data to be used in simulation, to be passed as pointer in user routines
        end subroutine
        
    end interface
    
    
    contains 
        
end module usr_interface_mod
    
