! Define datatypes for saving parameters relevant for different parts of the program
! These are primarily connected to reading of parameters, but some are added during 
! the program running. 
    
module types_mod
use umat_mod
use usr_interface_mod
use iso_c_binding
    implicit none
    
    private	! Make private by default
    public strl
    ! Main categories
    public glob_typ, sim_typ, opt_typ
    ! Simulation sub-categories
    public mesh1d_typ, exp_typ, iter_typ, err_typ, outp_typ, init_typ, atp_mr_typ, usr_sim_typ
    ! Optimization sub-categories
    public start_typ, end_cond_typ, usr_opt_typ
    
    ! Other categories used internally
    public evectmp_typ, fdata_typ
    
    integer, parameter :: strl = 200

    ! == Global settings == 
    type glob_typ
        ! Run-type options
        integer, allocatable            :: run_type(:)          ! Integers describing what do run (single sim, optimization, parameter space etc.)
        
        ! Material model options
        character(len=80)               :: cmname='noname'      ! Material model name (abaqus default length)
        logical                         :: nlgeom=.false.       ! Finite strains? 
        integer                         :: nstatv=1             ! Number of state variables
        double precision, allocatable   :: ipar_init(:)        ! Initial material parameters (len = numparam)
        integer, allocatable            :: ipar_optim(:)       ! Material parameters t.b. optimized (len = numoptim)
        integer, allocatable            :: ipar_sctype(:)      ! Scaling type for parameters (lin, log, inv) (len = numoptim)
        double precision, allocatable   :: ipar_min(:)         ! Minimum values for optimized param (len = numoptim)
        double precision, allocatable   :: ipar_max(:)         ! Maximum values for optimized param (len = numoptim)
        double precision, allocatable   :: xvar_sets(:,:)       ! Sets of xvar values for different parameters, only used if ipar_sets is not set. 
        double precision, allocatable   :: ipar_sets(:,:)       ! Sets of mpar values for different parameters (only optimized parameters)
        character(len=strl)             :: umat_lib=''          ! Name of umat library to use
        procedure(umat_template),pointer,nopass:: umat_address  ! Pointer address to material routine
        
        ! User scaling of material parameters
        logical                         :: user_scale=.false.   ! Logical if a user supplied scaling subroutine is used
        character(len=strl)             :: user_scale_lib=''    ! Name/path to user scale library to use (should contain xvar_to_mpar and ipar_to_xvar)
        procedure(user_xvar_to_mpar_template),pointer,nopass:: xvar_to_mpar_addr ! Pointer address to user_xvar_to_mpar subroutine
        procedure(user_ipar_to_xvar_template),pointer,nopass:: ipar_to_xvar_addr ! Pointer address to user_ipar_to_xvar subroutine
        procedure(user_xvar_to_ipar_template),pointer,nopass:: xvar_to_ipar_addr ! Pointer address to user_xvar_to_ipar subroutine
        
        ! Other global settings
        double precision                :: num_grad_pert=1.e-8  ! Pertubation for calculating numerical gradient of objective function
	    logical				            :: error_history=.false.! If error history (per simulation call) should be saved to file <inputfile>.err
        integer                         :: opt_resnr=0          ! Set to number > 0 if results should be written during optimization. 
                                                                ! Typically used to make very custom error evalutions via external script
        
        ! Non-settable options
        character(len=strl)             :: outname              ! Base name for output files, determined from the input file name!
        integer                         :: resnr=0              ! Results should be saved to file with nr resnr if resnr>0
        
    end type glob_typ

    ! == Simulation settings ==
    ! External command
    type ext_cmd_typ
        character(len=strl)             :: script=''            ! Script which defines the external command
        character(len=strl)             :: logfile=''           ! Logfile to append stdout and stderr to if output results requested. No output if no file specified. 
    end type ext_cmd_typ
    
    ! Mesh settings
    type mesh1d_typ
        double precision, allocatable   :: node_pos(:)          ! Radial position of each node
        double precision, allocatable   :: node_pos_undef(:)    ! Initial nodal position (for continued simulation<0 node_pos may be deformed coords. node_pos_undef is updated on each import in atp)
        double precision                :: h0=1.d0              ! Height of considered gauge section
        integer                         :: ngp=1                ! Number of gauss points in each element
        logical                         :: abaqus_bbar=.false.  ! Use the abaqus b-bar method for linear elements with 2 gauss points
        integer                         :: element_order=1      ! Which element order (1st and 2nd supported, 1st default)
    end type mesh1d_typ
    
    ! Experiment settings
    type exp_typ
        double precision, allocatable   :: ctrl(:,:)            ! 
        character(len=strl)             :: exp_file             ! Experiment data file containing
        integer, allocatable            :: exp_info(:)          ! Gives the columns in the experimental data file if used, otherwise put to zero
        double precision, allocatable   :: exp_scale(:)         ! Scale factor for each column in the experiment data
        logical                         :: intpres_extstrn=.false.! Special control type for atp that gives strain control on outer surface as well as ensuring zero external pressure
    end type exp_typ
    
    ! Iteration settings
    type iter_typ
        ! General iteration settings
        integer                         :: max = 20                 ! Maximum number of iterations
        double precision                :: tol = 1.d-6              ! Iteration tolerance
        
        ! Adaptive time stepping settings
        double precision, allocatable   :: time_incr(:,:)           ! (step number, dtmain, dtmin, dt0, dtmax)*
        integer                         :: nconv_incr = 4           ! How many converged increments before increasing time step
        double precision                :: dt_incr = 1.5d0          ! How much to increase time step with if converged nconv_incr times
        
        ! Line search iteration settings
        integer                         :: ls_start = 3             ! Number of iterations before line search kicks in if error>error_old
        double precision                :: ls_alpha_factor = 0.5d0  ! How much the reduction factor is reduced every time error still is > error_old
        
    end type iter_typ 
    
    ! Error settings
    type err_typ
        integer                         :: err_norm_met=0       ! Method for scaling errors (i.e. no scaling, max-min value in cycle, max-min value in simulation)
        double precision, allocatable   :: err_scale(:,:)       ! Weight factors for different steps (step, axial, torsion, radial_inner, radial_outer)*
        double precision, allocatable   :: err_steps(:)         ! Which steps to calculate error for. error_steps=0 (default) saves for all
        double precision, allocatable   :: err_scale_ctrl(:,:)  ! Weight factors for error on controlled quantity (step data type) (only applicable for abs(ctrl)>1)
        integer                         :: error_type = 1       ! Error type (0=user subroutine, 1=square error sum, 2=cyclic error)
        
        ! Error counting variables, set by code
        integer                         :: hist_rows = -1           ! Used to know the size of the error_hist variable from previous simulations
        
        ! User subroutine error fields (error_type=0)
        character(len=strl)             :: error_lib=''             ! Path or name of user subroutine for error calculation
        procedure(user_error_template),pointer,nopass:: error_address    ! Pointer address to error calculation subroutine
        double precision, allocatable   :: user_settings(:)         ! Settings for user error calculation
        
        ! Square error fields (error_type=1)
        ! No fields required for this
        
        ! Cyclic error fields (error_type=2)
        double precision                :: nstep_cyc_err_calc=1.d0  ! What difference in the step column builds up one cycle?
        double precision                :: nstep_cyc_initial=0.d0   ! What value in the step column before the first cycle?
        double precision                :: err_type_scale=0.5       ! Weight of cycle average error versus cycle shape error
        
    end type err_typ
    
    ! Initial condition settings
    type init_typ
        integer                         :: cont_analysis=0      !=0 imply no continued analysis, >0 indicates continued from analysis number cont_analysis (and hence initial conditions for u should be used, initial conditions for statev always applied)
        double precision, allocatable   :: statev_init(:,:)     !Initial state variables (normally initialized to zeros)
        double precision                :: temp_init=0.d0       !Initial temperature, default is zero
    end type init_typ
    
    ! Result output settings
    type outp_typ
        logical                         :: result_onlymain=.true.   ! True if only results for main increments, otherwise for all increments
        logical                         :: result_inclexp=.false.   ! True if experiment data should be included in result file, false (default) otherwise
        double precision, allocatable   :: result_steps(:)          ! Which steps to write out result file for, if save_results=.true. 
                                                                    ! If result_steps=0 then all steps are saved (default)
                                                                    ! If result_steps<0 then no result file is made
        character(len=strl)             :: dbl_format=''            ! Output format for float values in result file (default set in check_input_mod.f90)
        integer                         :: log_output=1             ! Amount of output to write to log file (0: Only errors, 1: simulation summarizing messages, 2: full output)
        
        ! Additional output
        integer, allocatable            :: output_nodes(:)          ! Which nodes to output from (default is all) (Only applicable to atp simulation, ignored for mps)
        logical                         :: ur=.false.               ! Should radial displacements be output? (ignored for mps)
        integer, allocatable            :: output_elems(:)          ! Which elements to output from. Always all integration points. (Only applicable to atp simulation, ignored for mps)
        logical                         :: stress=.false.           ! Should stress be output? (ignored for mps)
        integer, allocatable            :: stress_comp(:)           ! Which stress components to output. nlgeom: (11,22,33,12,32,31,13,21,32), else (11,22,33,12,13,23)
        logical                         :: strain=.false.           ! Should strain be output? (ignored for mps)
        integer, allocatable            :: strain_comp(:)           ! Which strain components to output. nlgeom: (11,22,33,12,32,31), else (11,22,33,12,13,23)
        logical                         :: dfgrd=.false.            ! Should deformation gradient be output? (ignored for mps)
        integer, allocatable            :: dfgrd_comp(:)            ! Which deformation gradient components to output. Always: (11,22,33,12,32,31,13,21,32)
        logical                         :: statev=.false.           ! Should state variables be output?
        integer, allocatable            :: statev_comp(:)           ! Which state variable components to output.
        logical                         :: statev_norm=.false.      ! Should the norm of all state variables be written to output
        ! If the corresponding *_comp are set by user, * should be set to true!
        
    end type outp_typ
    
    ! ATP Material Removal settings
    type atp_mr_typ
        double precision                :: time_relx = 1.d0         ! Time for relaxation to zero load before material removal
        double precision                :: time_remesh = 1.d0       ! Time for remesh step (always done in one step, but time affects rate dependent materials)
        integer                         :: geom_iter_max = 5        ! Max number of iterations for finding the correct node positions
        double precision                :: node_pos_tol=1.d-6       ! Tolerance for node positions
        integer                         :: max_num_refinements=5    ! Maximum number of refinements per depth
        integer                         :: max_recursion_depth=10   ! Maximum times the recursive subroutine is called (technically not recursion depth)
    end type atp_mr_typ
    
    
        ! External user simulation subroutine
    type usr_sim_typ
        character(len=strl)             :: lib = ''             ! Name of dll/so library file for user simulation
        double precision, allocatable   :: user_data(:)         ! User simulation settings
        procedure(user_sim_template),pointer,nopass:: sim_addr  ! Pointer address to user simulation routine
    end type
    
    
    ! Values calculated during simulation setup
    type sim_setup_typ
        double precision, allocatable   :: expdata_array(:,:)   ! Data array with the contents of the exp_data file
        double precision, allocatable   :: disp_error_scale(:,:)! Error scaling for displacement
        double precision, allocatable   :: load_error_scale(:,:)! Error scaling for load
        integer, allocatable            :: stprows(:)           ! Rows in the experiment on which new steps start
        double precision, allocatable   :: steps(:)             ! Float values for each step
    end type sim_setup_typ
    
    ! End values after each simulation
    type end_res_typ
        double precision, allocatable   :: stress_end(:,:)      !Stress at end of simulation
        double precision, allocatable   :: strain_end(:,:)      !Strain at end of simulation
        double precision, allocatable   :: dfgrd_end(:,:)       !Deformation gradient at end of simulation
        double precision, allocatable   :: statev_end(:,:)      !State variables at end of simulation
        double precision, allocatable   :: u_end(:)             !Displacements for each dof at the end of the simulation
        double precision                :: disp_end(4)          !Simulated "strain" (eps_z, phi, cstr_i, cstr_o) at end of simulation
        double precision                :: load_end(4)          !Simulated loads (Fz, Torque, p_i, p_o) at end of simulation
        double precision                :: h0_true              !True initial height (=h0 usually, but if adjusted in prev cont analysis it can be different)
    end type end_res_typ
    
        
    ! Main simulation settings collection
    type sim_typ     ! Simulation settings
        ! General parameters
        integer                         :: stype                ! Simulation type
        integer                         :: n_errors=-1          ! Number of error calculated for current simulation (internal) (Used to allocate if needed)

        type(ext_cmd_typ)               :: ext_cmd              ! External simulation settings
        type(mesh1d_typ)                :: mesh1d               ! Mesh specification
        type(exp_typ)                   :: exp                  ! Experiment info
        type(iter_typ)                  :: iter
        type(err_typ)                   :: err                  ! Objective function scaling values
        type(init_typ)                  :: init                 ! Initial conditions
        type(outp_typ)                  :: outp                 ! Output settings
        type(atp_mr_typ)                :: atp_mr               ! ATP material removal settings
        
        type(sim_setup_typ)             :: sim_setup            ! Data found during the setup_simulation procedure
        type(end_res_typ)               :: end_res              ! Resulting conditions
        
        type(usr_sim_typ)               :: usr_sim              ! User simulation settings (stype=-1)

    end type sim_typ

    ! == Optimization settings ==
    
    type start_typ
        integer                         :: algorithm                ! Number of the algorithm to use
        double precision, allocatable   :: initial_step(:)          ! Initial step for algorithm
        integer                         :: num_sets=1               ! Maximum number of parameter sets (can be less if fewer are given as input)
        integer                         :: fixed_seed=0             ! Fixed seed, only used if user specifies it and thus set use_fixed_seed=.true. implictly
        logical                         :: use_fixed_seed=.false.   ! Is set true if option *fixed_seed is specified
    end type start_typ
    
    type end_cond_typ
        double precision                :: stopval=0.d0     ! Value of objective causing termination
        double precision                :: ftol_rel=0.d0    ! Relative tolerance on objective change to cause termination
        double precision                :: ftol_abs=0.d0    ! Absolute tolerance on objective change to cause termination
        double precision, allocatable   :: xtol_rel(:)      ! Relative tolerance on parameter change to cause termination
        double precision, allocatable   :: xtol_abs(:)      ! Absolute tolerance on parameter change to cause termination
        integer                         :: maxeval=0        ! (Approximate) maximum number of function evaluation
        double precision                :: maxtime=0.d0     ! (Approximate) maximum time in seconds for optimization
    end type end_cond_typ
    
    type usr_opt_typ
        character(len=strl)             :: lib = ''             ! Name of dll/so library file for user optimization
        procedure(user_opt_template),pointer,nopass:: opt_addr  ! Pointer address to user optimization routine
        double precision, allocatable   :: user_data(:)         ! User settings for optimization routine
    end type usr_opt_typ
    
    
    ! Optimization settings (specific to each optimization step)
    type opt_typ
        integer                         :: otype            ! Which optimization method to use for this stage
        type(start_typ)                 :: start            ! Set the parameters controlling the start of the optimization
        type(end_cond_typ)              :: end_cond         ! Set the parameters determining when the optimization stage should terminate
        type(usr_opt_typ)               :: usr_opt          ! User defined optimization (otype=-1)
    end type opt_typ

    type optinfo_typ    ! Information saved during optimization
        integer                         :: n_objfun_calls   ! Number of times objective function called
        double precision                :: min_objfun       ! The current lowest objective function value
        double precision                :: start_time       ! The starting time for current optimization run
    end type optinfo_typ
    
    type fdata_typ  !Data container for passing into objective function using f_data argument
        type(glob_typ)                  :: glob
        type(sim_typ), allocatable      :: sim(:)
        type(optinfo_typ)               :: optdata
    end type fdata_typ
    
    type evectmp_typ
        double precision, allocatable   :: evec_tmp(:)  ! Temporary storage for error vector
    end type evectmp_typ

    end module types_mod
   
