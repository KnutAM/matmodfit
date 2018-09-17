module run_types_mod
use output_mod
use types_mod
use convert_mpar_mod 
use simulation_mod

    implicit none
    
    private
    public  :: run_simulation       ! run_type=1: Simple simulation at initial material parameters
                                    ! Required input: f_data, mpar
    public  :: run_optimization     ! run_type=2: Standard optimization run
                                    ! Required input: f_data, opt, mpar
    public  :: run_optanalyzer      ! run_type=3: Analyze optimum point (correlation matrix)
                                    ! Required input: f_data, mpar
    
                                    
    contains
    
subroutine run_simulation(f_data, xvars)
    implicit none
    type(fdata_typ)                 :: f_data
    double precision                :: xvars(:,:)
    double precision, allocatable   :: ipar(:)  ! Variables
    double precision                :: error
    integer                         :: k1
    
    ! Allocate a actual material parameter vector to use for result file output
    allocate(ipar, source=f_data%glob%ipar_init)
    
    ! Simulation
    do k1=1,size(xvars,2)
        if (size(xvars,2)==1) then
            call write_output('Single simulation started', 'status')
        else
            call write_output('Single simulation of set '//int2str(k1)//'/'//int2str(size(xvars,2))//' started', 'status')
        endif
        
        ! Increment the simulation counter and write to result file
        f_data%glob%resnr = f_data%glob%resnr + 1
        call simulate(error, xvars(:,k1), f_data)
        
        ! Write simulation results to result file
        call xvar_to_ipar(ipar, xvars(:,k1), f_data%glob)
        call write_sim_results(error, ipar)
        
        call write_output('Single simulation finished', 'status')
    enddo
    
    
end subroutine

subroutine run_optimization(f_data, opt, xvars)
use optimization_mod 
    implicit none
    type(fdata_typ)                 :: f_data
    type(opt_typ)                   :: opt(:)
    double precision, allocatable   :: xvars(:,:)   ! Input parameters
    double precision, allocatable   :: ipars(:,:)   ! Input parameters
    double precision, allocatable   :: errors(:)
    integer                         :: resnr_tmp
    integer                         :: k1

    call write_output('Optimization started', 'status')
    
    ! Optimization
    resnr_tmp = f_data%glob%resnr
    f_data%glob%resnr = 0
    call optimize(errors, xvars, f_data, opt)
    f_data%glob%resnr = resnr_tmp
    
    ! Write results of optimization to result file
    allocate(ipars(size(f_data%glob%ipar_init), size(xvars,2)))
    do k1=1,size(xvars,2)
        call xvar_to_ipar(ipars(:,k1), xvars(:,k1), f_data%glob)
    enddo
    
    call write_opt_results(errors, ipars)
    
    call write_output('Optimization finished', 'status')
    
end subroutine

subroutine run_optanalyzer(f_data, xvars)
    implicit none
    type(fdata_typ)                 :: f_data
    double precision                :: xvars(:,:)
    double precision, allocatable   :: evec(:), dfdx(:), corr(:,:)
    double precision                :: error
    integer                         :: k1
    
    call write_output('Optimum analysis started', 'status')
    
    ! Optimum analyzer
    do k1=1,size(xvars,2)
        f_data%glob%resnr = f_data%glob%resnr + 1   ! Increment result number to ensure result datafiles are created.
        call optanalyzer(xvars(:,k1), f_data, error, dfdx, evec, corr)
        ! Write results to result file
        call write_opt_ana_results(error, corr)
    enddo
    
    call write_output('Optimum analysis finished', 'status')
    
    
end subroutine
    
end module
    
