program Penalized_MixedLogit
  use nrtype
  use IFPORT
  use GlobalModule, only     : K,N,J,RCDIM,iFree,           &
                               DeallocateGlobalVariables,   &
                               DeallocatePenaltyParameters, &
                               AllocatePenaltyParameters,   &
                               GetInputFile,InitializeParameters, &
                               ResultStructure,      &
                               ControlOptions
  use DataModule, only       : CreateData,RescaleData
  use LikelihoodModule, only : ComputeInitialGuess,  &
                               DefineFreeParameters, &
                               MinimizeBIC,          &
                               MaximizeLikelihood,   &
                               RunMonteCarlo
  use OutputModule, only     : SaveOutputs,ComputeStats,DefineFileNames,SaveData
#if USE_MPI==1
use GlobalModule, only : BroadcastParameters1,BroadcastParamters2,SendData
use LikelihoodModule, only : LikeFunc_slave
! MPI libraries
use mpi
#endif

  implicit none

  ! variables used to control command inputs
  integer(i4b), parameter         :: MaxArgLength=30  ! maximum length of command line arguments
  character(len=MaxArgLength)     :: InputFile        ! name of input file

  integer(i4b)           :: mode
  integer(i4b)           :: model
  real(dp), allocatable  :: x(:)
  real(dp)               :: LValue
  real(dp), allocatable  :: Grad(:),Hess(:,:)
  integer(i4b)           :: ierr
  integer(i4b)           :: IMC  ! seed for random number generator

  ! stats : structure containing results
  type(ResultStructure) :: stats

  ! MPI variables : USE_MPI=-1 version
!  integer(i4b)            :: ierr  ! mpi error flag
  integer(i4b)            :: pid   ! process id number
  integer(i4b), parameter :: MasterID=0  ! id number of master processor
  integer(i4b)            :: numprocs    ! number of processors

! Initialize MPI 
#if USE_MPI==1
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,pid,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,numprocs,ierr)
#else
  pid=0
  numprocs=1
#endif

  ! find input file name and initialize parameters
  !  1) get input file name from command line
  !  2) read in control flags and store in ControlOptions
  !  3) define problem size (N,J,K,RCDIM)
  !  4) read in parameters (b,CDiag,COffDiag)
  !  5) read in maximization options:  MaxOptions
  !  6) Define penalty parameters:     Penalty
  if (pid==MasterID) then
    call GetInputFile(InputFile)
    call InitializeParameters(InputFile)
    call DefineFileNames  ! define output file names
  end if

  if (ControlOptions%NMC>0) then
    model=2
    call RunMonteCarlo(ControlOptions%NMC,model)
  else
    ! Create xData and rescale
    IMC = 1
    call CreateData(N,J,K,RCDim,IMC)
    call RescaleData
    if (ControlOptions%SaveDataFlag==1) then
      call SaveData
    end if

    model=2    

    ! initialize penalty parameters
    call AllocatePenaltyParameters(model)

    ! Define iFree
    ! iFree = structure defining which parameters are variables and which are
    !         fixed parameters
    call DefineFreeParameters(model)

    ! x     = initial guess at maximizer of likelihood function
    allocate(x(iFree%nAll))
    call ComputeInitialGuess(model,x)

    allocate(Grad(iFree%nAll),Hess(iFree%nAll,iFree%nAll))

    if (ControlOptions%BICFlag==0) then
      ierr = 0
      stats%EstimationTime = timef()
      call MaximizeLikelihood(model,x,LValue,Grad,Hess,ierr)
      stats%EstimationTime = timef()
      call ComputeStats(x,LValue,Grad,Hess,N,Stats)
    
    elseif (ControlOptions%BICFlag==1) then
      ! Loop through values of penalty%lambda
      ! for each value, compute (X,L,BIC,nNonZero)
      ! choose option with smallest BIC
      call MinimizeBIC(model,x,LValue,Grad,Hess,Stats)
    end if
   
    if (ControlOptions%OutputFlag==1) then
      call SaveOutputs(model,x,LValue,Grad,Hess,Stats)
    end if

    deallocate(x)
    deallocate(Grad,Hess)
    call DeallocatePenaltyParameters
  end if  ! if MCFlag==1

  call DeallocateGlobalVariables
end program Penalized_MixedLogit

