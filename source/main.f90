program SparseDemand
! 1) Estimate Sparse quadratic demand system
!
! Revision history
! 22MAR2013 LN  updated and cleaned
!
  use nrtype
  use GlobalModule,     only : GetInputFile,InitializeParameters, &
                               ControlOptions,HHData,ComputeNMC,  &
                               DeallocateGlobalVariables,         &
                               MasterID,nWorkers
  use OutputModule,     only : DefineFileNames
  use LikelihoodModule, only : RunMonteCarlo

#if USE_MPI==1
  use GlobalModule,     only : BroadcastParameters
  ! MPI libraries
  use mpi
#endif

  implicit none

  ! random number seeds
  integer(i4b)      :: IMC1,IMC2

  ! variables used to control command inputs
  integer(i4b), parameter         :: MaxArgLength=30  ! maximum length of command line arguments
  character(len=MaxArgLength)     :: InputFile        ! name of input file

  ! MPI variables : USE_MPI==1 version
  integer(i4b)            :: ierr  ! mpi error flag
  integer(i4b)            :: pid   ! process id number

! Initialize MPI 
#if USE_MPI==1
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,pid,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,nWorkers,ierr)
  nWorkers = nWorkers - 1
#else
  pid=0
  nWorkers=1
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
  end if
 
  ! broadcast stuff
#if USE_MPI==1
  call BroadcastParameters(pid)
#endif

  if (CONTROLOPTIONS%MPIFLAG==1) then
    ! each processor, uses subset of data
    ! DefineFileNames
    if (pid==MasterID) then
      call DefineFileNames(pid)
    end if

    IMC1 = 1
    IMC2 = HHData%NMC

  ELSEIF (CONTROLOPTIONS%MPIFLAG==2) then
    ! each processor  uses subset of MC replications

    ! DefineFileNames
    call DefineFileNames(pid)

    ! compute IMC1 and IMC2 for current pid
    ! IMC1 = first random number seed for current processor
    ! IMC2 = final random number seed
    call ComputeNMC(pid,HHData%NMC,IMC1,IMC2)
  end if

  call RunMonteCarlo(IMC1,IMC2,pid)

  call DeallocateGlobalVariables

#if USE_MPI==1
  call mpi_finalize(ierr)
#endif
end program SparseDemand

