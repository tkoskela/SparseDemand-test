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
                               MasterID,nWorkers,                 &
                               CreateQuadRule
  use OutputModule,     only : DefineFileNames,CopyInputFilename
  use LikelihoodModule, only : RunMonteCarlo,AnalyseResults
  use IFPORT,           only : system

#if USE_MPI==1
  use GlobalModule,     only : BroadcastParameters
  ! MPI libraries
  use mpi
#endif

  implicit none

  ! random number seeds
  integer(i4b)      :: IMC1,IMC2

  ! variables used to control command inputs
  integer(i4b), parameter         :: MaxArgLength=200  ! maximum length of command line arguments
  character(len=MaxArgLength)     :: InputFile         ! name of input file
  character(len=MaxArgLength)     :: NewInputFile      ! name of copy of input file
  character(len=MaxArgLength)     :: InputFileString   ! string used to create copy
  integer(i4b)                    :: eflag             ! error flag for copying

  ! MPI variables : USE_MPI==1 version
  integer(i4b)            :: ierr  ! mpi error flag
  integer(i4b)            :: pid   ! process id number
  integer(i4b)            :: nprocs ! number of processes

  ! Data and time
  integer(i4b)            :: DateTime(8)

! Initialize MPI 
#if USE_MPI==1
  call mpi_init(ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,pid,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr)
  nWorkers = nprocs - 1
#else
  pid=0
  nWorkers=1
  nprocs = 1
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
    ! Copy input file to output directory 
    NewInputFile    = CopyInputFilename(InputFile)
    InputFileString = "cp " // trim(InputFile) // " " // trim(NewInputFile)
    eflag = system(InputFileString)
    call date_and_time(values = DateTime)
    print 71, "Parameter initialization complete. day = ",DateTime(3), &
              ' time = ',DateTime(5),':',DateTime(6)
71 format(a41,i2,a8,i2.2,a1,i2.2)
  end if
 
  ! broadcast stuff
#if USE_MPI==1
  call date_and_time(values = DateTime)
  print 101,'Process ',pid,': start BroadcastParameters: ', &
            'day = ',DateTime(3),' time = ',DateTime(5),':',DateTime(6)
101 format(a8,i4,a30,a6,i2,a8,i2.2,a1,i2.2)
  if (nprocs>1) then
    call BroadcastParameters(pid)
  end if
  print 101,'Process ',pid,': finish BroadcastParameters: ', &
            'day = ',DateTime(3),' time = ',DateTime(5),':',DateTime(6)
#endif
  ! create integration rule
  call CreateQuadRule(pid,nprocs)

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
    if (nWorkers==0) then
      IMC1 = 1
      IMC2 = HHData%NMC
    else 
      call ComputeNMC(pid,HHData%NMC,IMC1,IMC2)
    end if
  end if

  if (ControlOptions%TestLikeFlag<=5) then 
    ! Estimate model
    if (nWorkers==0) then
      print *,'Error. RunMontecarlo only works when nWorkers>0'
      stop
    end if
    call RunMonteCarlo(IMC1,IMC2,pid)
  elseif (ControlOptions%TestLikeFlag==6) then
    ! analyse results
    IMC1 = 1
    if (pid==MasterID) then
      call AnalyseResults(IMC1)
      print *,"AnalyseResults complete."
    end if
  end if
  call DeallocateGlobalVariables
  print *,"DeallocateGlobalVariables complete."
#if USE_MPI==1
  call mpi_finalize(ierr)
  print *,"mpi_finalize complete."
#endif
end program SparseDemand

