module GlobalModule
  use nrtype
  use PropertyList   ! tools for loading parameters from file
  implicit none

  type(Property)    :: PropList   ! list of keys in input file

  type ResultStructure
    real(dp)     :: BIC
    real(dp)     :: MinEigHess
    real(dp)     :: MaxGrad
    integer(i4b) :: nNonZero
    integer(i4b) :: N
    real(sp)     :: EstimationTime
    real(dp)     :: lambda
  end type

  type FlagStructure
    integer(i4b) :: TestLikeFlag
    integer(i4b) :: OutputFlag
    integer(i4b) :: TestIntegrationFlag
    integer(i4b) :: ScalingFlag         ! 1 to rescale data
    integer(i4b) :: SaveDataFlag
    integer(i4b) :: BICFlag
    integer(i4b) :: NMC
  end type

  type QuadRuleType
    integer(i4b)              :: flag
    integer(i4b), allocatable :: n(:)
    integer(i4b)              :: nall
    real(dp),     allocatable :: nodes(:,:)
    real(dp),     allocatable :: weights(:)
  end type

  type DataSelectorType
    integer(i4b), allocatable :: b(:)
    integer(i4b), allocatable :: e(:)
  end type

  type SelectFreeType
    integer(i4b), allocatable :: b(:)
    integer(i4b), allocatable :: bX(:)
    integer(i4b)              :: nb 
    integer(i4b), allocatable :: CDiag(:)
    integer(i4b), allocatable :: CDiagX(:)
    integer(i4b)              :: nCDiag
    integer(i4b), allocatable :: COffDiag(:)
    integer(i4b), allocatable :: COffDiagX(:)
    integer(i4b)              :: nCOffDiag
    integer(i4b)              :: nAll
  end type

  type PenaltyStructureType
    integer(i4b)              :: method 
    integer(i4b)              :: nx   ! size of structural parameters
    integer(i4b)              :: nx1  ! size of unpenalized structural parameters
    integer(i4b)              :: nx2  ! size of penalized structural parameters
    integer(i4b)              :: nxP  ! size of all parameters
    integer(i4b), allocatable :: x_index1(:)
    integer(i4b), allocatable :: x_index2(:)
    integer(i4b), allocatable :: xP_index1(:)
    integer(i4b), allocatable :: xP_index2(:)
    integer(i4b), allocatable :: xP_index_xPlus(:)
    integer(i4b), allocatable :: xP_index_xMinus(:)
    integer(i4b), allocatable :: xPlus_index(:)
    integer(i4b), allocatable :: xMinus_index(:)
    real(dp)                  :: lambda
    integer(i4b)              :: nLambda
    real(dp),     allocatable :: VectorLambda(:)
    real(dp)                  :: MinLambda,MaxLambda
  end type

  type MaxStructure
    character(len=40) :: MaxIter(2)
    ! MaxAlgorithm = 1   E04WDF  : Dense problem
    ! MaxAlgorithm = 2   E04VHF  : Sparse problem
    integer(i4b)      :: Algorithm(2)
  end type

  type(FlagStructure)         :: ControlOptions
  type(MaxStructure)          :: MaxOptions

  type(PenaltyStructureType)  :: Penalty
  type(DataSelectorType)      :: iXData
  type(QuadRuleType)          :: IntRule
  type(SelectFreeType)        :: iFree

  ! Size of problem
  integer(i4b)                :: K,N,J,RCDIM

  ! Parameters
  real(dp), allocatable       :: b(:),CDiag(:),COffDiag(:)

  ! Data
  real(dp), allocatable       :: xData(:,:,:)

  ! Data labels
  character(len=200), allocatable :: XLabels(:)
  character(len=200), allocatable :: CDiagLabels(:)
  character(len=200), allocatable :: COffDiagLabels(:)

  ! Scaling parameters
  real(dp), allocatable       :: ScalingMatrix(:,:)  ! NewX = X*ScalingMatrix 
  character(len=200)          :: ScalingFile

  ! Output directory
  character(len=200)          :: OutDir
contains
! 1) subroutine AllocateGlobalVariables
! 2) subroutine AllocatePenaltyParameters
! 3) subroutine DeallocateGlobalVariables
! 4) subroutine DeallocatePenaltyParameters
! 5) subroutine GetInputFile(InputFile)
! 6) subroutine InitializeParameters(InputFile)
! 7) subroutine ReadParameters
! 8)
!------------------------------------------------------------------------------
  subroutine AllocateGlobalVariables
    implicit none
    allocate(b(K),CDiag(RCDIM),COffDiag(RCDIM*(RCDIM-1)/2))
    allocate(xData(K,J,N))
    allocate(XLabels(K),CDiagLabels(K),COffDiagLabels(RCDIM*(RCDIM-1)/2))

    ! matrix to rescale data
    ! NewX = X*ScalingMatrix
    allocate(ScalingMatrix(K,K))

    ! Allocate memory for iXData
    ! iXData%b = column indexes for columns in xData corresponding to b
    ! iXData%e = column indexes for columns in xData corresponding to e
    allocate(iXData%b(K),iXData%e(RCDIM))
  end subroutine AllocateGlobalVariables

!------------------------------------------------------------------------------
  subroutine AllocatePenaltyParameters(model)
    implicit none

    integer(i4b), intent(in) :: model
    integer(i4b) :: method,nx,nx1,nx2,nxP
    character(len=PL_VAL_LEN)       :: cTemp      ! temporary character string
    integer(i4b) :: ErrFlag
    real(dp)     :: lambda
 
    ! method = 1   max L(x1,x2) + P(x2)  subject to   x2 = xPLus-xMinus
    !                                                 0 <= xPlus 
    !                                                 0 <= xMinus
    ErrFlag = GetVal(PropList,'method',cTemp)
    read(cTemp,'(i2)') method
   
    if (model==1) then
      nx  = K
    elseif (model==2) then
      nx  = K + RCDIM
    end if
    
    ErrFlag = GetVal(PropList,'nx1',cTemp)
    read(cTemp,'(i4)') nx1

    nx2 = nx-nx1
    nxP = nx + 2*nx2

    ErrFlag = GetVal(PropList,'lambda',cTemp)
    read(cTemp,'(g11.4)') lambda

    ! number of values of lambda to try
    ErrFlag = GetVal(PropList,'nLambda',cTemp)
    read(cTemp,'(i4)') Penalty%nLambda
  
    ! minimum value of lambda
    ErrFlag = GetVal(PropList,'MinLambda',cTemp)
    read(cTemp,'(g11.4)') Penalty%MinLambda

    ! max value of lambda
    ErrFlag = GetVal(PropList,'MaxLambda',cTemp)
    read(cTemp,'(g11.4)') Penalty%MaxLambda

    Penalty%method    = method 
    Penalty%nx        = nx
    Penalty%nx1       = nx1
    Penalty%nx2       = nx2
    Penalty%nxP       = nxP
    Penalty%lambda    = lambda
    
    ! x_index1 = elements that are NOT penalized
    allocate(Penalty%x_index1(nx1))
    allocate(Penalty%xP_index1(nx1))

    Penalty%x_index1 = (/1:nx1/)
    Penalty%xP_index1 = Penalty%x_index1
 
    ! x_index2 = elements that are penalized
    allocate(Penalty%x_index2(nx2),Penalty%xP_index2(nx2))
    Penalty%x_index2  = nx1+(/1:nx2/)
    Penalty%xP_index2 = Penalty%x_index2
    
    allocate(Penalty%xP_index_xPlus(nx2))
    Penalty%xP_index_xPlus = nx + (/1:nx2/)
  
    allocate(Penalty%xP_index_xMinus(nx2))
    Penalty%xP_index_xMinus = nx+nx2+(/1:nx2/)

    allocate(Penalty%xPlus_index(nx2))
    allocate(Penalty%xMinus_index(nx2))
    Penalty%xPlus_index = (/1:nx2/)
    Penalty%xMinus_index = nx2+(/1:nx2/)

    allocate(Penalty%VectorLambda(Penalty%nLambda))
    Penalty%VectorLambda = Penalty%MinLambda + &
                           (Penalty%MaxLambda-Penalty%MinLambda)*dble((/0:Penalty%nLambda-1/))/dble(Penalty%nLambda-1)
  end subroutine AllocatePenaltyParameters

!------------------------------------------------------------------------------
  subroutine DeallocateGlobalVariables
    implicit none

    deallocate(ScalingMatrix)

    deallocate(iXData%b,iXData%e)
    if (allocated(iFree%b)) then
      deallocate(iFree%b,iFree%bX)
      deallocate(iFree%CDiag,iFree%CDiagX)
      if (allocated(iFree%COffDiag)) then
        deallocate(iFree%COffDiag,iFree%COffDiagX)
      end if
    end if
    deallocate(b,CDiag,COffDiag)
    deallocate(xData)
    deallocate(XLabels,CDiagLabels,COffDiagLabels)
    deallocate(IntRule%n,IntRule%nodes,IntRule%weights)
  end subroutine DeallocateGlobalVariables
!------------------------------------------------------------------------------
  subroutine DeallocatePenaltyParameters
    implicit none
    deallocate(Penalty%x_index1,Penalty%xP_index1)
    deallocate(Penalty%x_index2,Penalty%xP_index2)
    deallocate(Penalty%xP_index_xPlus)
    deallocate(Penalty%xP_index_xMinus)
    deallocate(Penalty%xPlus_index)
    deallocate(Penalty%xMinus_index)
    deallocate(Penalty%VectorLambda)
  end subroutine DeallocatePenaltyParameters

!------------------------------------------------------------------------------
subroutine GetInputFile(InputFile)
    implicit none
 
    character(LEN=*), intent(out)    :: InputFile              ! name of parameter input file
    integer(i4b)                     :: nArgs                  ! number of cmd. line args.
    integer, parameter               :: MaxArgLength=128       ! max. cmd. line arg. length
    character(LEN=MaxArgLength)      :: ProgramName            ! name of this executable

    ! read in parameter file name from command line argument
    nArgs = IARGC()       ! get number of arguments
    if (nArgs==1) then
      CALL GETARG( 1, InputFile)
    elseif (nArgs ==0) then
      InputFile='inputs/A1.prop'
      print *, 'Assuming default parameter file: <inputs/A1.prop>'
      print *, 'For a different parameter file, include filename as a command line option.'
      CALL GETARG( 0, ProgramName )
      PRINT *, 'For example:  ', TRIM( ProgramName ), ' <InputPropertyFile>'
      print *, ' Replace <InputPropertyFile> with the name of the input file.'
    ENDIF
  end subroutine GetInputFile

!------------------------------------------------------------------------------
subroutine InitializeParameters(InputFile)
    use IFPORT       ! intel fortran portability library
    implicit none
    character(len=*), intent(inout) :: InputFile  ! Name of input file
    integer(i4b)                    :: ErrFlag    ! error flag when reading inputs
    character(len=PL_VAL_LEN)       :: cTemp      ! temporary character string
    character(len=PL_VAL_LEN), allocatable :: cTemp1(:) ! temporary character vector
    real(dp)                        :: dTemp      ! temporary double precision real
    integer(i4b), allocatable       :: iTemp(:)   ! temporary integer array
    integer(i4b)                    :: nRead      ! number of values read
    logical                         :: DirExists  ! true if dir exists
    integer(i4b)                    :: i1

    ! read parameters from property file
    ErrFlag = LoadPropertyFile(PropList,InputFile)
    if (ErrFlag /= PL_OK) then
      print *,'Error. Could not load ',trim(InputFile)
      stop
    end if

    ! name of output directory
    ErrFlag = GetVal(PropList,'OutDir',OutDir)
    if (ErrFlag /= PL_OK) then
      print *, 'No output directory specified in input file.'
      stop
    end if

    ! check whether outputdirectory exists.
    ! If not, create.
    inquire(DIRECTORY=OutDir,exist=DirExists)
    if (.not. DirExists) then
      print *, 'Warning: Output directory does not exists.'
      print *, 'Creating ',trim(OutDir),'.'
      DirExists = makedirqq(OutDir)
      if (.not. DirExists) then
        print *, 'Error: Failed to create output directory ', trim(OutDir), '.'
        stop
      end if
    end if

   ! outputFlag: 0 do not save output
   !             1 save output
   ErrFlag = GetVal(PropList,'OutputFlag',cTemp)
   read(cTemp,'(i2)') ControlOptions%OutputFlag

   ! TestLikeFlag:  1 test likelihood
   !                0 maximise likelihood
   ErrFlag = GetVal(PropList,'TestLikeFlag',cTemp)
   read(cTemp,'(i2)') ControlOptions%TestLikeFlag

   ! TestIntegrationFlag : 0 do not test
   !                       1 test accuracy of integration
   ErrFlag = GetVal(PropList,'TestIntegrationFlag',cTemp)
   read(cTemp,'(i2)') ControlOptions%TestIntegrationFlag

   ! ScalingFlag :  0 do not rescale data
   !                1 rescale data
   ErrFlag = GetVal(PropList,'ScalingFlag',cTemp)
   read(cTemp,'(i2)') ControlOptions%ScalingFlag

  ! SaveDataFlag : 0 do not save data
  !                1 save data to file
  ErrFlag = GetVal(PropList,'SaveDataFlag',cTemp)
  read(cTemp,'(i2)') ControlOptions%SaveDataFlag

  ! BICFlag : 0 do not minimize BIC
  !           1 minimize BIC
  ErrFlag = GetVal(PropList,'BICFlag',cTemp)
  read(cTemp,'(i2)') ControlOptions%BICFlag

  ! NMC : 0 no MC repetitions
  !     > 0 number of MC repetitions
  ErrFlag = GetVal(PropList,'NMC',cTemp)
  read(cTemp,'(i4)') ControlOptions%NMC

  ! Define problem size:  (N,J,K,RCDIM)
  ! N     = sample size
  ! J     = number of options in choice set
  ! K     = number of regressors
  ! RCDIM = dimension of random coefficients
  ErrFlag = GetVal(PropList,'N',cTemp)
  read(cTemp,'(i5)') N
  ErrFlag = GetVal(PropList,'J',cTemp)
  read(cTemp,'(i5)') J
  ErrFlag = GetVal(PropList,'K',cTemp)
  read(cTemp,'(i5)') K
  ErrFlag = GetVal(PropList,'RCDIM',cTemp)
  read(cTemp,'(i5)') RCDIM

  ! allocate memory for (b,CDiag,COffDiag,xData,iXData)
  call AllocateGlobalVariables
  
  ! read in (b,CDiag,COffDiag)
  call ReadParameters

  ! set maximization options
  ErrFlag = GetVal(PropList,'MaxIter',cTemp)
  allocate(cTemp1(2))
  read(cTemp,'(2A6)') cTemp1
  do i1=1,2
    MaxOptions%MaxIter(i1) = 'Major Iterations Limit = ' // trim(cTemp1(i1))
  end do
  deallocate(cTemp1)

  ErrFlag = GetVal(PropList,'MaxAlgorithm',cTemp)
  read(cTemp,'(2i2)') MaxOptions%Algorithm

  ! define integration rule
  call DefineIntegrationRule

  end subroutine InitializeParameters

!------------------------------------------------------------------------------
! read in initial values for (b,CDiag,COffDiag)
  subroutine ReadParameters
    implicit none
    integer(i4b)       :: bUnit,CDiagUnit,COffDiagUnit
    character(len=200) :: bFile,CDiagFile,COffDiagFile
    integer(i4b)       :: i1

    ! define iXData
    iXData%b = (/1:K/)
    iXData%e = (/1:RCDIM/)

    ! b: mean impact of each regressor
    bUnit = 310
    bFile = 'inputs/b.raw'
    open(UNIT = bUnit, &
         FILE = bFile, &
         ACTION = 'read')
    b = 0.0d0
    do i1=1,K
      read(bUnit,'(g11.4)') b(i1)
    end do
    close(bUnit)

    ! read in diagonal elements of cholesky decomposition of sigma
    CDiagUnit = 321
    CDiagFile = 'inputs/CDiag.raw'
    open(UNIT = CDiagUnit, &
         FILE = CDiagFile, &
         ACTION = 'read')
    CDiag = 0.0d0
    do i1=1,RCDIM
      read(CDiagUnit,'(g11.4)') CDiag(i1)
    end do
    close(CDiagUnit)
 
    ! read in off-diagonal elements of sigma
    COffDiagUnit = 332
    COffDiagFile = 'inputs/COffDiag.raw'
    open(UNIT = COffDiagUnit, &
         FILE = COffDiagFile, &
         ACTION = 'read')
    COffDiag = 0.0d0
    do i1=1,RCDIM*(RCDIM-1)/2
      read(COffDiagUnit,'(g11.4)') COffDiag(i1)
    end do
    close(COffDiagUnit)
  end subroutine ReadParameters

subroutine DefineIntegrationRule
  implicit none
  integer(i4b) :: ErrFlag
  character(len=PL_VAL_LEN)       :: cTemp      ! temporary character string

  ! integration rule
  ! 1 = Gauss hermite
  ! 2 = Sparse Hermite
  ! 3 = Pseudo Monte Carlo
  ! 4 = Quasi Monte Carlo
  ! 5 = Taylor expand around zero
  ! 6 = Laplace exansions
  ErrFlag = GetVal(PropList,'IntegrationFlag',cTemp)
  read(cTemp,'(i2)') IntRule%flag

  allocate(IntRule%n(RCDIM))
  IntRule%n = 2
  ErrFlag = GetVal(PropList,'nQuad',cTemp)
  ! Currently, prop file contains at most 16 values for nQuad
  ! if RCDIM>16, nQuad(i1) = 2 for i1>16
  if (RCDIM<=16) then
    read(cTemp,'(<RCDIM>I4)') IntRule%n
  else
    read(cTemp,'(16I4)') IntRule%n(1:16)
  end if

  if (IntRule%flag==1) then
    IntRule%nAll = product(IntRule%n)
  elseif (IntRule%flag==3 .or. IntRule%flag==4) then
    ErrFlag = GetVal(PropList,'nQuadAll',cTemp)
    read(cTemp,'(I6)') IntRule%nAll
  end if

  allocate(IntRule%nodes(IntRule%nAll,RCDIM))
  allocate(IntRule%weights(IntRule%nAll))
  call DefineIntegrationNodes(IntRule%flag,IntRule%n,IntRule.nAll, &
                              IntRule%nodes,IntRule%weights)

 end subroutine DefineIntegrationRule 

subroutine DefineIntegrationNodes(flag,n1,nAll,nodes,weights)
  use tools_module, only : kron1
  use nr, only : gauher
  implicit none
  ! 1 : Gauss-Hermite               1 : Working
  ! 2 : Sparse-Hermite              2 : NOT working
  ! 3 : pseudo-Montecarlo           3 : working 
  ! 4 : quasi-montecarlo            4 : NOT working
  ! 5 : taylor expand around 0      5 : NOT working
  ! 6 : Laplace                     6 : NOT working
  integer(i4b), intent(in)  :: flag   ! flag to select integration rule
  integer(i4b), intent(in)  :: n1(:)  ! number of nodes in each 1D rule
  integer(i4b), intent(in)  :: nall   ! total number of nodes in full multi-D rule
  real(dp),     intent(out) :: nodes(:,:),weights(:)
  integer(i4b)              :: i1
  real(dp), allocatable     :: e1(:),e2(:)
  real(dp), allocatable     :: temp1(:),temp2(:)
  real(dp), allocatable     :: x(:),w(:)

  ! variables needed to generate random numbers
  integer(i4b)              :: genid,subid,lseed,lstate,ifail
  integer(i4b), allocatable :: seed(:),state(:)
  integer(i4b)              :: mode,lr,ldx,ldc
  real(dp)                  :: mu(RCDIM)
  real(dp), allocatable     :: C(:,:),R(:)

  nodes   = 0.0d0
  weights = 1.0d0
  if (flag==1) then
    ! 0 = Gauss-Hermite
    do i1=1,RCDim
      allocate(x(n1(i1)))
      allocate(w(n1(i1)))
      call gauher(x,w)
      if (i1==1) then
        allocate(e1(nAll/n1(i1)))
        e1 = 1.0d0
        call kron1(x,e1,nodes(:,i1))
        call kron1(w,e1,weights)
        deallocate(e1)
     else if (i1>1 .and. i1<RCDIM) then
        allocate(e1(product(n1(1:i1-1))),e2(product(n1(i1+1:RCDIM))))
        allocate(temp1(product(n1(i1:RCDim))))
        allocate(temp2(nAll))
        e1=1.0d0
        e2=1.0d0
        temp1=0.0d0
        call kron1(x,e2,temp1)
        call kron1(e1,temp1,nodes(:,i1))
        call kron1(w,e2,temp1)
        call kron1(e1,temp1,temp2)
        weights = weights*temp2
        deallocate(e1,e2)
        deallocate(temp1,temp2)
      else if (i1==RCDIM) then
        allocate(e1(product(n1(1:i1-1))))
        allocate(temp1(nAll))
        e1 = 1.0d0
        call kron1(e1,x,nodes(:,i1))
        call kron1(e1,w,temp1)
        weights = weights*temp1
        deallocate(e1,temp1)
      end if
      deallocate(x,w)
    end do

  nodes=sqrt(2.0)*nodes
  weights=weights*pi_d**dble(-RCDim/2)
elseif (flag==2) then
!  % 1 = Sparse_hermite
!  SparseLevel = 2  ;
!  nQuadAll = sparse_grid_herm_size ( RCDim, SparseLevel );
!  [weights,nodes]=sparse_grid_herm(RCDim,SparseLevel,nQuadAll);
!  nodes = sqrt(2)*nodes';
!  weights = (pi^(-RCDim/2))*weights';
elseif (flag==3) then
  ! pseudo-MonteCarlo
  
  ! initialize random number generator
  genid = 3   ! Mersenne twister algorithm
  subid = 0   ! not used when genid==3
  lseed = 624 ! number of seeds needed for genid==3
  allocate(seed(lseed))
  seed = 1
  lstate = 633  ! min value for genid==3
  allocate(state(lstate))
  ifail = 0
  call G05KFF(genid,subid,seed,lseed,state,lstate,ifail)

  ! generate random numbers
  mode  = 2
  ifail = 0
  ldc   = RCDIM
  ldx   = nall
  mu    = 0.0d0
  allocate(C(LDC,RCDIM))
  C     = 0.0d0
  do i1=1,RCDIM
    C(i1,i1) =1.0d0
  end do
  LR = RCDIM*(RCDIM+1)+1
  allocate(R(LR))
  call G05RZF(mode,nall,RCDIM,MU,C,RCDIM,R,LR,state, &
              nodes,nAll,ifail)
  weights = 1.0d0/dble(nAll)
  deallocate(seed)
  deallocate(state)
  deallocate(C,R)
elseif (flag==5) then
  ! Laplace expansion
  !nodes    = 0;
  !weights  = 0;
end if

end subroutine DefineIntegrationNodes

end module GlobalModule
