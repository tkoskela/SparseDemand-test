! Revision history
! 2015AUG14  LN  add ReadWriteParameters and CopyParameters

module GlobalModule
  use ConstantsModule
  use PropertyList   ! tools for loading parameters from file
  implicit none

  type(Property)    :: PropList   ! list of keys in input file

  ! structure storing filenames
  type FilenameStructure
    character(len=99) :: mue,C,D
    character(len=99) :: InvCDiag,InvCOffDiag
    character(len=99) :: BC_beta,BC_CDiag,BC_COffDiag
    character(len=99) :: BD_beta,BD_CDiag,BD_COffDiag
    character(len=99) :: BD_Z,BC_Z
    character(len=99) :: sigp
    character(len=99) :: BasePriceFile
    character(len=99) :: TaxParmsFile
  end type
  type(FilenameStructure) :: ParmFiles

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
    integer(i4b) :: OutputFlag
    integer(i4b) :: SaveDataFlag
    integer(i4b) :: SimulateData
    integer(i4b) :: TestLikeFlag
    integer(i4b) :: TestIntegrationFlag
    integer(i4b) :: ComputeHessFlag
    integer(i4b) :: BICFlag
    integer(i4b) :: MPIFlag   ! 2 = each processor, 1 MC rep
                              ! 1 = each processor, subset of data
    integer(i4b) :: HotStart  ! 1 = load previous results from file
  end type

  type QuadRuleType
    integer(i4b), allocatable :: nQuad(:)   ! number of nodes in each dim.
    real(dp),     allocatable :: nodes(:,:) ! nAll x dim
    real(dp),     allocatable :: weights(:) ! nAll x 1
    integer(i4b)              :: nall
    integer(i4b)              :: flag       ! define integration rule
  end type

  type QuadRuleMultiType
    type(QuadRuleType), allocatable :: rule(:)
    integer(i4b)                    :: flag
    integer(i4b)                    :: nAll
    integer(i4b)                    :: nquadRaw
  end type

  type DataStructure
    integer(i4b)                :: NMC       ! number of MC replications
    integer(i4b)                :: N         ! number of households
    integer(i4b)                :: NSim      ! number of households in simulation
    integer(i4b)                :: M         ! number of markets
    real(dp), allocatable       :: e(:,:)    ! bliss points
    real(dp), allocatable       :: q(:,:)    ! q(1:d1,i1) = nonzero elements of quantity for i1
    real(dp), allocatable       :: p(:,:)    ! p(:,i1)    = prices for household i1
    real(dp), allocatable       :: expenditure(:)
    real(dp), allocatable       :: utility(:)
    integer(i4b),   allocatable :: market(:) ! market id
    integer(i4b),   allocatable :: iNonZero(:,:) ! iNonZero(1:d1,i1) = indexes of nonzero elements of q for i1
    integer(i4b),   allocatable :: iZero(:,:)    ! iZero(1:d3,i1)    = indexes of zero elements of q for i1
    integer(i4b),   allocatable :: nNonZero(:)   ! nNonZero(i1)      = number of goods with nonzero demand for i1
    character(20),  allocatable :: ColumnLabels(:) ! (J x 1) labels for columns of B
    real(dp),       allocatable :: eta(:,:)
    character(200)              :: RawDataFile
    character(8)                :: RawDataFormat
    integer(i4b),   allocatable :: HHID(:)
    integer(i4b),   allocatable :: shopid(:)
    integer(i4b),   allocatable :: fascia(:)
    character(len=12), allocatable :: fasciaChar(:)
    integer(i4b),   allocatable :: internet(:)
    integer(i4b),   allocatable :: SmallStore(:)
    integer(i4b),   allocatable :: date(:),day(:),month(:),week(:)
    integer(i4b)                :: nRawVars    ! number of variables in raw data file
    integer(i4b)                :: EstimationSeed
    integer(i4b)                :: SimulationSeed
    real(dp),       allocatable :: em_weights(:,:) ! weights for em algorithm
  end type

  type ParmsStructure
    character(Len=128)    :: file      ! name of file to save parms
    integer(i4b)          :: unit      ! unit number of parms file

    integer(i4b)          :: model    ! 1 = baseline, 2 = random coefficients
    integer(i4b)          :: K,J      ! dimensions of utility matrix
    integer(i4b)          :: nBC       ! size(BC)


    real(dp), allocatable :: B(:,:)   ! K x J matrix
    real(dp), allocatable :: D(:)     ! J x 1 vector, sqrt of diagonal of B for columns <=K
    real(dp), allocatable :: BC(:)     ! off-diagonal elements of B

    real(dp), allocatable :: MuE(:)      ! mean of e
    real(dp), allocatable :: mue_month(:,:) ! seasonal variation in mue
    real(dp), allocatable :: SIG(:,:)    ! variance of e
    real(dp), allocatable :: CSig(:,:)   ! cholesky decomp of SIG
    real(dp), allocatable :: InvC(:,:)   ! inverse of C
    real(dp), allocatable :: InvCDiag(:)    ! K x 1 sqrt of diag of InvC
    real(dp), allocatable :: InvCOffDiag(:) ! (K-1)*(K-2)/2 off-diagonal elements of InvC

    real(dp), allocatable :: GradBC(:,:)          ! (K x n) gradient of B w.r.t. C
    real(dp), allocatable :: GradBD(:,:)          ! (K x J) gradient of B w.r.t. D
    real(dp), allocatable :: GradInvCOffDiag(:,:) ! (K x n) gradient of InvC w.r.t. OffDiag
    real(dp), allocatable :: GradInvCDiag(:,:)    ! (K x K) gradient of InvC w.r.t. Diag

    ! Random B parameters
    !   B      = SphereToMatrix(C,D)
    !   log(D) = BD_Z*BD_beta + BD_C * eta
    !   C      = pi*normcdf(yC)
    !   yC     = BC_Z*BC_beta  + BC_C * eta
    !
    !   size(BD_Z)    = (J x BD_Z_DIM)
    !   size(BD_beta) = (BD_Z_DIM x 1)
    !   size(CD)    = (J x dim_eta)
    !   size(eta)  = (dim_eta x 1)
    !   size(yc)    = K*(K-1)/2 + (J-K) * (K-1)
    !               = nBC
    !   size(BC_Z)  = (nBC x BC_Z_DIM)
    !
    integer(i4b)          :: dim_eta      ! dimension of random coefficients in
                                          ! (BC,BD)
    integer(i4b)          :: BD_z_dim,BC_z_dim  ! dimension of product characteristics
    integer(i4b)          :: nBC_COffDiag   ! size of BC_COffDiag
    integer(i4b)          :: nBD_COffDiag   ! size of BD_COffDiag

    real(dp), allocatable :: BC_beta(:)   ! slope of y(c) w.r.t. etaZ
    real(dp), allocatable :: BD_beta(:)   ! log(D) = BD_beta* BD_z + BD_C* eta
    real(dp), allocatable :: BC_C(:,:)    ! matrix of factor loadings in BC: (nBC x dim_eta)
    real(dp), allocatable :: BD_C(:,:)    ! matrix of factor loadings in BD: (J x dim_eta)
    real(dp), allocatable :: BC_z(:,:)    ! product characteristics used in C    (nBC x BC_DIM_Z)
    real(dp), allocatable :: BD_z(:,:)    ! product characteristics used in D    (J   x BD_DIM_Z)

    real(dp), allocatable :: BC_CDiag(:)    ! diagonal elements of BC_C      (nBC x 1)
    real(dp), allocatable :: BC_COffDiag(:) ! off-diagonal elements of BC_C  (nBC_COffDiag  x 1)
    real(dp), allocatable :: BD_CDiag(:)    ! diagonal elements of BD_C      (J x 1)
    real(dp), allocatable :: BD_COffDiag(:) ! off-diagonal elements of BD_C  (nBD_COffDiag  x 1)

    real(dp), allocatable :: GradBD_C_CDiag(:,:)
    real(dp), allocatable :: GradBD_C_COffDiag(:,:)
    real(dp), allocatable :: GradBC_C_CDiag(:,:)
    real(dp), allocatable :: GradBC_C_COffDiag(:,:)

    real(dp), allocatable :: SigP(:,:)

    ! Month dummies
    real(dp), allocatable :: BD_month(:,:)

    ! (lower,upper) bounds on parameters
    real(dp), allocatable :: D_L(:),D_H(:)
    real(dp), allocatable :: BC_L(:),BC_H(:)
    real(dp), allocatable :: MUE_L(:),MUE_H(:)
    real(dp), allocatable :: InvCDiag_L(:),InvCDiag_H(:)
    real(dp), allocatable :: InvCOffDiag_L(:),InvCOffDiag_H(:)

    real(dp)              :: BC_lo,BC_hi  ! BC(i1) = bc_lo + (bc_hi-bc_lo) * pi
                                          !                 * normcdf(y)
    real(dp)              :: InvCDiag_LO  ! BL(InvCDiag) = InvCDiag_LO
    real(dp)              :: InvCDiag_HI  ! BU(InvCDiag) = x + InvCDiag_HI
    real(dp)              :: InvCOffDiag_LO !BL(InvCOffDiag) = pi * InvCOffDiag_LO
    real(dp)              :: InvCOffDiag_HI !BL(InvCOffDiag) = pi * InvCOffDiag_HI
    real(dp)              :: BC_beta_lo   ! lower bound
    real(dp)              :: BC_beta_hi   ! upper bound
    real(dp)              :: BC_CDiag_lo   ! lower bound
    real(dp)              :: BC_CDiag_hi   ! upper bound
    real(dp)              :: BD_beta_lo   ! lower bound
    real(dp)              :: BD_beta_hi   ! upper bound
    real(dp)              :: BD_CDiag_lo   ! lower bound
    real(dp)              :: BD_CDiag_hi   ! upper bound

    real(dp)              :: BD_month_lo  ! lower bound
    real(dp)              :: BD_month_hi  ! upper bound

    ! parameters used to determine details of analysis of results
    integer(i4b) :: nPrices_plot  ! number of prices to plot when demand plotting
  end type

  !  structure containing indexes
  !  defines map between structural parameters and x
  !  x = argument of likelihood function
  type SelectFreeType
    integer(i4b), allocatable :: D(:)
    integer(i4b), allocatable :: xD(:)
    integer(i4b), allocatable :: BC(:)
    integer(i4b), allocatable :: xBC(:)
    integer(i4b), allocatable :: MuE(:)
    integer(i4b), allocatable :: xMuE(:)
    integer(i4b), allocatable :: mue_month(:)
    integer(i4b), allocatable :: xmue_month(:)
    integer(i4b), allocatable :: InvCDiag(:)
    integer(i4b), allocatable :: xInvCDiag(:)
    integer(i4b), allocatable :: InvCOffDiag(:)
    integer(i4b), allocatable :: xInvCOffDiag(:)

    integer(i4b), allocatable :: BD_beta(:)
    integer(i4b), allocatable :: xBD_beta(:)
    integer(i4b), allocatable :: BD_CDiag(:)
    integer(i4b), allocatable :: xBD_CDiag(:)
    integer(i4b), allocatable :: BD_COffDiag(:)
    integer(i4b), allocatable :: xBD_COffDiag(:)
    integer(i4b), allocatable :: BC_beta(:)
    integer(i4b), allocatable :: xBC_beta(:)
    integer(i4b), allocatable :: BC_CDiag(:)
    integer(i4b), allocatable :: xBC_CDiag(:)
    integer(i4b), allocatable :: BC_COffDiag(:)
    integer(i4b), allocatable :: xBC_COffDiag(:)

    integer(i4b), allocatable :: BD_month(:)
    integer(i4b), allocatable :: xBD_month(:)

    ! first free element in each vector
    !  e.g. iFree%mue = (/iFree%mue1:parms%k/)
    integer(i4b)              :: mue1,mue_month1,invcdiag1,invcoffdiag1
    integer(i4b)              :: bd_beta1,bd_cdiag1,bd_coffdiag1
    integer(i4b)              :: bc_beta1,bc_cdiag1,bc_coffdiag1
    integer(i4b)              :: bd_month1

    integer(i4b)              :: nD
    integer(i4b)              :: nBC
    integer(i4b)              :: nMuE
    integer(i4b)              :: nMue_month
    integer(i4b)              :: nInvCDiag
    integer(i4b)              :: nInvCOffDiag

    integer(i4b)              :: nBD_beta
    integer(i4b)              :: nBD_CDiag
    integer(i4b)              :: nBD_COffDiag

    integer(i4b)              :: nBD_month

    integer(i4b)              :: nBC_beta
    integer(i4b)              :: nBC_CDiag
    integer(i4b)              :: nBC_COffDiag

    integer(i4b)              :: nAll

    integer(i4b)              :: flagD
    integer(i4b)              :: flagBC
    integer(i4b)              :: flagMUE
    integer(i4b)              :: flagMUE_month
    integer(i4b)              :: flagInvCDiag
    integer(i4b)              :: flagInvCOffDiag

    integer(i4b)              :: flagBC_beta
    integer(i4b)              :: flagBC_CDiag
    integer(i4b)              :: flagBC_COffDiag
    integer(i4b)              :: flagBD_beta
    integer(i4b)              :: flagBD_CDiag
    integer(i4b)              :: flagBD_COffDiag
    integer(i4b)              :: flagBD_month

    integer(i4b)              :: nPerIter  ! number free per iteration
    integer(i4b)              :: RandFlag  ! 1 for random directions
                                           ! 0 for non-random directions
    integer(i4b)              :: seed

    character(len=20), allocatable :: xlabels(:)
    integer(i4b)              :: nHess0,nHess1 ! used by ComputeHess2
                                               ! determine block of hessian
                                               ! to compute
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
    ! MaxAlgorithm = 1   E04WDF  : Dense problem
    ! MaxAlgorithm = 2   E04VHF  : Sparse problem
    integer(i4b)       :: Algorithm
    real(dp)           :: AbsTol     ! absolute tolerance for D01ESF Bayesian computations
    real(dp)           :: RelTol     ! relative tolerance for D01ESF Bayesian computations
    integer(i4b)           :: MaxLevel   ! maximum level for D01ESF Bayesian computations
    character(len=200) :: OptionsFile   ! file with E04WDF options
    character(len=200) :: BasisFile
    character(len=200) :: BackupBasisFile
    character(len=200) :: OldBasisFile
    integer(i4b)       :: SaveBasis  ! 1 to save basis info in BasisFile
    integer(i4b)       :: LoadBasis  ! 1 to load basis info from OldBasisFile
    integer(i4b)       :: SaveBasisFreq  ! save basis every i1 iterations
    real(dp)           :: DeltaX     ! finite difference for gradient approx.
    real(dp)           :: em_tol     ! tolerance for EM convergence
  end type

  type BayesType
    integer(i4b)    :: nAll
    real(dp), allocatable :: x(:,:)   ! x(:,i1) = draw number i1 from density of x
    real(dp), allocatable :: w(:)     ! integration weights
    real(dp), allocatable :: prior(:) ! prior evaluated at x(:,i1)

    real(dp)              :: BD_beta_lo,BD_beta_hi
    real(dp)              :: BD_month_lo,BD_month_hi
    real(dp)              :: BD_CDiag_lo,BD_CDiag_hi
    real(dp)              :: BD_COffDiag_lo,BD_COffDiag_hi

    real(dp)              :: BC_beta_lo,BC_beta_hi
    real(dp)              :: BC_CDiag_lo,BC_CDiag_hi
    real(dp)              :: BC_COffDiag_lo,BC_COffDiag_hi

    real(dp)              :: MUE_lo,MUE_hi
    real(dp)              :: MUE_month_lo,MUE_month_hi
    real(dp)              :: InvCDiag_lo,InvCDiag_hi
    real(dp)              :: InvCOffDiag_lo,InvCOffDiag_hi
  end type

  type(BayesType)             :: bayes ! information for Bayesian computations

  type(FlagStructure)         :: ControlOptions
  type(MaxStructure)          :: MaxOptions

  type(PenaltyStructureType)  :: Penalty
  type(QuadRuleMultiType), allocatable :: RandomE(:)  ! RandomE(d) = rule for d-dimensional rule for e
  type(QuadRuleType),  allocatable     :: RandomB(:)   ! MC rule: RandomB(i1) = quad rule of household i1
                                                   ! Gauss rule: RandomB(1) = quad rule
  type(SelectFreeType)        :: iFree
  type(ParmsStructure)        :: parms,parms0  ! parms  = current parameters
                                               ! parms0 = true parmaeters
  type(DataStructure)         :: HHData

  ! Output directory
  character(len=200)          :: OutDir
  character(len=200)          :: InputDir

  character(len=5), parameter :: OutParmsDir = 'parms'
  character(len=5), parameter :: TestFileName = 'a.txt'

  real(dp)                    :: inf
  real(dp)                    :: small

  integer(i4b), parameter     :: MasterID=0
  integer(i4b)                :: nWorkers
  integer(i4b), allocatable   :: G05State(:) ! random number generator state vector

  ! Structure holding gradient of DensityFunc2 w.r.t. parameters
  type DensityGradType
    ! (S1,mu1Tilde,SigmaTilde12,V,C2,Q,COmega11)
    ! z1 = inv(S1) * VT * p1 + S1*VT*q1 - mu1Tilde - SigmaTilde12 * inv(C2^T) * Q^T * x
    real(dp), allocatable :: S1(:,:)  ! gradient w.r.t. S1  : (nx x d1)
    real(dp), allocatable :: mu1Tilde(:,:) ! gradient w.r.t. nu1 : (nx x d1)
    real(dp), allocatable :: SigmaTilde12(:,:,:) ! gradient w.r.t. SigmaTilde12
                                             ! (nx x d1 x d2)
    real(dp), allocatable :: VT(:,:,:)       ! (nx x nV x nV)  (nV = size(V))
    real(dp), allocatable :: C2(:,:,:)       ! (nx x d2 x d2)
    real(dp), allocatable :: Q(:,:,:)        ! (nx x d2 x d2)
    real(dp), allocatable :: COmega11(:,:,:) ! (nx x d1 x d1)
  end type

contains
! 1) 141 subroutine AllocateGlobalVariables
! 2) 183 subroutine AllocatePenaltyParameters
! 3) 279 subroutine ComputeNMC
! 4) 312 subroutine DeallocateGlobalVariables
! 5) 342 subroutine DeallocatePenaltyParameters
! 6) 354 subroutine GetInputFile(InputFile)
! 7) 377 subroutine InitializeParameters(InputFile)
! 8) 479 subroutine ReadParameters
! 9) 547 subroutine DefineIntegrationRule
! 10) 589 subroutine DefineIntegrationNodes
! 11) 767 subroutine BroadCastParameters
!------------------------------------------------------------------------------
  subroutine AllocateGlobalVariables
    implicit none
    integer(i4b)    :: K,J,nc

    K = parms%K
    J = parms%J
    parms%nBC = K*(K-1)/2 + (K-1)*(J-K)
    parms%nBC_COffDiag = parms%dim_eta*(parms%dim_eta-1)/2 &
      + (parms%nBC - parms%dim_eta) * (parms%dim_eta-1)
    parms%nBD_COffDiag = parms%dim_eta*(parms%dim_eta-1)/2 &
      + (parms%J - parms%dim_eta) * (parms%dim_eta-1)

    ! allocate memory for parms
    call AllocateParms(parms)

  end subroutine AllocateGlobalVariables

  subroutine AllocateParms(LocalParms)
    implicit none
    type(ParmsStructure), intent(inout) :: LocalParms

    integer(i4b) :: K,J,nBC,nC

    K   = LocalParms%K
    J   = LocalParms%J
    nBC = LocalParms%nBC
    nC  = K*(K-1)/2

    ! B  (K x J)      matrix in utility
    ! D  (J x 1)      norm of each column of B
    ! BC (nBC x 1)  off-diagonal elements of B
    allocate(LocalParms%B(K,J),LocalParms%D(J))
    allocate(LocalParms%BC(LocalParms%nBC))

    ! mean and variance of e
    allocate(LocalParms%MuE(K))
    allocate(LocalParms%mue_month(K,12))
    allocate(LocalParms%SIG(K,K))
    allocate(LocalParms%CSig(K,K))
    allocate(LocalParms%InvC(K,K))
    allocate(LocalParms%InvCDiag(K))
    allocate(LocalParms%InvCOffDiag(nC))

    ! gradient of B and InvC w.r.t. (phi,D,phiC,DC)
    allocate(LocalParms%GradBC(LocalParms%K,LocalParms%nBC))
    allocate(LocalParms%GradBD(LocalParms%K,LocalParms%J))
    allocate(LocalParms%GradInvCOffDiag(LocalParms%K,nC))
    allocate(LocalParms%GradInvCDiag(LocalParms%K,LocalParms%K))

    if (LocalParms%model>=2) then

      ! parameters used to define BC_C and its gradient
      ! size(BC_C)              = nBC x dim_eta
      ! size(BC_CDiag)          = nBC x 1
      ! size(BC_COffDiag)       = nBC_COffDiag x 1
      ! size(GradBC_C_CDiag)    = dim_eta x nBC   (tranpose to get BC_C)
      ! size(GradBC_C_COffDiag) = dim_eta x nBC_COffDiag
      allocate(LocalParms%BC_beta(LocalParms%BC_z_dim))
      allocate(LocalParms%BC_CDiag(LocalParms%nBC))
      allocate(LocalParms%BC_COffDiag(LocalParms%nBC_COffDiag))
      allocate(LocalParms%GradBC_C_CDiag(LocalParms%dim_eta,LocalParms%nBC))
      allocate(LocalParms%GradBC_C_COffDiag(LocalParms%dim_eta,LocalParms%nBC_COffDiag))
      allocate(LocalParms%BC_C(LocalParms%nBC,LocalParms%dim_eta))
      allocate(LocalParms%BC_Z(LocalParms%nBC,LocalParms%BC_z_dim))

      ! parameters used to define BD_C and its gradient
      ! size(BD_C)              = J x dim_eta
      ! size(BD_CDiag)          = J x 1
      ! size(BD_COffDiag)       = nBD_COffDiag x 1
      ! size(GradBD_C_CDiag)    = dim_eta x J   (tranpose to get BD_C)
      ! size(GradBD_C_COffDiag) = dim_eta x nBD_COffDiag
      allocate(LocalParms%BD_beta(LocalParms%BD_z_dim))
      allocate(LocalParms%BD_month(LocalParms%J,12))
      allocate(LocalParms%BD_CDiag(LocalParms%J))
      allocate(LocalParms%BD_COffDiag(LocalParms%nBD_COffDiag))
      allocate(LocalParms%GradBD_C_CDiag(LocalParms%dim_eta,LocalParms%J))
      allocate(LocalParms%GradBD_C_COffDiag(LocalParms%dim_eta,LocalParms%nBD_COffDiag))
      allocate(LocalParms%BD_C(LocalParms%J,LocalParms%dim_eta))
      allocate(LocalParms%BD_Z(LocalParms%J,LocalParms%BD_z_dim))

      allocate(LocalParms%sigp(LocalParms%J,LocalParms%J))
    end if
  end subroutine AllocateParms

  ! allocate memory for data
  subroutine AllocateHHData
    implicit none

    allocate(HHData%q(parms%K,HHData%N))
    allocate(HHData%p(parms%J,HHData%N))
    allocate(HHData%e(parms%K,HHData%N))
    allocate(HHData%expenditure(HHData%N))
    allocate(HHData%utility(HHData%N))
    allocate(HHData%market(HHData%N))
    allocate(HHData%iNonZero(parms%K,HHData%N))
    allocate(HHData%iZero(parms%J,HHData%N))
    allocate(HHData%nNonZero(HHData%N))
    allocate(HHData%ColumnLabels(parms%J))
    allocate(HHData%HHID(HHData%N))
    allocate(HHData%date(HHData%N))
    allocate(HHData%day(HHData%N))
    allocate(HHData%month(HHData%N))
    allocate(HHData%week(HHData%N))
    allocate(HHData%shopid(HHData%N))
    allocate(HHData%fascia(HHData%N))
    allocate(HHData%fasciaChar(HHData%N))
    allocate(HHData%internet(HHData%N))
    allocate(HHData%SmallStore(HHData%N))

    if (parms%model>=2) then
      allocate(HHData%eta(parms%dim_eta,HHData%N))
    end if
  end subroutine AllocateHHData

  ! allocate memory for data
  subroutine AllocateLocalHHData(LocalHHData)
    implicit none
    type(DataStructure), intent(inout) :: LocalHHData

    allocate(LocalHHData%q(parms%K,LocalHHData%N))
    allocate(LocalHHData%p(parms%J,LocalHHData%N))
    allocate(LocalHHData%e(parms%K,LocalHHData%N))
    allocate(LocalHHData%expenditure(LocalHHData%N))
    allocate(LocalHHData%utility(LocalHHData%N))
    allocate(LocalHHData%market(LocalHHData%N))
    allocate(LocalHHData%iNonZero(parms%K,LocalHHData%N))
    allocate(LocalHHData%iZero(parms%J,LocalHHData%N))
    allocate(LocalHHData%nNonZero(LocalHHData%N))
    allocate(LocalHHData%ColumnLabels(parms%J))
    allocate(LocalHHData%HHID(LocalHHData%N))
    allocate(LocalHHData%date(LocalHHData%N))
    allocate(LocalHHData%day(LocalHHData%N))
    allocate(LocalHHData%month(LocalHHData%N))
    allocate(LocalHHData%week(LocalHHData%N))
    allocate(LocalHHData%shopid(LocalHHData%N))
    allocate(LocalHHData%fascia(LocalHHData%N))
    allocate(LocalHHData%fasciaChar(LocalHHData%N))
    allocate(LocalHHData%internet(LocalHHData%N))
    allocate(LocalHHData%SmallStore(LocalHHData%N))

    if (parms%model>=2) then
      allocate(LocalHHData%eta(parms%dim_eta,LocalHHData%N))
    end if
end subroutine AllocateLocalHHData

!------------------------------------------------------------------------------
  subroutine AllocatePenaltyParameters(pid)
#if USE_MPI==1
    use mpi
#endif
    implicit none

    integer(i4b), intent(in) :: pid
    integer(i4b) :: method,nx,nx1,nx2,nxP
    character(len=PL_VAL_LEN)       :: cTemp      ! temporary character string
    integer(i4b) :: ErrFlag,ierr,ix
    real(dp)     :: lambda

    ! method = 1   max L(x1,x2) + P(x2)  subject to   x2 = xPLus-xMinus
    !                                                 0 <= xPlus
    !                                                 0 <= xMinus

    ! read data from parameter file
    !    (if pid==MasterID)
    if (pid==MasterID) then
      ErrFlag = GetVal(PropList,'method',cTemp)
      read(cTemp,'(i2)') method

      ErrFlag = GetVal(PropList,'nx1',cTemp)
      read(cTemp,'(i4)') nx1

      ErrFlag = GetVal(PropList,'lambda',cTemp)
      read(cTemp,'(f11.0)') lambda

      ! number of values of lambda to try
      ErrFlag = GetVal(PropList,'nLambda',cTemp)
      read(cTemp,'(i4)') Penalty%nLambda

      ! minimum value of lambda
      ErrFlag = GetVal(PropList,'MinLambda',cTemp)
      read(cTemp,'(f11.0)') Penalty%MinLambda

      ! max value of lambda
      ErrFlag = GetVal(PropList,'MaxLambda',cTemp)
      read(cTemp,'(f11.0)') Penalty%MaxLambda
    endif

! broadcast information from MasterID to all processors
#if USE_MPI==1
  call mpi_bcast(method,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(nx1,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(Penalty%nlambda,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(nx1,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(lambda,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(Penalty%MinLambda,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(Penalty%MaxLambda,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr)
#endif

    if (parms%model==1) then
      nx  = iFree%nD + iFree%nBC + iFree%nMUE + iFree%nINVCDiag +  &
            iFree%nINVCOffDiag
    else if  (parms%model>=2) then
      nx = iFree%nBC_beta + iFree%nBC_CDiag + iFree%nBC_COffDiag &
         + iFree%nBD_beta + iFree%nBD_CDiag + iFree%nBD_COffDiag &
         + iFree%nMuE + iFree%nMue_month + iFree%nInvCDiag + iFree%nInvCOffDiag     &
         + iFree%nBD_month
    end if

    nx2 = nx-nx1
    nxP = nx + 2*nx2

    Penalty%method    = method
    Penalty%nx        = nx
    Penalty%nx1       = nx1
    Penalty%nx2       = nx2
    Penalty%nxP       = nxP
    Penalty%lambda    = lambda

    ! x_index1 = elements that are NOT penalized
    allocate(Penalty%x_index1(nx1))
    allocate(Penalty%xP_index1(nx1))

    Penalty%x_index1 = (/(ix,ix=1,nx1)/)
    Penalty%xP_index1 = Penalty%x_index1

    ! x_index2 = elements that are penalized
    allocate(Penalty%x_index2(nx2),Penalty%xP_index2(nx2))
    Penalty%x_index2  = nx1+(/(ix,ix=1,nx2)/)
    Penalty%xP_index2 = Penalty%x_index2

    allocate(Penalty%xP_index_xPlus(nx2))
    Penalty%xP_index_xPlus = nx + (/(ix,ix=1,nx2)/)

    allocate(Penalty%xP_index_xMinus(nx2))
    Penalty%xP_index_xMinus = nx+nx2+(/(ix,ix=1,nx2)/)

    allocate(Penalty%xPlus_index(nx2))
    allocate(Penalty%xMinus_index(nx2))
    Penalty%xPlus_index = (/(ix,ix=1,nx2)/)
    Penalty%xMinus_index = nx2+(/(ix,ix=1,nx2)/)

    allocate(Penalty%VectorLambda(Penalty%nLambda))
    Penalty%VectorLambda = Penalty%MinLambda + &
                           (Penalty%MaxLambda-Penalty%MinLambda)*dble((/(ix,ix=0,Penalty%nLambda-1)/))/dble(Penalty%nLambda-1)
  end subroutine AllocatePenaltyParameters

subroutine ComputeNMC(pid,NMC,IMC1,IMC2)
  implicit none
  integer(i4b), intent(in)  :: pid,NMC
  integer(i4b), intent(out) :: IMC1,IMC2
  integer(i4b)              :: NMC1,NMC2
  integer(i4b)              :: N1,N2,ix
  integer(i4b), allocatable :: j1(:),j2(:)
  ! NMC = total number of reps
  ! NMC1  = number of reps per processor for first N1 processors
  ! NMC2  = number of reps per processor for second N2 processors
  ! N1*NMC1 + N2*NMC2 = NMC
  ! NMC2 = NMC1+1
  ! N1 + N2 = NPROC
  NMC1 = NMC/(nWorkers+1)  ! 8/3 = 2
  NMC2 = NMC+1
  N1 = NMC-(nWorkers+1)*NMC1
  N2 = nWorkers+1-N1

  allocate(j1(N1),j2(N2))
!  J1 = (/1:1+(N1-1)*NMC1:NMC1/)
  J1 = (/(ix,ix=1,(N1-1)*NMC1+1,NMC1)/)
  if (pid+1<=N1) then
    IMC1 = j1(pid+1)
    IMC2 = IMC1+NMC1-1;
  end if
!  J2 = (/1:1+(N2-1)*NMC2:NMC2/)+N1*NMC1
  J2 = (/(ix,ix=1,NMC2*(N2-1)+1,NMC2)/) + NMC1*N1
  if (pid+1>N1) then
    IMC1 = j2(pid+1-N1)
    IMC2 = IMC1+NMC2-1
  end if
  deallocate(j1,j2)
end subroutine ComputeNMC

subroutine DeallocateParms(LocalParms)
  implicit none
  type(ParmsStructure), intent(inout) :: LocalParms

  if (allocated(LocalParms%B)) then
    deallocate(LocalParms%B)
  end if

  if (allocated(LocalParms%D)) then
    deallocate(LocalParms%D)
  end if

  if (allocated(LocalParms%BC)) then
    deallocate(LocalParms%BC)
  end if

  if (allocated(LocalParms%MUE)) then
    deallocate(LocalParms%MUE)
  end if

  if (allocated(LocalParms%MUE_month)) then
    deallocate(LocalParms%MUE_month)
  end if

  if (allocated(LocalParms%SIG)) then
    deallocate(LocalParms%SIG)
  end if

  if (allocated(LocalParms%csig)) then
    deallocate(LocalParms%csig)
  end if

  if (allocated(LocalParms%invc)) then
    deallocate(LocalParms%invc)
  end if

  if (allocated(LocalParms%InvCDiag)) then
    deallocate(LocalParms%InvCDiag)
  end if

  if (allocated(LocalParms%InvCOffDiag)) then
    deallocate(LocalParms%InvCOffDiag)
  end if

  if (allocated(LocalParms%GradBC)) then
    deallocate(LocalParms%GradBC)
  end if

  if (allocated(LocalParms%GradBD)) then
    deallocate(LocalParms%GradBD)
  end if

  if (allocated(LocalParms%GradInvCDiag)) then
    deallocate(LocalParms%GradInvCDiag)
  end if

  if (allocated(LocalParms%GradInvCOffDiag)) then
    deallocate(LocalParms%GradInvCOffDiag)
  end if

  if (allocated(LocalParms%BC_beta)) then
    deallocate(LocalParms%BC_beta)
  end if
  if (allocated(LocalParms%BC_CDiag)) then
    deallocate(LocalParms%BC_CDiag)
  end if
  if (allocated(LocalParms%BC_COffDiag)) then
    deallocate(LocalParms%BC_COffDiag)
  end if
  if (allocated(LocalParms%GradBC_C_CDiag)) then
    deallocate(LocalParms%GradBC_C_CDiag)
  end if
  if (allocated(LocalParms%GradBC_C_COffDiag)) then
    deallocate(LocalParms%GradBC_C_COffDiag)
  end if
  if (allocated(LocalParms%BC_C)) then
    deallocate(LocalParms%BC_C)
  end if
  if (allocated(LocalParms%BC_Z)) then
    deallocate(LocalParms%BC_Z)
  end if
  if (allocated(LocalParms%BD_beta)) then
    deallocate(LocalParms%BD_beta)
  end if
  if (allocated(LocalParms%BD_month)) then
    deallocate(LocalParms%BD_month)
  end if
  if (allocated(LocalParms%BD_CDiag)) then
    deallocate(LocalParms%BD_CDiag)
  end if
  if (allocated(LocalParms%BD_COffDiag)) then
    deallocate(LocalParms%BD_COffDiag)
  end if
  if (allocated(LocalParms%GradBD_C_CDiag)) then
    deallocate(LocalParms%GradBD_C_CDiag)
  end if
  if (allocated(LocalParms%GradBD_C_COffDiag)) then
    deallocate(LocalParms%GradBD_C_COffDiag)
  end if
  if (allocated(LocalParms%BD_C)) then
    deallocate(LocalParms%BD_C)
  end if
  if (allocated(LocalParms%BD_Z)) then
    deallocate(LocalParms%BD_Z)
  end if
  if (allocated(LocalParms%sigp)) then
    deallocate(LocalParms%sigp)
  end if

end subroutine DeallocateParms

subroutine DeallocateLocalHHData(LocalHHData)
  implicit none
  type(DataStructure), intent(inout) :: LocalHHData

  deallocate(LocalHHData%q,LocalHHData%p,LocalHHData%iNonZero,LocalHHData%iZero,LocalHHData%nNonZero)
  deallocate(LocalHHData%expenditure,LocalHHData%utility)
  deallocate(LocalHHData%HHID,LocalHHData%shopid,LocalHHData%date,LocalHHData%day)
  deallocate(LocalHHData%month,LocalHHData%week)
  deallocate(LocalHHData%fascia,LocalHHData%internet,LocalHHData%SmallStore)
  deallocate(LocalHHData%fasciaChar)

  deallocate(LocalHHData%market,LocalHHData%e)
  deallocate(LocalHHData%ColumnLabels)
  if (allocated(LocalHHData%eta)) then
    deallocate(LocalHHData%eta)
  end if
end subroutine DeallocateLocalHHData

!------------------------------------------------------------------------------
subroutine DeallocateGlobalVariables
  implicit none
  integer(i4b) i1

  !-----------------------------------------------------------------------------
  ! Deallocate parms
  call DeallocateParms(parms)


  !--------------------------------------------------------------------------
  ! deallocate HHData
  deallocate(HHData%q,HHData%p,HHData%iNonZero,HHData%iZero,HHData%nNonZero)
  deallocate(HHData%HHID,HHData%shopid,HHData%date,HHData%day)
  deallocate(HHData%month,HHData%week)
  deallocate(HHData%fascia,HHData%internet,HHData%SmallStore)
  deallocate(HHData%fasciaChar)
  deallocate(HHData%market,HHData%e)
  deallocate(HHData%expenditure)
  deallocate(HHData%utility)
  deallocate(HHData%ColumnLabels)
  if (allocated(HHData%eta)) then
    deallocate(HHData%eta)
  end if


  ! Deallocate integration rule
  deallocate(RandomE)

  ! deallocate RandomB
  if (parms%model>=2) then
    deallocate(RandomB)
  end if

  call DeallocateIFree(iFree)

end subroutine DeallocateGlobalVariables

subroutine DeallocateIFree(LocalIFree)
  implicit none
  type(SelectFreeType), intent(inout) :: LocalIFree

  !------------------------------------------------------------------------------
  ! Deallocate iFree
  if (allocated(LocalIFree%D)) then
    deallocate(LocalIFree%D,LocalIFree%xD)
  end if
  if (allocated(LocalIFree%BC)) then
    deallocate(LocalIFree%BC,LocalIFree%xBC)
  end if
  if (allocated(LocalIFree%MuE)) then
    deallocate(LocalIFree%MuE,LocalIFree%xMuE)
  end if
  if (allocated(LocalIFree%MuE_month)) then
    deallocate(LocalIFree%MuE_month,LocalIFree%xMuE_month)
  end if
  if (allocated(LocalIFree%InvCDiag)) then
    deallocate(LocalIFree%InvCDiag,LocalIFree%xInvCDiag)
  end if
  if (allocated(LocalIFree%InvCOffDiag)) then
    deallocate(LocalIFree%InvCOffDiag,LocalIFree%xInvCOffDiag)
  end if

  if (allocated(LocalIFree%BD_beta)) then
    deallocate(LocalIFree%BD_beta)
    deallocate(LocalIFree%xBD_beta)
  end if

  if (allocated(LocalIFree%BC_beta)) then
    deallocate(LocalIFree%BC_beta)
    deallocate(LocalIFree%xBC_beta)
  end if

  if (allocated(LocalIFree%BD_CDiag)) then
    deallocate(LocalIFree%BD_CDiag)
    deallocate(LocalIFree%xBD_CDiag)
  end if

  if (allocated(LocalIFree%BD_COffDiag)) then
    deallocate(LocalIFree%BD_COffDiag)
    deallocate(LocalIFree%xBD_COffDiag)
  end if

  if (allocated(LocalIFree%BC_CDiag)) then
    deallocate(LocalIFree%BC_CDiag)
    deallocate(LocalIFree%xBC_CDiag)
  end if

  if (allocated(LocalIFree%BC_COffDiag)) then
    deallocate(LocalIFree%BC_COffDiag)
    deallocate(LocalIFree%xBC_COffDiag)
  end if

  if (allocated(LocalIFree%xlabels)) then
    deallocate(LocalIFree%xlabels)
  end if

  if (allocated(LocalIFree%BD_month)) then
    deallocate(LocalIFree%BD_month)
    deallocate(LocalIFree%xBD_month)
  end if
end subroutine DeallocateIFree

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
#ifndef  __GFORTRAN__
     use IFPORT       ! intel fortran portability library
#endif
    implicit none
    character(len=*), intent(inout) :: InputFile  ! Name of input file
    integer(i4b)                    :: ErrFlag    ! error flag when reading inputs
    character(len=PL_VAL_LEN)       :: cTemp      ! temporary character string
    character(len=PL_VAL_LEN), allocatable :: cTemp1(:) ! temporary character vector
    real(dp)                        :: dTemp      ! temporary double precision real
    integer(i4b), allocatable       :: iTemp(:)   ! temporary integer array
    integer(i4b)                    :: nRead      ! number of values read
    integer(i4b)                    :: i1,DirStatus
    character(Len=100)              :: TempDir,TempFile,cmd_string
    character(len=20)  :: fmt1

    ! read parameters from property file
    ErrFlag = LoadPropertyFile(PropList,InputFile)
    if (ErrFlag /= PL_OK) then
      print *,'Error. Could not load ',trim(InputFile)
      stop
    end if

    inf = huge(1.0d0)
!    small = 1.0d-50
    ErrFlag = GetVal(PropList,'small',cTemp)
    read(cTemp,'(f12.0)') small

    ! name of output directory
    ErrFlag = GetVal(PropList,'OutDir',OutDir)
    if (ErrFlag /= PL_OK) then
      print *, 'No output directory specified in input file.'
      stop
    end if

    ! check whether outputdirectory exists.
    ! If not, create.
    TempFile = trim(OutDir) // '/' // trim(TestFileName)
    open(unit=85,file=TempFile,action='write',iostat=DirStatus)
    close(85)
    if (DirStatus .ne. 0) then
      print *, 'Warning: Output directory does not exist.'
      print *, 'Creating ',trim(OutDir),'.'
      cmd_string = 'mkdir ' // trim(OutDir)
      DirStatus = system(cmd_string)
      if (DirStatus .ne. 0) then
        print *, 'Error: Failed to create output directory ', trim(OutDir), '.'
        stop
      end if
    end if

    ! name of input directory
    ErrFlag = GetVal(PropList,'InputDir',InputDir)
    if (ErrFlag /= PL_OK) then
      print *, 'No input directory specified in input file.'
      stop
    end if

   ! Raw data file
   ErrFlag = GetVal(PropList,'RawData_FILE',HHData%RawDataFile)

   ErrFlag = GetVal(PropList,'RawDataFormat',HHData%RawDataFormat)

   ! Estimation and simulation seeds
   ErrFlag = GetVal(PropList,'EstimationSeed',cTemp)
   read(cTemp,'(i4)') HHData%EstimationSeed

   ErrFlag = GetVal(PropList,'SimulationSeed',cTemp)
   read(cTemp,'(i4)') HHData%SimulationSeed

   ! file containing base price information to use in analysis
   ErrFlag = GetVal(PropList,'BasePriceFile',cTemp)
   ParmFiles%BasePriceFile = trim(InputDir) // '/' // trim(cTemp)

   ! file containing parameters for tax experiments
   ErrFlag = GetVal(PropList,'TaxParmsFile',cTemp)
   ParmFiles%TaxParmsFile = trim(InputDir) // '/' // trim(cTemp)

   ! number of prices to plot when plotting demand
   ErrFlag = GetVal(PropList,'nPrices_plot',cTemp)
   read(cTemp,'(i4)') parms%nPrices_plot

   ! Parameter filenames
   ErrFlag = GetVal(PropList,'MUE_FILE',ParmFiles%MUE)
   ErrFlag = GetVal(PropList,'C_FILE',ParmFiles%C)
   ErrFlag = GetVal(PropList,'D_FILE',ParmFiles%D)
   ErrFlag = GetVal(PropList,'INVCDIAG_FILE',ParmFiles%InvCDiag)
   ErrFlag = GetVal(PropList,'INVCOFFDIAG_FILE',ParmFiles%InvCOffDiag)
   ErrFlag = GetVal(PropList,'BC_BETA_FILE',ParmFiles%BC_beta)
   ErrFlag = GetVal(PropList,'BC_CDIAG_FILE',ParmFiles%BC_CDiag)
   ErrFlag = GetVal(PropList,'BC_COFFDIAG_FILE',ParmFiles%BC_COffDiag)
   ErrFlag = GetVal(PropList,'BD_BETA_FILE',ParmFiles%BD_beta)
   ErrFlag = GetVal(PropList,'BD_CDIAG_FILE',ParmFiles%BD_CDiag)
   ErrFlag = GetVal(PropList,'BD_COFFDIAG_FILE',ParmFiles%BD_COffDiag)
   ErrFlag = GetVal(PropList,'SIGP_FILE',ParmFiles%sigp)
   ErrFlag = GetVal(PropList,'BD_Z_FILE',ParmFiles%BD_Z)
   ErrFlag = GetVal(PropList,'BC_Z_FILE',ParmFiles%BC_Z)

   ! outputFlag: 0 do not save output
   !             1 save output
   ErrFlag = GetVal(PropList,'OutputFlag',cTemp)
   read(cTemp,'(i2)') ControlOptions%OutputFlag

   ! TestLikeFlag:  1 test likelihood
   !                0 maximise likelihood
   ErrFlag = GetVal(PropList,'TestLikeFlag',cTemp)
   read(cTemp,'(i2)') ControlOptions%TestLikeFlag

   ! ComputeHessFlag: 1 use ComputeHess2(...)
   !                  2 load GradLHH from file and then compute hess
   ErrFlag = GetVal(PropList,'ComputeHessFlag',cTemp)
   read(cTemp,'(i2)') ControlOptions%ComputeHessFlag

   ! TestIntegrationFlag : 0 do not test
   !                       1 test accuracy of integration
   ErrFlag = GetVal(PropList,'TestIntegrationFlag',cTemp)
   read(cTemp,'(i2)') ControlOptions%TestIntegrationFlag

  ! SaveDataFlag : 0 do not save data
  !                1 save data to file
  ErrFlag = GetVal(PropList,'SaveDataFlag',cTemp)
  read(cTemp,'(i2)') ControlOptions%SaveDataFlag

  ! SimulateDataFlag : 0 load data from disk
  !                    1 simulate data
  ErrFlag = GetVal(PropList,'SimulateDataFlag',cTemp)
  read(cTemp,'(i2)') ControlOptions%SimulateData

  ! BICFlag : 0 do not minimize BIC
  !           1 minimize BIC
  ErrFlag = GetVal(PropList,'BICFlag',cTemp)
  read(cTemp,'(i2)') ControlOptions%BICFlag

  ! MPIFlag : 2 = each processor, subset of MC reps
  !           1 = each processor, subset of data
  ErrFlag = GetVal(PropList,'MPIFlag',cTemp)
  read(cTemp,'(i2)') ControlOptions%MPIFlag

  ! HotStart 1 = load parameters from previous results
  !          0 = use default parameter values
  ErrFlag = GetVal(PropList,'HotStart',cTemp)
  read(cTemp,'(i2)') ControlOptions%HotStart

  ! set maximization options
  ErrFlag = GetVal(PropList,'MaxAlgorithm',cTemp)
  read(cTemp,'(i2)') MaxOptions%Algorithm

  ! Absolute tolerance for D01ESF Bayesian estimator
  ErrFlag = GetVal(PropList,'AbsoluteTolerance',cTemp)
  read(cTemp,'(f12.0)') MaxOptions%AbsTol

  ! Relative tolerance for D01ESF Bayesian estimator
  ErrFlag = GetVal(PropList,'RelativeTolerance',cTemp)
  read(cTemp,'(f12.0)') MaxOptions%RelTol

  ! Maximum level for D01ESF Bayesian estimator
  ErrFlag = GetVal(PropList,'MaxLevel',cTemp)
  read(cTemp,'(i2)') MaxOptions%MaxLevel

  ErrFlag = GetVal(PropList,'E04OptionsFile',cTemp)
  MaxOptions%OptionsFile = trim(InputDir) // '/' // trim(cTemp)

  ErrFlag = GetVal(PropList,'SaveBasisFlag',cTemp)
  read(cTemp,'(i2)') MaxOptions%SaveBasis

  ErrFlag = GetVal(PropList,'LoadBasisFlag',cTemp)
  read(cTemp,'(i2)') MaxOptions%LoadBasis

  ErrFlag = GetVal(PropList,'SaveBasisFreq',cTemp)
  read(cTemp,'(i5)') MaxOptions%SaveBasisFreq

  ErrFlag = GetVal(PropList,'BasisFile',cTemp)
  MaxOptions%BasisFile = trim(OutDir) // '/' // trim(cTemp)

  ErrFlag = GetVal(PropList,'OldBasisFile',cTemp)
  MaxOptions%OldBasisFile = trim(OutDir) // '/' // trim(cTemp)

  ErrFlag = GetVal(PropList,'BackupBasisFile',cTemp)
  MaxOptions%BackupBasisFile = trim(OutDir) // '/' // trim(cTemp)

  ErrFlag = GetVal(PropList,'DeltaX',cTemp)
  read(cTemp,'(f12.0)') MaxOptions%DeltaX

  ErrFlag = GetVal(PropList,'em_tol',cTemp)
  read(cTemp,'(f12.0)') MaxOptions%em_tol

  ! NMC : 0 no MC repetitions
  !     > 0 number of MC repetitions
  ErrFlag = GetVal(PropList,'NMC',cTemp)
  read(cTemp,'(i4)') HHData%NMC

  ! Define problem size:  (N,J,K)
  ! N     = sample size
  ! NSim  = sample size for simulations
  ! M     = number of markets
  ! J     = number of products
  ! K     = rank of demand system
  ErrFlag = GetVal(PropList,'N',cTemp)
  read(cTemp,'(i7)') HHData%N
  ErrFlag = GetVal(PropList,'NSIM',cTemp)
  read(cTemp,'(i7)') HHData%NSIM
  ErrFlag = GetVal(PropList,'M',cTemp)
  read(cTemp,'(i5)') HHData%M

  ErrFlag = GetVal(PropList,'nRawVars',cTemp)
  read(cTemp,'(i4)') HHData%nRawVars

  ! file for parms output
  TempDir = trim(OutDir) // '/' // trim(OutParmsDir)
  TempFile = trim(TempDir) // '/' // trim(TestFileName)
  open(unit=85,file=TempFile,action='write',iostat=DirStatus)
  close(85)
  if (DirStatus .ne. 0) then
    print *, 'Warning: ' // trim(OutParmsDir) // ' directory does not exist.'
    print *, 'Creating ',trim(TempDir),'.'
    cmd_string = 'mkdir ' // trim(tempDir)
    DirStatus = system(trim(cmd_string))
    if (DirStatus .ne. 0) then
      print *, 'Error: Failed to create output directory ', trim(TempDir), '.'
      stop
    end if
  end if
  ! output file for parameters
  parms%file  = trim(TempDir) // '/parms.csv'
  parms%unit = 51

  ErrFlag = GetVal(PropList,'J',cTemp)
  read(cTemp,'(i5)') parms%J
  ErrFlag = GetVal(PropList,'K',cTemp)
  read(cTemp,'(i5)') parms%K

  ErrFlag = GetVal(PropList,'model',cTemp)
  read(cTemp,'(i3)') parms%model

  if (parms%model>=2) then
    ErrFlag = GetVal(PropList,'BC_z_dim',cTemp)
    read(cTemp,'(i5)') parms%BC_z_dim
    if (parms%BC_z_dim==0) then
      parms%BC_z_dim = parms%J*(parms%K-1)-(parms%K*(parms%K-1))/2
    end if

    ErrFlag = GetVal(PropList,'BD_z_dim',cTemp)
    read(cTemp,'(i5)') parms%BD_z_dim
    parms%BD_z_dim = merge(parms%BD_z_dim,parms%J,parms%BD_z_dim>0)

    ErrFlag = GetVal(PropList,'BC_lo',cTemp)
    read(cTemp,'(f12.0)') parms%BC_lo

    ErrFlag = GetVal(PropList,'BC_hi',cTemp)
    read(cTemp,'(f12.0)') parms%BC_hi

    ErrFlag = GetVal(PropList,'InvCDiag_LO',cTemp)
    read(cTemp,'(f12.0)') parms%InvCDiag_LO

    ErrFlag = GetVal(PropList,'InvCDiag_HI',cTemp)
    read(cTemp,'(f12.0)') parms%InvCDiag_HI

    ErrFlag = GetVal(PropList,'InvCOffDiag_LO',cTemp)
    read(cTemp,'(f12.0)') parms%InvCOffDiag_LO

    ErrFlag = GetVal(PropList,'InvCOffDiag_HI',cTemp)
    read(cTemp,'(f12.0)') parms%InvCOffDiag_HI

    ErrFlag = GetVal(PropList,'BC_beta_lo',cTemp)
    read(cTemp,'(f12.0)') parms%BC_beta_lo

    ErrFlag = GetVal(PropList,'BC_beta_hi',cTemp)
    read(cTemp,'(f12.0)') parms%BC_beta_hi

    ErrFlag = GetVal(PropList,'BC_CDiag_lo',cTemp)
    read(cTemp,'(f12.0)') parms%BC_CDiag_lo

    ErrFlag = GetVal(PropList,'BC_CDiag_hi',cTemp)
    read(cTemp,'(f12.0)') parms%BC_CDiag_hi

    ErrFlag = GetVal(PropList,'BD_beta_lo',cTemp)
    read(cTemp,'(f12.0)') parms%BD_beta_lo

    ErrFlag = GetVal(PropList,'BD_beta_hi',cTemp)
    read(cTemp,'(f12.0)') parms%BD_beta_hi

    ErrFlag = GetVal(PropList,'BD_CDiag_lo',cTemp)
    read(cTemp,'(f12.0)') parms%BD_CDiag_lo

    ErrFlag = GetVal(PropList,'BD_CDiag_hi',cTemp)
    read(cTemp,'(f12.0)') parms%BD_CDiag_hi

    ! RandomB parameters: dimension of random elements of C and D
    ErrFlag = GetVal(PropList,'dim_eta',cTemp)
    read(cTemp,'(i3)') parms%dim_eta

    ! Month dummy parameters: upper and lower bounds
    ErrFlag = GetVal(PropList,'BD_month_lo',cTemp)
    read(cTemp,'(f12.0)') parms%BD_month_lo
    ErrFlag = GetVal(PropList,'BD_month_hi',cTemp)
    read(cTemp,'(f12.0)') parms%BD_month_hi

  end if

  ! allocate memory for (parms,QuadRule,iFree) and HHData
  call AllocateGlobalVariables
  call AllocateHHData

  ! Model 1:
  !     read in (parms%BC,parms%BD,parms%MuE,parms%InvCDiag,parms%InvCOffDiag)
  !
  ! Model 2:
  !     read in (betaC,CCDiag,CCOffDiag,betaD,CDDiag,CDOffDiag)
  !          and (parms%MuE,parms%InvCDiag,parms%InvCOffDiag)
  call ReadParameters

  ! Read in integration rule parameters
  allocate(RandomE(parms%K))
  ErrFlag = GetVal(PropList,'IntegrationFlag',cTemp)
  read(cTemp,'(i2)') RandomE(1)%flag
  RandomE(1:parms%K)%flag = RandomE(1)%flag

  ErrFlag = GetVal(PropList,'nQuadAll',cTemp)
  write(fmt1,'(a1,i2,a3)') '(',parms%K,'i4)'
  read(cTemp,fmt1) RandomE(1:parms%k)%nAll

  ErrFlag = GetVal(PropList,'nQuad',cTemp)
  write(fmt1,'(a1,i2,a3)') '(',parms%K,'i3)'
  read(cTemp,fmt1) RandomE(1:parms%k)%nQuadRaw

  if (parms%model>=2) then
    allocate(RandomB(1))
    ! Define RandomB rule
    ErrFlag = GetVal(PropList,'RandomB_flag',cTemp)
    read(cTemp,'(i2)') RandomB(1)%flag

    ErrFlag = GetVal(PropList,'RandomB_nall',cTemp)
    read(cTemp,'(i6)') RandomB(1)%nall

    allocate(RandomB(1)%nQuad(parms%dim_eta),source=3)
    ErrFlag = GetVal(PropList,'RandomB_nQuad',cTemp)
    nRead = min(parms%dim_eta,10)
    write(fmt1,'(a1,i2,a3)') '(',nread,'i3)'
    read(cTemp,fmt1) RandomB(1)%nQuad(1:nRead)

  end if ! (parms%model>=2) then

  ! initialize random number generator for:
  !   1) pseudo-MC integration rule
  !   2) simulate data for monte carlo study
  !   this subroutine sets G05State
  call InitializeG05

  ! read in flags to select free parameters
  if (parms%model==1) then
    ErrFlag = GetVal(PropList,'FreeFlagD',cTemp)
    read(cTemp,'(i2)') iFree%flagD
    ErrFlag = GetVal(PropList,'FreeFlagBC',cTemp)
    read(cTemp,'(i2)') iFree%flagBC
  end if
  ErrFlag = GetVal(PropList,'FreeFlagMUE',cTemp)
  read(cTemp,'(i2)') iFree%flagMUE
  ErrFlag = GetVal(PropList,'FreeFlagMUE_month',cTemp)
  read(cTemp,'(i2)') iFree%flagMUE_month

  ErrFlag = GetVal(PropList,'FreeFlagInvCDiag',cTemp)
  read(cTemp,'(i2)') iFree%flagInvCDiag
  ErrFlag = GetVal(PropList,'FreeFlagInvCOffDiag',cTemp)
  read(cTemp,'(i2)') iFree%flagInvCOffDiag

  ! read in (mue1,invcdiag1,invcoffdiag1)
  ErrFlag = GetVal(PropList,'free_mue1',cTemp)
  read(cTemp,'(i3)') iFree%mue1

  ErrFlag = GetVal(PropList,'free_mue_month1',cTemp)
  read(cTemp,'(i3)') iFree%mue_month1

  ErrFlag = GetVal(PropList,'free_cdiag1',cTemp)
  read(cTemp,'(i3)') iFree%invcdiag1

  ErrFlag = GetVal(PropList,'free_coffdiag1',cTemp)
  read(cTemp,'(i3)') iFree%invcoffdiag1

  if (parms%model>=2) then
    ErrFlag = GetVal(PropList,'FreeFlagBC_beta',cTemp)
    read(cTemp,'(i2)') iFree%flagBC_beta

    ErrFlag = GetVal(PropList,'FreeFlagBD_beta',cTemp)
    read(cTemp,'(i2)') iFree%flagBD_beta

    ErrFlag = GetVal(PropList,'FreeFlagBC_CDiag',cTemp)
    read(cTemp,'(i2)') iFree%flagBC_CDiag

    ErrFlag = GetVal(PropList,'FreeFlagBD_CDiag',cTemp)
    read(cTemp,'(i2)') iFree%flagBD_CDiag

    ErrFlag = GetVal(PropList,'FreeFlagBC_COffDiag',cTemp)
    read(cTemp,'(i2)') iFree%flagBC_COffDiag

    ErrFlag = GetVal(PropList,'FreeFlagBD_COffDiag',cTemp)
    read(cTemp,'(i2)') iFree%flagBD_COffDiag

    ErrFlag = GetVal(PropList,'FreeFlagBD_month',cTemp)
    read(cTemp,'(i2)') iFree%flagBD_month

    ! Read in (bc_beta1,bc_cdiag1,bc_coffdiag1)
    ! read in (bd_beta1,bd_cdiag1,bd_coffdiag1)
    ErrFlag = GetVal(PropList,'free_bc_beta1',cTemp)
    read(ctemp,'(i3)') iFree%bc_beta1

    ErrFlag = GetVal(PropList,'free_bc_cdiag1',cTemp)
    read(ctemp,'(i3)') iFree%bc_cdiag1

    ErrFlag = GetVal(PropList,'free_bc_coffdiag1',cTemp)
    read(ctemp,'(i3)') iFree%bc_coffdiag1

    ErrFlag = GetVal(PropList,'free_bd_beta1',cTemp)
    read(ctemp,'(i3)') iFree%bd_beta1

    ErrFlag = GetVal(PropList,'free_bd_cdiag1',cTemp)
    read(ctemp,'(i3)') iFree%bd_cdiag1

    ErrFlag = GetVal(PropList,'free_bd_coffdiag1',cTemp)
    read(ctemp,'(i3)') iFree%bd_coffdiag1

    ErrFlag = GetVal(PropList,'free_bd_month1',cTemp)
    read(ctemp,'(i3)') iFree%bd_month1

  end if

  ErrFlag = GetVal(PropList,'nPerIter',cTemp)
  read(cTemp,'(i2)') iFree%nPerIter

  ErrFlag = GetVal(PropList,'FreeRandFlag',cTemp)
  read(cTemp,'(i2)') iFree%RandFlag

  ErrFlag = GetVal(PropList,'FreeSeed',cTemp)
  read(cTemp,'(i10)') iFree%seed

  ! Determine which block of hessian to compute
  ! used in ComputeHess2
  ErrFlag = GetVal(PropList,'nHess0',ctemp)
  read(ctemp,'(i4)') iFree%nHess0
  ErrFlag = GetVal(PropList,'nHess1',ctemp)
  read(ctemp,'(i4)') iFree%nHess1

end subroutine InitializeParameters


!------------------------------------------------------------------------------
! read in initial values for (D,phi,MuE,InvC_D,InvC_phi)
! and for                    (betaD,CDDiag,CDOffDiag)
!                            (betaC,CCDiag,CCOffDiag)
  subroutine ReadParameters
    use NewTools, only : SphereToMatrix,MatrixInverse
    implicit none
    integer(i4b)       :: tempunit
    character(len=200) :: tempfile
    character(len=20)  :: fmt1

    integer(i4b), allocatable :: tempindex(:)
    integer(i4b)              :: row,col
    real(dp),     allocatable :: temp2(:,:)
    integer(i4b)              :: i1,n

    tempunit = 65

    if (parms%model==1) then
      ! D: norm of columns of B
      tempfile = trim(InputDir) // '/D.raw'
      open(UNIT = tempunit, &
           FILE = tempfile, &
           ACTION = 'read')
      parms%D = 0.0d0

      do i1=1,parms%J
        read(tempunit,'(g25.16)') parms%D(i1)
      end do
      close(tempunit)

      ! C : off-diagonal elements of B
      tempfile = trim(InputDir) // '/C.raw'
      open(UNIT = tempunit,  &
           FILE = tempfile,  &
           ACTION = 'read')
       parms%BC = 0.0d0
      do i1=1,parms%nBC
        read(tempunit,'(g25.16)') parms%BC(i1)
      end do
      close(tempunit)

      call SphereToMatrix(parms%BC,parms%D,parms%K,parms%J,parms%B)
    end if ! if (model==1)

    ! parms%MuE : mean of e
    tempfile = trim(InputDir) // '/' // trim(ParmFiles%MUE)
    open(UNIT = tempunit,  &
         FILE = tempfile,  &
         ACTION = 'read')
    parms%MUE = 0.0d0
    do i1=1,parms%K
      read(tempunit,'(g25.16)') parms%MUE(i1)
    end do
    close(tempunit)

    parms%mue_month = 0.0d0

    ! parms%InvCDiag : diagonal elements of inverse of C = chol(sig)
    tempfile = trim(InputDir) // '/' // trim(ParmFiles%InvCDiag)
    open(UNIT = tempunit,  &
         FILE = tempfile,  &
         ACTION = 'read')
    parms%InvCDiag = 0.0d0
    do i1=1,parms%K
      read(tempunit,'(g25.16)') parms%InvCDiag(i1)
    end do
    close(tempunit)

    ! parms%InvCOffDiag : off-diagonal elements of inverse of C = chol(sig)
    tempfile = trim(InputDir) // '/' // trim(ParmFiles%InvCOffDiag)
    open(UNIT = tempunit,  &
         FILE = tempfile,  &
         ACTION = 'read')
    parms%InvCOffDiag = 0.0d0
    do i1=1,parms%K*(parms%K-1)/2
      read(tempunit,'(g25.16)') parms%InvCOffDiag(i1)
    end do
    close(tempunit)
    call SphereToMatrix(parms%InvCOffDiag,parms%InvCDiag,parms%K,parms%K,parms%InvC)
    ! compute CSig = inv(InvC)
    parms%InvC = transpose(parms%InvC)
    call MatrixInverse(parms%InvC,parms%CSig,'Lower triangular')
    parms%sig  = matmul(parms%CSig,transpose(parms%CSig))

    if (parms%model>=2) then
      ! read in (BC_beta,BD_beta)
      tempfile = trim(InputDir) // '/' // trim(ParmFiles%BC_Beta)
      open(UNIT = tempunit,  &
           FILE = tempfile,  &
           ACTION = 'read')
      parms%BC_beta = 0.0d0
      do i1=1,parms%BC_z_dim
        read(tempunit,'(g25.16)') parms%BC_beta(i1)
      end do
      close(tempunit)

      tempfile = trim(InputDir) //  '/' // trim(ParmFiles%BD_Beta)
      open(UNIT = tempunit,  &
           FILE = tempfile,  &
           ACTION = 'read')
      parms%BD_beta = 0.0d0
      do i1=1,parms%BD_z_dim
        read(tempunit,'(g25.16)') parms%BD_beta(i1)
      end do
      close(tempunit)

      ! read (BC_CDiag,BD_CDiag)
      tempfile = trim(InputDir) // '/' // trim(ParmFiles%BC_CDiag)
      open(UNIT = tempunit,  &
           FILE = tempfile,  &
           ACTION = 'read')
      parms%BC_CDiag = 0.0d0
      do i1=1,parms%nBC
        read(tempunit,'(g25.16)') parms%BC_CDiag(i1)
      end do
      close(tempunit)

      tempfile = trim(InputDir) // '/' // trim(ParmFiles%BD_CDiag)
      open(UNIT = tempunit,  &
           FILE = tempfile,  &
           ACTION = 'read')
      parms%BD_CDiag = 0.0d0
      do i1=1,parms%J
        read(tempunit,'(g25.16)') parms%BD_CDiag(i1)
      end do
      close(tempunit)

      ! BD_month: default = 0.0
      parms%BD_month = 0.0d0

      ! read (BC_COffDiag,BD_COffDiag)
      tempfile = trim(InputDir) // '/' // trim(ParmFiles%BC_COffDiag)
      open(UNIT = tempunit,  &
           FILE = tempfile,  &
           ACTION = 'read')
      parms%BC_COffDiag = 0.0d0
      do i1=1,parms%nBC_COffDiag
        read(tempunit,'(g25.16)') parms%BC_COffDiag(i1)
      end do
      close(tempunit)

      tempfile = trim(InputDir) // '/' // trim(ParmFiles%BD_COffDiag)
      open(UNIT = tempunit,  &
           FILE = tempfile,  &
           ACTION = 'read')
      parms%BD_COffDiag = 0.0d0
      do i1=1,parms%nBD_COffDiag
        read(tempunit,'(g25.16)') parms%BD_COffDiag(i1)
      end do
      close(tempunit)

      ! Compute parms%BD_C
      ! size(BD_C) = (J x dim_eta)
      ! SphereToMatrix creates upper triangular matrix of size
      !    (dim_eta x J)
      ! BD_C is lower triangular:  we need the transpose as below.
      !     size(BD_C) = (J x dim_eta)
      allocate(temp2(parms%dim_eta,parms%J))
      call SphereToMatrix(parms%BD_COffDiag,parms%BD_CDiag, &
                          parms%dim_eta,parms%J,temp2, &
                          parms%GradBD_C_COffDiag,parms%GradBD_C_CDiag)
      parms%BD_C = transpose(temp2)
      deallocate(temp2)

      ! Compute parms%BC_C
      ! size(BC_C) = (nBC x dim_eta)
      ! SphereToMatrix creates upper triangular matrix of size
      !     (dim_eta x nBC)
      ! size(BC_C) = (nBC x dim_eta)
      ! after call, compute transpose to create lower triangular matrix
      allocate(temp2(parms%dim_eta,parms%nBC))
      call SphereToMatrix(parms%BC_COffDiag,parms%BC_CDiag, &
                          parms%dim_eta,parms%nBC,temp2, &
                          parms%GradBC_C_COffDiag,parms%GradBC_C_CDiag)
      parms%BC_C = transpose(temp2)
      deallocate(temp2)

      ! Load BD_Z (J x BD_Z_DIM) matrix
      tempfile = trim(InputDir) // "/" // trim(ParmFiles%BD_Z)
      open(unit=tempunit, &
           file=tempfile, &
           action='read')
      parms%BD_Z = 0.0d0

      write(fmt1,'(a1,i2,a7)') '(',parms%BD_Z_DIM,'g25.16)'
      do i1=1,parms%J
        read(tempunit,fmt1) parms%BD_Z(i1,:)
        !read(tempunit,1390) parms%BD_Z(i1,:)
      end do
!1390 format(<parms%BD_Z_DIM>f25.16)
      close(tempunit)

      ! Load parms%BC_Z   (nBC x BC_z_dim)
      !   c(j,k) = BC_Z * BC_beta + eta * BC_C
      parms%BC_z = 0.0d0
      tempfile = trim(InputDir) // "/" // trim(ParmFiles%BC_Z)
      open(unit=tempunit, &
           file=tempfile, &
           action='read')
      write(fmt1,'(a1,i3,a7)') '(',parms%BC_Z_DIM,'g25.16)'
      do i1=1,parms%nBC
        read(tempunit,fmt1) parms%BC_Z(i1,:)
      end do
!1420 format(<parms%BC_Z_DIM>f25.16)
      close(tempunit)
    end if ! if (model>=2) then

    tempfile = trim(InputDir) // '/' // trim(ParmFiles%sigp)
    open(UNIT = tempunit, &
         FILE = tempfile, &
         ACTION = 'read')

    parms%sigp = 0.0d0
    do i1=1,parms%J*parms%J
      row = (i1-1)/parms%J + 1
      col = i1-parms%J*(row-1)
      read(tempunit,897) parms%sigp(row,col)
    end do
    897 format(f25.0)
    close(tempunit)

end subroutine ReadParameters

subroutine CreateQuadRule(pid,nprocs)
  implicit none
  integer(i4b), intent(in) :: pid,nprocs
  integer(i4b)             :: N1,R,nskip,ndraws
  logical                  :: MC_EFlag,MC_BFlag

  ! Flags for MC vs Gausian rules
  MC_EFlag = (RandomE(1)%flag==2 .or. RandomE(1)%flag==3 .or. RandomE(1)%flag==7)
  MC_BFlag = .FALSE.
  if (parms%model>=2) then
    MC_BFlag = (RandomB(1)%flag==2 .or. RandomB(1)%flag==3 .or. RandomB(1)%flag==7)
  end if

  if (nprocs>1) then
    if (pid>MasterID) then
      N1 = HHData%N / nworkers
      R  = HHData%N - nworkers*N1
      ! ndraws = number of draws per worker
      ! nskip  = number of draws to skip ahead to ensure
      !          independent sequences across processors
      if (.not. MC_EFlag .and. .not. MC_BFlag) then
        ! Both rules gaussian
        ndraws = 0
      else if (MC_EFlag .and. .not. MC_BFlag) then
        ! Only RandomE is MC rule
        ndraws = RandomE(1)%nall * ((parms%K+1)*parms%K)/2
      else if (.not. MC_EFlag .and. MC_BFlag) then
        ! Only RandomB is MC rule
        ndraws = RandomB(1)%nall * parms%dim_eta
      else if (MC_EFlag .and. MC_BFlag) then
        ! Both rules are MC rules
        ndraws = RandomE(1)%nall * ((parms%K+1)*parms%K)/2 &
               + RandomB(1)%nall * parms%dim_eta
      end if
      nskip = ndraws*((pid-1)*N1+min(pid-1,R))
      N1 = merge(N1+1,N1,pid<=R)
      call DefineIntegrationRule(N1,pid,nskip,MC_EFlag,MC_BFlag)
    end if
  else
    N1    = HHData%N
    nskip = 0
    call DefineIntegrationRule(N1,pid,nskip,MC_EFlag,MC_BFlag)
  end if
end subroutine CreateQuadRule

subroutine DefineIntegrationRule(N1,pid,nskip,MC_EFlag,MC_BFlag)
  use nag_library, only : G05KJF
  implicit none
  integer(i4b), intent(in)   :: N1,pid,nskip
  logical,      intent(in)   :: MC_EFlag,MC_BFlag
  integer(i4b)               :: ErrFlag,flag,nall
  integer(i4b)               :: n,hh,k1
  integer(i4b), allocatable  :: nQuad(:)

  ! integration rule
  ! 0 = Gauss hermite            0 : not working
  ! 1 = Sparse Hermite           1 : not working
  ! 2 = Pseudo Monte Carlo       2 : not working
  ! 3 = Quasi Monte Carlo        3 : not working
  ! 6 = Gauss-Legendre on (-1,1) 6 : working
  ! 7 = pseudo-MC on (-1,1)      7 : working

  ! skip ahead random number generator to ensure all rules are independent
  call G05KJF(nskip,G05State,ErrFlag)

  if (.not. MC_EFlag) then
    ! RandomE rule is Gaussian rule
    do k1=1,parms%K
      ! Tensor product rules
      ! same rule for all households
      allocate(RandomE(k1)%rule(1))
      allocate(RandomE(k1)%rule(1)%nQuad(k1))
      RandomE(k1)%rule(1)%nQuad = RandomE(k1)%nQuadRaw
      nall = product(RandomE(k1)%rule(1)%nQuad)
      allocate(RandomE(k1)%rule(1)%nodes(nall,k1))
      allocate(RandomE(k1)%rule(1)%weights(nall))
      call DefineIntegrationNodes(k1,RandomE(k1)%flag,RandomE(k1)%rule(1)%nQuad,nall, &
                                  RandomE(k1)%rule(1)%nodes,RandomE(k1)%rule(1)%weights)
    end do
  else
    ! allocate memory for MC RandomE rule
    do k1=1,parms%K
      allocate(RandomE(k1)%rule(N1))
    end do
  end if

  if (.not. MC_BFlag) then
    ! RandomB rule is Gaussian rule
    RandomB(1)%nall =product(RandomB(1)%nQuad)
    allocate(nquad(parms%dim_eta))
    nquad = RandomB(1)%nquad
    allocate(RandomB(1)%nodes(RandomB(1)%nall,parms%dim_eta))
    allocate(RandomB(1)%weights(RandomB(1)%nall))
    call DefineIntegrationNodes(parms%dim_eta,RandomB(1)%flag, &
                                nquad, &
                                RandomB(1)%nall, &
                                RandomB(1)%nodes,RandomB(1)%weights)
  else
    ! allocate memory for RandomB MC rule
    flag = RandomB(1)%flag
    nall = RandomB(1)%nall
    allocate(nquad(parms%dim_eta))
    nquad = RandomB(1)%nquad
    deallocate(RandomB)
    allocate(RandomB(N1))
    RandomB%flag = flag
    RandomB%nall = nall
  end if


  if (MC_EFlag) then
    ! RandomE is MC rule
    do hh=1,N1
      do k1=1,parms%K
        ! Monte Carlo rules
        ! different rule for each household
        nall = RandomE(k1)%nall
        allocate(RandomE(k1)%rule(hh)%nodes(nall,k1))
        allocate(RandomE(k1)%rule(hh)%weights(nall))
        call DefineIntegrationNodes(k1,RandomE(k1)%flag,RandomE(k1)%rule(hh)%nQuad,nall, &
                                    RandomE(k1)%rule(hh)%nodes,RandomE(k1)%rule(hh)%weights)
      end do

      if (MC_BFlag) then
        ! RandomB is also MC rule
        nall = RandomB(hh)%nall
        allocate(RandomB(hh)%nodes(nall,parms%dim_eta))
        allocate(RandomB(hh)%weights(nall))
        call DefineIntegrationNodes(parms%dim_eta,RandomB(hh)%flag, &
                                    nquad, &
                                    nall, &
                                    RandomB(hh)%nodes,RandomB(hh)%weights)
      end if
    end do
  elseif (.not. MC_EFlag .and. MC_BFlag) then
    ! Only Random B is MC rule
    do hh=1,N1
      nall = RandomB(hh)%nall
      allocate(RandomB(hh)%nodes(nall,parms%dim_eta))
      allocate(RandomB(hh)%weights(RandomB(hh)%nall))
      call DefineIntegrationNodes(parms%dim_eta,RandomB(hh)%flag, &
                                  nquad, &
                                  RandomB(hh)%nall, &
                                  RandomB(hh)%nodes,RandomB(hh)%weights)
    end do
  end if !if (MC_EFlag) then
end subroutine DefineIntegrationRule

subroutine DefineIntegrationNodes(d,flag,n1,nAll,nodes,weights)
  use ToolsModule, only : kron1
  use GaussianQuadrature, only : GaussianQuadratureRule, hermite, legendre
  use nag_library, only : G05KFF,G05SAF,G05RZF
  implicit none
  ! 0 : Gauss-Hermite               1 : Working
  ! 1 : Sparse-Hermite              2 : NOT working
  ! 2 : pseudo-Montecarlo           3 : working
  ! 3 : quasi-montecarlo            4 : NOT working
  ! 6 : Gauss-Legendre on (-1,1)    6 : working
  ! 7 : pseudo-MC on (-1,1)         7 : working
  integer(i4b), intent(in)  :: d      ! dimension of integration rule
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
  real(dp)                  :: mu(d)
  real(dp), allocatable     :: C(:,:),R(:)

  nodes   = 0.0d0
  weights = 1.0d0
  if (flag==0) then
    ! 0 = Gauss-Hermite
    do i1=1,d
      allocate(x(n1(i1)))
      allocate(w(n1(i1)))
      call GaussianQuadratureRule(hermite, x, w)
      if (i1==1) then
        allocate(e1(nAll/n1(i1)))
        e1 = 1.0d0
        call kron1(x,e1,nodes(:,i1))
        call kron1(w,e1,weights)
        deallocate(e1)
     else if (i1>1 .and. i1<d) then
        allocate(e1(product(n1(1:i1-1))),e2(product(n1(i1+1:d))))
        allocate(temp1(product(n1(i1:d))))
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
      else if (i1==d) then
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
  weights=weights*pi_d**dble(-d/2)
elseif (flag==1) then
!  % 1 = Sparse_hermite
!  SparseLevel = 2  ;
!  nQuadAll = sparse_grid_herm_size ( RCDim, SparseLevel );
!  [weights,nodes]=sparse_grid_herm(RCDim,SparseLevel,nQuadAll);
!  nodes = sqrt(2)*nodes';
!  weights = (pi^(-RCDim/2))*weights';
elseif (flag==2) then
  ! pseudo-MonteCarlo
  ! generate random numbers
  mode  = 2
  ifail = 0
  ldc   = d
  ldx   = nall
  mu    = 0.0d0
  allocate(C(LDC,d))
  C     = 0.0d0
  do i1=1,d
    C(i1,i1) =1.0d0
  end do
  LR = d*(d+1)+1
  allocate(R(LR))
  call G05RZF(mode,nall,d,MU,C,d,R,LR,G05State, &
              nodes,nAll,ifail)
  weights = 1.0d0/dble(nAll)
  deallocate(C,R)
elseif (flag==6) then
  ! Gauss-Legendre on (-1,1)
    do i1=1,d
      allocate(x(n1(i1)))
      allocate(w(n1(i1)))
      call GaussianQuadratureRule(legendre, x, w)
      if (i1==1) then
        allocate(e1(nAll/n1(i1)))
        e1 = 1.0d0
        call kron1(x,e1,nodes(:,i1))
        call kron1(w,e1,weights)
        deallocate(e1)
     else if (i1>1 .and. i1<d) then
        allocate(e1(product(n1(1:i1-1))),e2(product(n1(i1+1:d))))
        allocate(temp1(product(n1(i1:d))))
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
      else if (i1==d) then
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
elseif (flag==7) then
  ! 7 : pseudo-MC on (-1,1)         7 : working

  ! generate random numbers
  ifail = 0
  allocate(x(nall*d))
  call G05SAF(nall*d,G05State,x,ifail)
  x = 2.0d0*x-1.0d0
  nodes = reshape(x,(/nall,d/))
  weights = 1.0d0/dble(nall)
  deallocate(x)
end if

end subroutine DefineIntegrationNodes

#if USE_MPI==1
subroutine BroadcastParameters(pid)
  use mpi
  !use IFPORT
  implicit none
  integer(i4b), intent(in)  :: pid
  integer(i4b)              :: ierr(30),i1,i2,nerr,ierr_barrier
  integer(i4b), allocatable :: ierr_quad(:),ierr_B(:)
  integer(i4b)              :: n  ! size of matrix for RandomB and RandomD
  character(len=20)         :: fmt1
  ! Broadcast OutDir and Control Flags
  ierr = 0
  !call sleep(10)
  call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
  call mpi_bcast(OutDir,len(OutDir),MPI_CHARACTER,MasterID,MPI_COMM_WORLD,ierr(1))
  call mpi_bcast(ControlOptions%OutputFlag,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(2))
  call mpi_bcast(ControlOptions%TestLikeFlag,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(3))
  call mpi_bcast(ControlOptions%TestIntegrationFlag,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(4))
  call mpi_bcast(ControlOptions%SaveDataFlag,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(5))
  call mpi_bcast(ControlOptions%BICFlag,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(6))
  call mpi_bcast(ControlOptions%MPIFlag,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(7))
  call mpi_bcast(MaxOptions%Algorithm,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(8))
  call mpi_bcast(MaxOptions%AbsTol,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(9))
  call mpi_bcast(MaxOptions%RelTol,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(10))
  call mpi_bcast(MaxOptions%MaxLevel,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(11))

  call mpi_bcast(small,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(12))
  call mpi_bcast(inf,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(13))

  ! broadcast dimensions of problem size
  call mpi_bcast(parms%J,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(14))
  call mpi_bcast(parms%K,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(15))
  call mpi_bcast(parms%model,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(16))
  if (parms%model>=2) then
    call mpi_bcast(parms%BC_z_dim,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(17))
    call mpi_bcast(parms%BD_z_dim,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(18))
    call mpi_bcast(parms%dim_eta,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(19))
    call mpi_bcast(parms%BC_lo,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(20))
    call mpi_bcast(parms%BC_hi,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%BC_beta_lo,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(20))
    call mpi_bcast(parms%BC_beta_hi,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%BC_CDiag_lo,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(20))
    call mpi_bcast(parms%BC_CDiag_hi,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%BD_beta_lo,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(20))
    call mpi_bcast(parms%BD_beta_hi,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%BD_CDiag_lo,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(20))
    call mpi_bcast(parms%BD_CDiag_hi,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%InvCDiag_LO,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%InvCDiag_HI,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%InvCOffDiag_LO,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%InvCOffDiag_HI,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%BD_month_hi,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%BD_month_lo,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
  end if

  ! parameters defining counterfactuals for analysis of results
  call mpi_bcast(parms%nPrices_plot,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(21))

  ! data information
  call mpi_bcast(HHData%NMC,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(22))
  call mpi_bcast(HHData%N,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(23))
  call mpi_bcast(HHData%NSim,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(23))
  call mpi_bcast(HHData%M,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(24))

  ! allocate memory for (b,CDiag,COffDiag,xData,iXData)
  !  (this has already been done by pid==MasterID)
  if (pid>MasterID) then
    call AllocateGlobalVariables
  end if

  print 1568,'ierr',pid,ierr(1:24)
1568 format(a4,i4,24i3)

  ! broadcast the other parameters in parms
  call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
  call mpi_bcast(parms%file,len(parms%file),MPI_CHARACTER,MASTERID,MPI_COMM_WORLD,ierr(24))
  call mpi_bcast(parms%unit,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(24))
  call BroadcastParms(parms,pid)

  ! broadcast information on integration rule
  if (pid>masterID) then
    allocate(RandomE(parms%K))
  end if
  call mpi_bcast(RandomE%flag,parms%K,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(25))
  call mpi_bcast(RandomE%nAll,parms%K,MPI_INteger,MasterID,MPI_COMM_WORLD,ierr(26))
  call mpi_bcast(RandomE%nquadRaw,parms%K,MPI_INteger,MasterID,MPI_COMM_WORLD,ierr(26))

  if (pid>MasterID) allocate(RandomB(1))
  if (pid>MasterID) allocate(RandomB(1)%nquad(parms%dim_eta))
  call mpi_bcast(RandomB%flag,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(25))
  call mpi_bcast(RandomB%nAll,1,MPI_INteger,MasterID,MPI_COMM_WORLD,ierr(26))
  call mpi_bcast(RandomB(1)%nQuad,parms%dim_eta,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr(26))

  if (pid==MasterID) n = size(G05State)
  call mpi_bcast(n,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr(26))
  if (pid>MasterID) allocate(G05State(n))
  call mpi_bcast(G05State,n,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(26))

  call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)

end subroutine BroadcastParameters

subroutine BroadcastParms(LocalParms,pid)
  use mpi
  implicit none
  type(ParmsStructure), intent(inout) :: LocalParms
  integer(i4b), intent(in)            :: pid
  integer(i4b)                        :: i1,nerr,ierr_barrier
  integer(i4b), allocatable           :: ierr(:)
  character(len=20)                   :: fmt1

  allocate(ierr(2*LocalParms%J+3*LocalParms%K+12+2*LocalParms%dim_eta &
                +LocalParms%BD_Z_DIM+LocalParms%BC_Z_DIM+11))
  ierr = 0
  ! broadcast b,CDiag,COffDiag
  call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)

  call mpi_bcast(LocalParms%D,LocalParms%J,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(1))
  call mpi_bcast(LocalParms%BC,LocalParms%nBC,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(2))

  do i1=1,LocalParms%J
    call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
    call mpi_bcast(LocalParms%B(:,i1),LocalParms%K,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(2+i1))
  end do
  nerr = 2+LocalParms%J

  call mpi_bcast(LocalParms%mue,LocalParms%K,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+1))
  call mpi_bcast(LocalParms%InvCDiag,LocalParms%K,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+2))
  call mpi_bcast(LocalParms%InvCOffDiag,LocalParms%K*(LocalParms%K-1)/2,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+3))
  nerr = nerr+3

  do i1=1,12
    call mpi_bcast(LocalParms%mue_month(:,i1),LocalParms%K,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+1))
  end do
  nerr=nerr+12

  do i1=1,LocalParms%k
    call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
    call mpi_bcast(LocalParms%InvC(:,i1),LocalParms%K,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+i1))
    call mpi_bcast(LocalParms%CSig(:,i1),LocalParms%K,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+Localparms%K+i1))
    call mpi_bcast(LocalParms%Sig(:,i1),LocalParms%K,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+2*Localparms%K+i1))
  end do
  nerr = nerr+3*LocalParms%K

  if (LocalParms%model>=2) then
    call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
    call mpi_bcast(LocalParms%BD_beta,LocalParms%BD_z_dim,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+1))
    call mpi_bcast(LocalParms%BD_CDiag,LocalParms%J,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+2))
    call mpi_bcast(LocalParms%BD_COffDiag,LocalParms%nBD_COffDiag,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+3))
    do i1=1,LocalParms%J
      call mpi_bcast(LocalParms%BD_month(i1,:),12,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+3))
    end do
    nerr = nerr+3+LocalParms%J
    do i1=1,LocalParms%dim_eta
      call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
      call mpi_bcast(LocalParms%BD_C(:,i1),LocalParms%J,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+i1))
    end do
    nerr = nerr+LocalParms%dim_eta

    call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
    call mpi_bcast(LocalParms%BC_beta,LocalParms%BC_z_dim,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+1))
    call mpi_bcast(LocalParms%BC_CDiag,LocalParms%nBC,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+2))
    call mpi_bcast(LocalParms%BC_COffDiag,LocalParms%nBC_COffDiag,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+3))
    do i1=1,LocalParms%dim_eta
      call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
      call mpi_bcast(LocalParms%BC_C(:,i1),LocalParms%nBC,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+3+i1))
    end do
    nerr = nerr + 3 + LocalParms%dim_eta
    do i1=1,LocalParms%BD_Z_dim
      call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
      call mpi_bcast(LocalParms%BD_Z(:,i1),LocalParms%J,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+i1))
    end do
    nerr = nerr + LocalParms%BD_Z_dim
    do i1=1,LocalParms%BC_Z_dim
      call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
      call mpi_bcast(LocalParms%BC_Z(:,i1),LocalParms%nBC,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+i1))
    end do
    nerr = nerr + LocalParms%BC_Z_DIM
  end if

  write(fmt1,'(a8,i3,a3)') '(a11,i4,',nerr,'i3)'
  print fmt1,'parms_ierr',pid,ierr
!1677 format(a11,i4,<nerr>i3)
end subroutine BroadcastParms

subroutine BroadcastIFree(pid)
  use mpi
  implicit none
  integer(i4b), intent(in) :: pid
  integer(i4b) :: ierr
  integer(i4b), allocatable :: temp1(:),temp2(:)

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagD,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagBC,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagMUE,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagMUE_month,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagInvCDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagInvCOffDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)

  call mpi_bcast(iFree%flagBC_beta,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagBC_CDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagBC_COffDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagBD_beta,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagBD_CDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagBD_COffDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagBD_month,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nPerIter,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%RandFlag,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%seed,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)

  call mpi_bcast(iFree%MUE1,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%MUE_month1,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%InvCDiag1,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%InvCOffDiag1,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%BC_beta1,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%BC_CDiag1,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%BC_COffDiag1,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%BD_beta1,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%BD_CDiag1,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%BD_COffDiag1,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%BD_month1,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)

  call mpi_bcast(iFree%nD,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nBC,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nInvCDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nInvCOffDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nMuE,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nMuE_month,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)

  call mpi_bcast(iFree%nBD_beta,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nBD_CDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nBD_COffDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nBD_month,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)

  call mpi_bcast(iFree%nBC_beta,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nBC_CDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nBC_COffDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)

  call mpi_bcast(iFree%nAll,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  print *,"pid = ",pid," broadcast of size iFree complete."
  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if (iFree%nD>0) then
    allocate(temp1(iFree%nD))
    allocate(temp2(iFree%nD))
    if (pid>MasterID .and. .not. allocated(iFree%d)) then
      allocate(iFree%D(iFree%nD))
      allocate(iFree%xD(iFree%nD))
    end if
    if (pid==MasterID) then
      temp1 = iFree%D
      temp2 = iFree%xD
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp1,iFree%nD,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp2,iFree%nD,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    iFree%D = temp1
    iFree%xD = temp2
    deallocate(temp1,temp2)
  end if

  if (iFree%nBC>0) then
    allocate(temp1(iFree%nBC))
    allocate(temp2(iFree%nBC))
    if (pid>MasterID .and. .not. allocated(iFree%BC)) then
      allocate(iFree%BC(iFree%nBC))
      allocate(iFree%xBC(iFree%nBC))
    end if
    if (pid==MasterID) then
      temp1 = iFree%BC
      temp2 = iFree%xBC
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp1,iFree%nBC,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp2,iFree%nBC,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    iFree%BC = temp1
    iFree%xBC = temp2
    deallocate(temp1,temp2)
  end if

  if (iFree%nMuE>0) then
    allocate(temp1(iFree%nMUE))
    allocate(temp2(iFree%nMUE))
    if (pid>MasterID .and. .not. allocated(iFree%MUE)) then
      allocate(iFree%MuE(iFree%nMuE))
      allocate(iFree%xMuE(iFree%nMuE))
    end if
    if (pid==MasterID) then
      temp1 = iFree%MUE
      temp2 = iFree%xMUE
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp1,iFree%nMuE,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp2,iFree%nMuE,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    iFree%MUE = temp1
    iFree%xMUE = temp2
    deallocate(temp1,temp2)
  end if

  if (iFree%nMuE_month>0) then
    allocate(temp1(iFree%nMUE_month))
    allocate(temp2(iFree%nMUE_month))
    if (pid>MasterID .and. .not. allocated(iFree%MUE_month)) then
      allocate(iFree%MuE_month(iFree%nMuE_month))
      allocate(iFree%xMuE_month(iFree%nMuE_month))
    end if
    if (pid==MasterID) then
      temp1 = iFree%MUE_month
      temp2 = iFree%xMUE_month
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp1,iFree%nMuE_month,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp2,iFree%nMuE_month,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    iFree%MUE_month = temp1
    iFree%xMUE_month = temp2
    deallocate(temp1,temp2)
  end if

  if (iFree%nInvCDiag>0) then
    allocate(temp1(iFree%nInvCDiag))
    allocate(temp2(iFree%nInvCDiag))
    if (pid>MasterID .and. .not. allocated(iFree%InvCDiag)) then
      allocate(iFree%InvCDiag(iFree%nInvCDiag))
      allocate(iFree%xInvCDiag(iFree%nInvCDiag))
    end if
    if (pid==MasterID) then
      temp1 = iFree%InvCDiag
      temp2 = iFree%xInvCDiag
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp1,iFree%nInvCDiag,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp2,iFree%nInvCDiag,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    iFree%InvCDiag = temp1
    iFree%xInvCDiag = temp2
    deallocate(temp1,temp2)
  end if

  if (iFree%nInvCOffDiag>0) then
    allocate(temp1(iFree%nInvCOffDiag))
    allocate(temp2(iFree%nInvCOffDiag))

    if (pid>MasterID .and. .not. allocated(iFree%InvCOffDiag)) then
      allocate(iFree%InvCOffDiag(iFree%nInvCOffDiag))
      allocate(iFree%xInvCOffDiag(iFree%nInvCOffDiag))
    end if
    if (pid==MasterID) then
      temp1 = iFree%InvCOffDiag
      temp2 = iFree%xInvCOffDiag
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp1,iFree%nInvCOffDiag,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp2,iFree%nInvCOffDiag,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    iFree%InvCOffDiag = temp1
    iFree%xInvCOffDiag = temp2
    deallocate(temp1,temp2)
  end if

  if (iFree%nBD_Beta>0) then
    allocate(temp1(iFree%nBD_beta))
    allocate(temp2(iFree%nBD_beta))

    if (pid>MasterID .and. .not. allocated(iFree%BD_beta)) then
      allocate(iFree%BD_beta(iFree%nBD_beta))
      allocate(iFree%xBD_beta(iFree%nBD_beta))
    end if
    if (pid==MasterID) then
      temp1 = iFree%BD_beta
      temp2 = iFree%xBD_beta
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp1,iFree%nBD_beta,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp2,iFree%nBD_beta,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    iFree%BD_beta = temp1
    iFree%xBD_beta = temp2
    deallocate(temp1,temp2)
  end if

  if (iFree%nBD_month>0) then
    allocate(temp1(iFree%nBD_month))
    allocate(temp2(iFree%nBD_month))

    if (pid>MasterID .and. .not. allocated(iFree%BD_month)) then
      allocate(iFree%BD_month(iFree%nBD_month))
      allocate(iFree%xBD_month(iFree%nBD_month))
    end if
    if (pid==MasterID) then
      temp1 = iFree%BD_month
      temp2 = iFree%xBD_month
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp1,iFree%nBD_month,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp2,iFree%nBD_month,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    iFree%BD_month = temp1
    iFree%xBD_month = temp2
    deallocate(temp1,temp2)
  end if

  if (iFree%nBD_CDiag>0) then
    allocate(temp1(iFree%nBD_CDiag))
    allocate(temp2(iFree%nBD_CDiag))

    if (pid>MasterID .and. .not. allocated(iFree%BD_CDiag)) then
      allocate(iFree%BD_CDiag(iFree%nBD_CDiag))
      allocate(iFree%xBD_CDiag(iFree%nBD_CDiag))
    end if
    if (pid==MasterID) then
      temp1 = iFree%BD_CDiag
      temp2 = iFree%xBD_CDiag
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp1,iFree%nBD_CDiag,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp2,iFree%nBD_CDiag,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    iFree%BD_CDiag = temp1
    iFree%xBD_CDiag = temp2
    deallocate(temp1,temp2)
  end if

  if (iFree%nBD_COffDiag>0) then
    allocate(temp1(iFree%nBD_COffDiag))
    allocate(temp2(iFree%nBD_COffDiag))
    if (pid>MasterID .and. .not. allocated(iFree%BD_COffDiag)) then
      allocate(iFree%BD_COffDiag(iFree%nBD_COffDiag))
      allocate(iFree%xBD_COffDiag(iFree%nBD_COffDiag))
    end if
    if (pid==MasterID) then
      temp1 = iFree%BD_COffDiag
      temp2 = iFree%xBD_COffDiag
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp1,iFree%nBD_COffDiag,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp2,iFree%nBD_COffDiag,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    iFree%BD_COffDiag = temp1
    iFree%xBD_COffDiag = temp2
    deallocate(temp1,temp2)
  end if

  if (iFree%nBC_Beta>0) then
    allocate(temp1(iFree%nBC_Beta))
    allocate(temp2(iFree%nBC_Beta))
    if (pid>MasterID .and. .not. allocated(iFree%BC_beta)) then
      allocate(iFree%BC_beta(iFree%nBC_beta))
      allocate(iFree%xBC_beta(iFree%nBC_beta))
    end if
    if (pid==MasterID) then
      temp1 = iFree%BC_beta
      temp2 = iFree%xBC_beta
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp1,iFree%nBC_beta,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp2,iFree%nBC_beta,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    iFree%BC_beta = temp1
    iFree%xBC_beta = temp2
    deallocate(temp1,temp2)
  end if

  if (iFree%nBC_CDiag>0) then
    allocate(temp1(iFree%nBC_CDiag))
    allocate(temp2(iFree%nBC_CDiag))
    if (pid>MasterID .and. .not. allocated(iFree%BC_CDiag)) then
      allocate(iFree%BC_CDiag(iFree%nBC_CDiag))
      allocate(iFree%xBC_CDiag(iFree%nBC_CDiag))
    end if
    if (pid==MasterID) then
      temp1 = iFree%BC_CDiag
      temp2 = iFree%xBC_CDiag
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp1,iFree%nBC_CDiag,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp2,iFree%nBC_CDiag,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    iFree%BC_CDiag = temp1
    iFree%xBC_CDiag = temp2
    deallocate(temp1,temp2)
  end if

  if (iFree%nBC_COffDiag>0) then
    allocate(temp1(iFree%nBC_COffDiag))
    allocate(temp2(iFree%nBC_COffDiag))
    if (pid>MasterID .and. .not. allocated(iFree%BC_COffDiag)) then
      allocate(iFree%BC_COffDiag(iFree%nBC_COffDiag))
      allocate(iFree%xBC_COffDiag(iFree%nBC_COffDiag))
    end if
    if (pid==MasterID) then
      temp1 = iFree%BC_COffDiag
      temp2 = iFree%xBC_COffDiag
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp1,iFree%nBC_COffDiag,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    call mpi_bcast(temp2,iFree%nBC_COffDiag,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
    iFree%BC_COffDiag = temp1
    iFree%xBC_COffDiag = temp2
    deallocate(temp1,temp2)
  end if

end subroutine BroadcastIFree

#endif

! Written :  2015AUG14 LN
subroutine ReadWriteParameters(LocalParms,LocalAction)
  use NewTools, only : MatrixToSphere,MatrixInverse
#ifndef  __GFORTRAN__
  use IFPORT       ! intel fortran portability library
#endif
  implicit none
  type(ParmsStructure), intent(inout) :: LocalParms
  character(LEN=*),     intent(in)    :: LocalAction

  integer(i4b)                        :: i1,FileStatus
  character(len=30)                   :: TempString
  character(len=30)                   :: fmt102,fmt103,fmt104,fmt107,fmt108,fmt109
  character(len=8)                    :: datestring
  character(len=200)                  :: parms_copy_file
  character(len=400)                  :: cmd_string
  ! formats for strings in parms.csv
  111 format(a13)   !'model,K,J,nBC'
  112 format(a9)    ! B(K x J)
  113 format(a11)   ! MUE(K x 1)
  114 format(a12)   ! CSIG (K x K)
  115 format(a25)   ! dim_eta,BC_Z_DIM,BD_Z_DIM
  116 format(a37)   ! nBC_COffDiag,nBD_COffDiag,BC_lo,BC_hi
  117 format(a21)   ! BC_Beta(BC_Z_Dim x 1)
  118 format(a21)   ! BD_beta(BD_Z_DIM x 1)
  119 format(a19)   ! BC_C(nBC x dim_eta)

  ! formats for data in parms.csv
  101 format(3(i4,','),i4)                               ! (model,K,J,nBC)
  write(fmt102,'(a1,i2,a15)') '(',LocalParms%J,'(f25.16,:,","))'
  write(fmt103,'(a1,i2,a15)') '(',LocalParms%K,'(f25.16,:,","))'
  write(fmt104,'(a1,i2,a15)') '(',LocalParms%K,'(f25.16,:,","))'
  105 format(2(i4,','),i4)                               ! (dim_eta,BC_Z_DIM,BD_Z_DIM)
  106 format(2(i4,','),f25.16,',',f25.16)                ! nBC_COffDiag,nBD_COffDiag,BC_LO,BC_HI
  write(fmt107,'(a1,i2,a15)') '(',LocalParms%BC_Z_DIM,'(f25.16,:,","))'
  write(fmt108,'(a1,i2,a15)') '(',LocalParms%BD_Z_DIM,'(f25.16,:,","))'
  write(fmt109,'(a1,i2,a15)') '(',LocalParms%DIM_ETA,'(f25.16,:,","))'
  110 format(12(f25.16,:,','))                           ! BD_month

  120 format(a17)      ! mue_month(K x 12)
  121 format(12(f25.16,:,','))

  if (LocalAction=='read') then
    call date_and_time(date=datestring)
    parms_copy_file = trim(OutDir) // "/parms/parms" // datestring // ".csv"
    cmd_string = "cp " // LocalParms%file // " " // parms_copy_file
    FileStatus = system(cmd_string)
  end if

  ! open file for reading or writing
  open(unit=LocalParms%unit,  &
       file=LocalParms%file,  &
       action=LocalAction)

  ! choose action
  select case (LocalAction)

  !-------------------------------------------------------------------------------------------
  ! Read parameters from file
  !-------------------------------------------------------------------------------------------
  case ('read')

    ! (model,K,J,nBC)
    read(LocalParms%unit,111) TempString
    print *,TempString
    read(LocalParms%unit,101) LocalParms%model,LocalParms%K,LocalParms%J,LocalParms%nBC

    ! B (K x J)
    read(LocalParms%unit,112) TempString
    print *,TempString
    do i1=1,LocalParms%K
      read(LocalParms%unit,fmt102) LocalParms%B(i1,:)
    end do

    ! update (D,BC): convert B to spherical coordinates (D,BC)
    call MatrixToSphere(LocalParms%B,LocalParms%D,LocalParms%BC)

    ! MUE   (K x 1)
    read(LocalParms%unit,113) TempString
    print *, TempString
    read(LocalParms%unit,fmt103) LocalParms%MUE

    ! MUE_month   (K x 12)
    read(LocalParms%unit,120) TempString
    print *, TempString
    do i1=1,LocalParms%K
      read(LocalParms%unit,121) LocalParms%MUE_month(i1,:)
    end do

    ! CSIG  (K x K)  Cholesky decomposition of sig
    read(LocalParms%unit,114) TempString
    print *,TempString
    do i1=1,LocalParms%K
      read(LocalParms%unit,fmt104) LocalParms%CSIG(i1,:)
    end do

    ! SIG = covariance matrix of epsilon
    ! InvC = inv(CSIG)
    ! (InvCDiag,InvCOffDiag) = spherical representation of tranpose(InvC)
    LocalParms%SIG = matmul(LocalParms%CSIG,transpose(LocalParms%CSIG))
    call MatrixInverse(LocalParms%CSIG,LocalParms%InvC,'Lower triangular')
    call MatrixToSphere(transpose(LocalParms%InvC),LocalParms%InvCDiag,LocalParms%InvCOffDiag)

    ! (dim_eta,BC_Z_DIM,BD_Z_DIM)
    read(LocalParms%unit,115) TempString
    print *,TempString
    read(LocalParms%unit,105) LocalParms%dim_eta,LocalParms%BC_Z_DIM,LocalParms%BD_Z_DIM

    ! (nBC_COffDiag,nBD_COffDiag)
    read(LocalParms%unit,116) TempString
    print *, TempString
    read(LocalParms%unit,106) LocalParms%nBC_COffDiag,LocalParms%nBD_COffDiag, &
                                LocalParms%BC_lo,LocalParms%BC_hi
    ! BC_beta  (BC_Z_DIM)
    read(LocalParms%unit,117) TempString
    print *, TempString
    read(LocalParms%unit,fmt107) LocalParms%BC_BETA

    ! BD_beta  (BD_Z_DIM)
    read(LocalParms%unit,118) TempString
    print *, TempString
    read(LocalParms%unit,fmt108) LocalParms%BD_BETA

    ! BC_C  (nBC x dim_eta)
    read(LocalParms%unit,119) TempString
    print *,TempString
    do i1=1,LocalParms%nBC
      read(LocalParms%unit,fmt109) LocalParms%BC_C(i1,:)
    end do

    ! BD_C  (J x dim_eta)
    read(LocalParms%unit,119) TempString
    print *,TempString
    do i1=1,LocalParms%J
      read(LocalParms%unit,fmt109) LocalParms%BD_C(i1,:)
    end do

   ! BD_month (J x 12)
   read(LocalParms%unit,119) TempString
   print *,TempString
   do i1=1,localparms%J
     read(LocalParms%unit,110) LocalParms%BD_month(i1,:)
   end do

    ! compute spherical coordinate representation of (BC_C,BD_C)
    ! BC_C = lower triangular. MatrixToSphere requires upper triangular
    ! BC_D = lower triangular. MatrixToSphere requires upper triangular
    call MatrixToSphere(transpose(LocalParms%BC_C),LocalParms%BC_CDiag,LocalParms%BC_COffDiag)
    call MatrixToSphere(transpose(LocalParms%BD_C),LocalParms%BD_CDiag,LocalParms%BD_COffDiag)

  !-------------------------------------------------------------------------------------------
  ! Write parameters to file
  ! -------------------------------------------------------------------------------------------
  case ('write')

    ! (model,K,J,nBC)
    write(LocalParms%unit,111) 'model,K,J,nBC'
    write(LocalParms%unit,101) LocalParms%model,LocalParms%K,LocalParms%J,LocalParms%nBC

    ! B  (K x J)
    write(LocalParms%unit,112) 'B(K x J)'
    do i1=1,LocalParms%K
      write(LocalParms%unit,fmt102) LocalParms%B(i1,:)
    end do

    ! MUE
    write(LocalParms%unit,113) 'MUE(K x 1)'
    write(LocalParms%unit,fmt103) LocalParms%MUE

    write(LocalParms%unit,120) 'mue_month(K x 12)'
    do i1=1,LocalParms%K
      write(LocalParms%unit,121) LocalParms%mue_month(i1,:)
    end do

    ! CSIG  (K x K)  Cholesky decomposition of sig
    write(LocalParms%unit,114) 'CSIG(K x K)'
    do i1=1,LocalParms%K
      write(LocalParms%unit,fmt104) LocalParms%CSIG(i1,:)
    end do

    ! (dim_eta,BC_Z_DIM,BD_Z_DIM)
    write(LocalParms%unit,115) 'dim_eta,BC_Z_DIM,BD_Z_DIM'
    write(LocalParms%unit,105) LocalParms%dim_eta,LocalParms%BC_Z_DIM,LocalParms%BD_Z_DIM

    ! (nBC_COffDiag,nBD_COffDiag,BC_lo,BC_hi)
    write(LocalParms%unit,116) 'nBC_COffDiag,nBD_COffDiag,BC_lo,BC_hi'
    write(LocalParms%unit,106) LocalParms%nBC_COffDiag,LocalParms%nBD_COffDiag, &
                                    LocalParms%BC_lo,LocalParms%BC_hi

    ! BC_beta  (BC_Z_DIM)
    write(LocalParms%unit,117) 'BC_BETA(BC_Z_DIM x 1)'
    write(LocalParms%unit,fmt107) LocalParms%BC_BETA

    ! BD_beta  (BD_Z_DIM)
    write(LocalParms%unit,118) 'BD_BETA(BD_Z_DIM x 1)'
    write(LocalParms%unit,fmt108) LocalParms%BD_BETA

    ! BC_C  (nBC x dim_eta)
    write(LocalParms%unit,119) 'BC_C(nBC x dim_eta)'
    do i1=1,LocalParms%nBC
      write(LocalParms%unit,fmt109) LocalParms%BC_C(i1,:)
    end do

    ! BD_C  (J x dim_eta)
    write(LocalParms%unit,119) 'BD_C(J x dim_eta)'
    do i1=1,LocalParms%J
      write(LocalParms%unit,fmt109) LocalParms%BD_C(i1,:)
    end do

   ! BD_month (J x 12)
   write(LocalParms%unit,119) 'BD_month(J x 12)'
   do i1=1,localparms%J
     write(LocalParms%unit,110) LocalParms%BD_month(i1,:)
   end do
  end select

  close(LocalParms%unit)

end subroutine ReadWriteParameters

subroutine CopyParameters(parms_in,parms_out)
  implicit none
  type(ParmsStructure), intent(in)    :: parms_in
  type(ParmsStructure), intent(inout) :: parms_out

  parms_out%file         = parms_in%file
  parms_out%unit         = parms_in%unit
  parms_out%model        = parms_in%model
  parms_out%K            = parms_in%K
  parms_out%J            = parms_in%J
  parms_out%nBC          = parms_in%nBC

  parms_out%dim_eta      = parms_in%dim_eta
  parms_out%BD_z_dim     = parms_in%BD_z_dim
  parms_out%BC_z_dim     = parms_in%BC_z_dim
  parms_out%nBC_COffDiag = parms_in%nBC_COffDiag
  parms_out%nBD_COffDiag = parms_in%nBD_COffDiag

  ! Allocate memory for parms_out
  call DeallocateParms(parms_out)
  call AllocateParms(parms_out)

  parms_out%B = parms_in%B
  parms_out%D = parms_in%D
  parms_out%BC = parms_in%BC
  parms_out%MUE = parms_in%MUE
  parms_out%MUE_month = parms_in%MUE_month
  parms_out%sig = parms_in%sig
  parms_out%CSIG = parms_in%CSIG
  parms_out%invc = parms_in%invc
  parms_out%InvCDiag = parms_in%InvCDiag
  parms_out%InvCOffDiag = parms_in%InvCOffDiag

  if (parms_out%model>=2) then
    parms_out%BC_beta     = parms_in%BC_beta
    parms_out%BC_C        = parms_in%BC_C
    parms_out%BC_z        = parms_in%BC_z
    parms_out%BC_CDiag    = parms_in%BC_CDiag
    parms_out%BC_COffDiag = parms_in%BC_COffDiag
    parms_out%GradBC_C_CDiag    = parms_in%GradBC_C_CDiag
    parms_out%GradBC_C_COffDiag = parms_in%GradBC_C_COffDiag

    parms_out%BD_beta     = parms_in%BD_beta
    parms_out%BD_C        = parms_in%BD_C
    parms_out%BD_z        = parms_in%BD_z
    parms_out%BD_CDiag    = parms_in%BD_CDiag
    parms_out%BD_COffDiag = parms_in%BD_COffDiag
    parms_out%BD_month    = parms_in%BD_month
    parms_out%GradBD_C_CDiag    = parms_in%GradBD_C_CDiag
    parms_out%GradBD_C_COffDiag = parms_in%GradBD_C_COffDiag

    parms_out%sigp = parms_in%sigp

  end if ! if (parms_out%model>=2) then

  parms_out%BC_lo        = parms_in%BC_lo
  parms_out%BC_hi        = parms_in%BC_hi
  parms_out%InvCDiag_lo  = parms_in%InvCDiag_lo
  parms_out%InvCDiag_hi  = parms_in%InvCDiag_hi
  parms_out%InvCOffDiag_lo  = parms_in%InvCOffDiag_lo
  parms_out%InvCOffDiag_hi  = parms_in%InvCOffDiag_hi
  parms_out%BC_beta_lo  = parms_in%BC_beta_lo
  parms_out%BC_beta_hi  = parms_in%BC_beta_hi
  parms_out%BC_CDiag_lo  = parms_in%BC_CDiag_lo
  parms_out%BC_CDiag_hi  = parms_in%BC_CDiag_hi
  parms_out%BD_beta_lo  = parms_in%BD_beta_lo
  parms_out%BD_beta_hi  = parms_in%BD_beta_hi
  parms_out%BD_CDiag_lo  = parms_in%BD_CDiag_lo
  parms_out%BD_CDiag_hi  = parms_in%BD_CDiag_hi
  parms_out%BD_month_lo = parms_in%BD_month_lo
  parms_out%BD_month_hi = parms_in%BD_month_hi

end subroutine CopyParameters

subroutine LoadBasePrice(p0)
  implicit none
  real(dp), intent(out) :: p0(:)
  integer(i4b) :: i1
  integer(i4b) :: BasePriceUnit

  BasePriceUnit = 100

  open(unit = BasePriceUnit, &
       file = ParmFiles%BasePriceFile, &
       action = 'read')

  do i1=1,parms%J
    read(BasePriceUnit,'(f25.0)') p0(i1)
  end do

  close(BasePriceUnit)
end subroutine LoadBasePrice

! Load tax parameters from file
! taxid    = id number for tax experiment
! taxlabel = description of tax experiment
! taxtype  = "ad valorem" or "excise"
! tax      = (J x ntax) for each tax vector of tax rates
subroutine LoadTaxParameters(taxid,taxlabel,taxtype,tax)
  implicit none
  integer(i4b),       allocatable,intent(out) :: taxid(:)
  character(len=100), allocatable,intent(out) :: taxlabel(:)
  character(len=10),  allocatable,intent(out) :: taxtype(:)
  real(dp),           allocatable,intent(out) :: tax(:,:)

  integer(i4b)       :: TaxParmsUnit

  character(len=1024)       :: buffer1,buffer2
  integer(i4b), allocatable :: itemp(:)
  integer(i4b)              :: i1,icomma,ntax
  character(len=20)         :: fmt1
  ! tax parameter file: (J+3 x ntax)  csv file
  ! row 1       taxid
  ! row 2       taxlabel
  ! row 3       tax type
  ! row 4 - J+3 tax rate on product 1 - J
  TaxParmsUnit = 45
  open(unit = TaxParmsUnit, &
       file = ParmFiles%TaxParmsFile, &
       action = "read")

  ! read in row 1: taxid (integers) and determine number of columns
  read(TaxParmsUnit,'(a)') buffer1

  i1 = 1
  icomma = 1
  allocate(itemp(20))
  itemp = 0
  do while (icomma>0)
    icomma = index(buffer1,",")
    if (icomma>0) then
      read(buffer1(1:(icomma-1)),'(i4)') itemp(i1)
      i1=i1+1
      buffer1=buffer1((icomma+1):len_trim(buffer1))
    else
      read(buffer1,'(i4)') itemp(i1)
    end if
  end do
  ntax = i1
  allocate(taxid(ntax))
  taxid = itemp(1:i1)

  ! Read in rows 2-3
  allocate(taxlabel(ntax),taxtype(ntax))

  read(TaxParmsUnit,'(a)') buffer1
  read(TaxParmsUnit,'(a)') buffer2

  do i1=1,ntax
    if (i1<ntax) then
      icomma       = index(buffer1,",")
      taxlabel(i1) = buffer1(1:(icomma-1))
      buffer1      = buffer1((icomma+1):len_trim(buffer1))

      icomma       = index(buffer2,",")
      taxtype(i1)  = buffer2(1:(icomma-1))
      buffer2      = buffer2((icomma+1):len_trim(buffer2))
    else
      taxlabel(i1) = buffer1
      taxtype(i1)  = buffer2
    end if
  end do

  ! read in (J x n1) tax rates
  allocate(tax(parms%J,ntax))
  tax = 1.0d0
  write(fmt1,'(a1,i2,a6)') '(',ntax,'g13.4)'
  do i1=1,parms%J
    read(TaxParmsUnit,fmt1) tax(i1,:)
  end do
! 69 format(<ntax>g13.4)

  deallocate(itemp)

  close(TaxParmsUnit)
end subroutine LoadTaxParameters

! initialize random number generator
subroutine  InitializeG05
  use nag_library, only : G05KFF
  implicit none
  integer(i4b)              :: ifail
  integer(i4b)              :: genid,subid,lseed,lstate
  integer(i4b), allocatable :: seed(:)

  ! initialize random number generator
  genid = 3   ! Mersenne twister algorithm
  subid = 0   ! not used when genid==3
  lseed = 624 ! number of seeds needed for genid==3
  allocate(seed(lseed))
  call SetSeed(seed)

  lstate = 0
  allocate(G05State(lstate))
  ifail=0
  call G05KFF(genid,subid,seed,lseed,G05State,lstate,ifail)
  deallocate(G05State)
  allocate(G05State(lstate))

  ifail = 0
  call G05KFF(genid,subid,seed,lseed,G05State,lstate,ifail)
  deallocate(seed)
end subroutine InitializeG05

! create random sequence of integers to use as seed for random number generator
subroutine SetSeed(seed)
  use nag_library, only : G05KFF,G05TLF
  implicit none
  integer(i4b), intent(out) :: seed(:)
  integer(i4b)              :: n
  integer(i4b)              :: a,b
  integer(i4b)              :: ifail

  integer(i4b) :: genid,subid,lseed,lstate
  integer(i4b), allocatable :: seed0(:),state(:)

  n = size(seed,1)
  genid = 1
  subid = 0
  ifail = 0
  lseed = 1
  allocate(seed0(lseed))
  seed0 = 14620580

  lstate = 0
  allocate(state(lstate))
  ifail = 0
  call G05KFF(genid,subid,seed0,lseed,state,lstate,ifail)
  deallocate(state)
  allocate(state(lstate))
  ifail = 0
  call G05KFF(genid,subid,seed0,lseed,state,lstate,ifail)

  a     = 1
  b     = 9999999
  call G05TLF(n,a,b,state,seed,ifail)

  deallocate(seed0,state)
end subroutine SetSeed

end module GlobalModule
