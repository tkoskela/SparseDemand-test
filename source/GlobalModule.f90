! Revision history
! 2015AUG14  LN  add ReadWriteParameters and CopyParameters

module GlobalModule
  use nrtype
  use PropertyList   ! tools for loading parameters from file
  implicit none

  type(Property)    :: PropList   ! list of keys in input file

  ! structure storing filenames
  type FilenameStructure
    character(len=99) :: mue,C,D
    character(len=99) :: InvCDiag,InvCOffDiag
    character(len=99) :: BC_beta,BC_CDiag,BC_COffDiag
    character(len=99) :: BD_beta,BD_CDiag,BD_COffDiag
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

  type QuadNode
    integer(i4b), allocatable :: nQuad(:)   ! number of nodes in each dim.
    real(dp),     allocatable :: nodes(:,:) ! nAll x dim
    real(dp),     allocatable :: weights(:) ! nAll x 1
    integer(i4b)              :: nall
  end type

  type QuadRuleType
    integer(i4b)                :: dMax
    integer(i4b), allocatable   :: flag(:)  ! one for each dimension
    integer(i4b), allocatable   :: nall(:)
    type(QuadNode), allocatable :: rule(:)
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
    integer(i4b),   allocatable :: internet(:)
    integer(i4b),   allocatable :: SmallStore(:)
    integer(i4b),   allocatable :: date(:),day(:),month(:),week(:)
    integer(i4b)                :: nRawVars    ! number of variables in raw data file
    integer(i4b)                :: EstimationSeed
    integer(i4b)                :: SimulationSeed
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
    real(dp)              :: BD_beta_lo   ! lower bound
    real(dp)              :: BD_beta_hi   ! upper bound
    
    ! parameters used to determine details of analysisof results
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

    integer(i4b)              :: nD 
    integer(i4b)              :: nBC
    integer(i4b)              :: nMuE
    integer(i4b)              :: nInvCDiag
    integer(i4b)              :: nInvCOffDiag

    integer(i4b)              :: nBD_beta
    integer(i4b)              :: nBD_CDiag
    integer(i4b)              :: nBD_COffDiag
    
    integer(i4b)              :: nBC_beta
    integer(i4b)              :: nBC_CDiag
    integer(i4b)              :: nBC_COffDiag

    integer(i4b)              :: nAll

    integer(i4b)              :: flagD
    integer(i4b)              :: flagBC
    integer(i4b)              :: flagMUE
    integer(i4b)              :: flagInvCDiag
    integer(i4b)              :: flagInvCOffDiag

    integer(i4b)              :: flagBC_beta
    integer(i4b)              :: flagBC_CDiag
    integer(i4b)              :: flagBC_COffDiag
    integer(i4b)              :: flagBD_beta
    integer(i4b)              :: flagBD_CDiag
    integer(i4b)              :: flagBD_COffDiag
    integer(i4b)              :: OneAtATime

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
    real(dp)           :: MaxLevel   ! maximum level for D01ESF Bayesian computations
    character(len=200) :: OptionsFile   ! file with E04WDF options
    character(len=200) :: BasisFile
    character(len=200) :: BackupBasisFile
    character(len=200) :: OldBasisFile
    integer(i4b)       :: SaveBasis  ! 1 to save basis info in BasisFile
    integer(i4b)       :: LoadBasis  ! 1 to load basis info from OldBasisFile
    integer(i4b)       :: SaveBasisFreq  ! save basis every i1 iterations 
  end type

  type BayesType
    integer(i4b)    :: nAll
    real(dp), allocatable :: x(:,:)   ! x(:,i1) = draw number i1 from density of x
    real(dp), allocatable :: w(:)     ! integration weights
    real(dp), allocatable :: prior(:) ! prior evaluated at x(:,i1)

    real(dp)              :: BD_beta_lo,BD_beta_hi
    real(dp)              :: BD_CDiag_lo,BD_CDiag_hi
    real(dp)              :: BD_COffDiag_lo,BD_COffDiag_hi

    real(dp)              :: BC_beta_lo,BC_beta_hi
    real(dp)              :: BC_CDiag_lo,BC_CDiag_hi
    real(dp)              :: BC_COffDiag_lo,BC_COffDiag_hi

    real(dp)              :: MUE_lo,MUE_hi
    real(dp)              :: InvCDiag_lo,InvCDiag_hi
    real(dp)              :: InvCOffDiag_lo,InvCOffDiag_hi
  end type

  type(BayesType)             :: bayes ! information for Bayesian computations

  type(FlagStructure)         :: ControlOptions
  type(MaxStructure)          :: MaxOptions

  type(PenaltyStructureType)  :: Penalty
  type(QuadRuleType)          :: IntRule
  type(QuadNode)              :: RandomB
  type(SelectFreeType)        :: iFree
  type(ParmsStructure)        :: parms,parms0  ! parms  = current parameters
                                               ! parms0 = true parmaeters
  type(DataStructure)         :: HHData

  ! Output directory
  character(len=200)          :: OutDir
  character(len=200)          :: InputDir

  real(dp)                    :: inf
  real(dp)                    :: small

  integer(i4b), parameter     :: MasterID=0
  integer(i4b)                :: nWorkers

  ! Structure holding gradient of DensityFunc2 w.r.t. parameters
  type DensityGradType
    ! (S1,nu1,Omega12,V,C2,Q,CPsi)
    ! z1 = inv(S1) * VT * p1 + S1*VT*q1 - nu1 - Omega12 * inv(C2^T) * Q^T * x
    real(dp), allocatable :: S1(:,:)  ! gradient w.r.t. S1  : (nx x d1)
    real(dp), allocatable :: nu1(:,:) ! gradient w.r.t. nu1 : (nx x d1)
    real(dp), allocatable :: Omega12(:,:,:) ! gradient w.r.t. Omega12
                                            ! (nx x d1 x d2)
    real(dp), allocatable :: VT(:,:,:)    ! (nx x nV x nV)  (nV = size(V))
    real(dp), allocatable :: C2(:,:,:)    ! (nx x d2 x d2)
    real(dp), allocatable :: Q(:,:,:)     ! (nx x d2 x d2)
    real(dp), allocatable :: CPsi(:,:,:)  ! (nx x d1 x d1)
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

    ! IntRule
    !    IntRule%rule(i1) contains details to integrate an i1 dimensional integral
    !                     over random bliss points
    !                     diff rule for each i1 <= K
    IntRule%dMax = parms%K
    allocate(IntRule%flag(parms%K))
    allocate(IntRule%nall(parms%K))
    allocate(IntRule%rule(parms%K))
   
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

    if (LocalParms%model==2) then

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
    allocate(HHData%internet(HHData%N))
    allocate(HHData%SmallStore(HHData%N))

    if (parms%model==2) then
      allocate(HHData%eta(parms%dim_eta,HHData%N))
    end if
  end subroutine AllocateHHData

  ! allocate memory for data
  subroutine AllocateLocalHHData(LocalHHData)
    implicit none
    type(DataStructure), allocatable,intent(inout) :: LocalHHData

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
    allocate(LocalHHData%internet(LocalHHData%N))
    allocate(LocalHHData%SmallStore(LocalHHData%N))

    if (parms%model==2) then
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
    integer(i4b) :: ErrFlag,ierr
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
    else if  (parms%model==2) then
      nx = iFree%nBC_beta + iFree%nBC_CDiag + iFree%nBC_COffDiag &
         + iFree%nBD_beta + iFree%nBD_CDiag + iFree%nBD_COffDiag &
         + iFree%nMuE + iFree%nInvCDiag + iFree%nInvCOffDiag
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

subroutine ComputeNMC(pid,NMC,IMC1,IMC2)
  implicit none
  integer(i4b), intent(in)  :: pid,NMC
  integer(i4b), intent(out) :: IMC1,IMC2
  integer(i4b)              :: NMC1,NMC2
  integer(i4b)              :: N1,N2
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
  J1 = (/1:1+(N1-1)*NMC1:NMC1/)
  if (pid+1<=N1) then
    IMC1 = j1(pid+1)
    IMC2 = IMC1+NMC1-1;
  end if
  J2 = (/1:1+(N2-1)*NMC2:NMC2/)+N1*NMC1
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
  deallocate(HHData%market,HHData%e)
  deallocate(HHData%expenditure)
  deallocate(HHData%utility)
  deallocate(HHData%ColumnLabels)
  if (allocated(HHData%eta)) then
    deallocate(HHData%eta)
  end if


  ! Deallocate integration rule
  deallocate(IntRule%flag)
  deallocate(IntRule%nall)
  
  do i1=1,parms%K
    deallocate(IntRule%rule(i1)%nQuad)
    deallocate(IntRule%rule(i1)%nodes)
    deallocate(IntRule%rule(i1)%weights)
  end do
  deallocate(IntRule%rule)
  
  ! deallocate RandomB
  if (parms%model==2) then
    deallocate(RandomB%nQuad)
    deallocate(RandomB%nodes)
    deallocate(RandomB%weights)
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
    character(Len=100)              :: TempDir

    
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
  read(cTemp,'(i3)') MaxOptions%SaveBasisFreq

  ErrFlag = GetVal(PropList,'BasisFile',cTemp)
  MaxOptions%BasisFile = trim(OutDir) // '/' // trim(cTemp)
  
  ErrFlag = GetVal(PropList,'OldBasisFile',cTemp)
  MaxOptions%OldBasisFile = trim(OutDir) // '/' // trim(cTemp)

  ErrFlag = GetVal(PropList,'BackupBasisFile',cTemp)
  MaxOptions%BackupBasisFile = trim(OutDir) // '/' // trim(cTemp)

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
  TempDir = trim(OutDir) // '/parms'
  inquire(DIRECTORY=TempDir,exist=DirExists)
    if (.not. DirExists) then
      print *, 'Warning: parms directory does not exist.'
      print *, 'Creating ',trim(TempDir),'.'
      DirExists = makedirqq(TempDir)
      if (.not. DirExists) then
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

  if (parms%model==2) then
    ErrFlag = GetVal(PropList,'BC_z_dim',cTemp)
    read(cTemp,'(i5)') parms%BC_z_dim
    ErrFlag = GetVal(PropList,'BD_z_dim',cTemp)
    read(cTemp,'(i5)') parms%BD_z_dim

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

    ErrFlag = GetVal(PropList,'BD_beta_lo',cTemp)
    read(cTemp,'(f12.0)') parms%BD_beta_lo
    
    ErrFlag = GetVal(PropList,'BD_beta_hi',cTemp)
    read(cTemp,'(f12.0)') parms%BD_beta_hi

    ! RandomB parameters: dimension of random elements of C and D
    ErrFlag = GetVal(PropList,'dim_eta',cTemp)
    read(cTemp,'(i3)') parms%dim_eta
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

  ! define integration rule
  ! Define IntRule and RandomB
  call DefineIntegrationRule

  ! read in flags to select free parameters
  if (parms%model==1) then
    ErrFlag = GetVal(PropList,'FreeFlagD',cTemp)
    read(cTemp,'(i2)') iFree%flagD
    ErrFlag = GetVal(PropList,'FreeFlagBC',cTemp)
    read(cTemp,'(i2)') iFree%flagBC
  end if
  ErrFlag = GetVal(PropList,'FreeFlagMUE',cTemp)
  read(cTemp,'(i2)') iFree%flagMUE
  ErrFlag = GetVal(PropList,'FreeFlagInvCDiag',cTemp)
  read(cTemp,'(i2)') iFree%flagInvCDiag
  ErrFlag = GetVal(PropList,'FreeFlagInvCOffDiag',cTemp)
  read(cTemp,'(i2)') iFree%flagInvCOffDiag

  if (parms%model==2) then
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
  end if

  ErrFlag = GetVal(PropList,'OneAtATime',cTemp)
  read(cTemp,'(i2)') iFree%OneAtATime

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
    integer(i4b)       :: unit_D,unit_C,unit_MUE,unit_invCDiag,unit_INVCOffDiag
    character(len=200) :: file_D,file_C,file_MUE,file_INVCDiag,file_INVCOffDiag
    integer(i4b)       :: unit_BC_beta,unit_BC_CDiag,unit_BC_COffDiag
    integer(i4b)       :: unit_BD_beta,unit_BD_CDiag,unit_BD_COffDiag
    character(len=200) :: file_BC_beta,file_BC_CDiag,file_BC_COffDiag
    character(len=200) :: file_BD_beta,file_BD_CDiag,file_BD_COffDiag
    integer(i4b)       :: unit_sigp
    character(len=200) :: file_sigp

    integer(i4b), allocatable :: index(:)
    integer(i4b)              :: row,col
    real(dp),     allocatable :: temp2(:,:)
    integer(i4b)              :: i1,n

    if (parms%model==1) then
      ! D: norm of columns of B
      unit_D = 471
      file_D = trim(InputDir) // '/D.raw'
      open(UNIT = unit_D, &
           FILE = file_D, &
           ACTION = 'read')
      parms%D = 0.0d0

      do i1=1,parms%J
        read(unit_D,'(g25.16)') parms%D(i1)
      end do
      close(unit_D)

      ! C : off-diagonal elements of B
      unit_C = 484
      file_C = trim(InputDir) // '/C.raw'
      open(UNIT = unit_C,  &
           FILE = file_C,  &
           ACTION = 'read')
       parms%BC = 0.0d0
      do i1=1,parms%nBC
        read(unit_C,'(g25.16)') parms%BC(i1)
      end do
      close(unit_C)

      call SphereToMatrix(parms%BC,parms%D,parms%K,parms%J,parms%B)
    end if ! if (model==1)

    ! parms%MuE : mean of e
    unit_MUE = 496
    !file_MUE = trim(InputDir) // '/MUE.raw'
    file_MUE = trim(InputDir) // '/' // trim(ParmFiles%MUE)
    open(UNIT = unit_MUE,  &
         FILE = file_MUE,  &
         ACTION = 'read')
    parms%MUE = 0.0d0
    do i1=1,parms%K
      read(unit_MUE,'(g25.16)') parms%MUE(i1)
    end do
    close(unit_MUE)

    ! parms%InvCDiag : diagonal elements of inverse of C = chol(sig)
    unit_InvCDiag = 508
    !file_InvCDiag = trim(InputDir) // '/INVCDiag.raw'
    file_InvCDiag = trim(InputDir) // '/' // trim(ParmFiles%InvCDiag)
    open(UNIT = unit_InvCDiag,  &
         FILE = file_InvCDiag,  &
         ACTION = 'read')
    parms%InvCDiag = 0.0d0
    do i1=1,parms%K
      read(unit_InvCDiag,'(g25.16)') parms%InvCDiag(i1)
    end do
    close(unit_InvCDiag)

    ! parms%InvCOffDiag : off-diagonal elements of inverse of C = chol(sig)
    unit_InvCOffDiag = 520
    !file_InvCOffDiag = trim(InputDir) // '/INVCOffDiag.raw'
    file_InvCOffDiag = trim(InputDir) // '/' // trim(ParmFiles%InvCOffDiag)
    open(UNIT = unit_InvCOffDiag,  &
         FILE = file_InvCOffDiag,  &
         ACTION = 'read')
    parms%InvCOffDiag = 0.0d0
    do i1=1,parms%K*(parms%K-1)/2
      read(unit_InvCOffDiag,'(g25.16)') parms%InvCOffDiag(i1)
    end do
    close(unit_InvCOffDiag)
    call SphereToMatrix(parms%InvCOffDiag,parms%InvCDiag,parms%K,parms%K,parms%InvC)
    ! compute CSig = inv(InvC)
    parms%InvC = transpose(parms%InvC)
    call MatrixInverse(parms%InvC,parms%CSig,'Lower triangular')
    parms%sig  = matmul(parms%CSig,transpose(parms%CSig))

    if (parms%model==2) then
      ! read in (BC_beta,BD_beta)
      unit_BC_beta = 618
      !file_BC_beta = trim(InputDir) // '/BC_beta.raw'
      file_BC_beta = trim(InputDir) // '/' // trim(ParmFiles%BC_Beta)
      open(UNIT = unit_BC_beta,  &
           FILE = file_BC_beta,  &
           ACTION = 'read')
      parms%BC_beta = 0.0d0
      do i1=1,parms%BC_z_dim
        read(unit_BC_beta,'(g25.16)') parms%BC_beta(i1)
      end do
      close(unit_BC_beta)
    
      unit_BD_beta = 749
      !file_BD_beta = trim(InputDir) //  '/BD_beta.raw'
      file_BD_beta = trim(InputDir) //  '/' // trim(ParmFiles%BD_Beta)
      open(UNIT = unit_BD_beta,  &
           FILE = file_BD_beta,  &
           ACTION = 'read')
      parms%BD_beta = 0.0d0
      do i1=1,parms%BD_z_dim
        read(unit_BD_beta,'(g25.16)') parms%BD_beta(i1)
      end do
      close(unit_BD_beta)
   
      ! read (BC_CDiag,BD_CDiag)
      unit_BC_CDiag = 761
      !file_BC_CDiag = trim(InputDir) // '/BC_CDiag.raw'
      file_BC_CDiag = trim(InputDir) // '/' // trim(ParmFiles%BC_CDiag)
      open(UNIT = unit_BC_CDiag,  &
           FILE = file_BC_CDiag,  &
           ACTION = 'read')
      parms%BC_CDiag = 0.0d0
      do i1=1,parms%nBC
        read(unit_BC_CDiag,'(g25.16)') parms%BC_CDiag(i1)
      end do
      close(unit_BC_CDiag)

      unit_BD_CDiag = 771
      !file_BD_CDiag = trim(InputDir) // '/BD_CDiag.raw'
      file_BD_CDiag = trim(InputDir) // '/' // trim(ParmFiles%BD_CDiag)
      open(UNIT = unit_BD_CDiag,  &
           FILE = file_BD_CDiag,  &
           ACTION = 'read')
      parms%BD_CDiag = 0.0d0
      do i1=1,parms%J
        read(unit_BD_CDiag,'(g25.16)') parms%BD_CDiag(i1)
      end do
      close(unit_BD_CDiag)

      ! read (BC_COffDiag,BD_COffDiag)
      unit_BC_COffDiag = 782
      !file_BC_COffDiag = trim(InputDir) // '/BC_COffDiag.raw'
      file_BC_COffDiag = trim(InputDir) // '/' // trim(ParmFiles%BC_COffDiag)
      open(UNIT = unit_BC_COffDiag,  &
           FILE = file_BC_COffDiag,  &
           ACTION = 'read')
      parms%BC_COffDiag = 0.0d0
      do i1=1,parms%nBC_COffDiag
        read(unit_BC_COffDiag,'(g25.16)') parms%BC_COffDiag(i1)
      end do
      close(unit_BC_COffDiag)

      unit_BD_COffDiag = 792
      !file_BD_COffDiag = trim(InputDir) // '/BD_COffDiag.raw'
      file_BD_COffDiag = trim(InputDir) // '/' // trim(ParmFiles%BD_COffDiag)
      open(UNIT = unit_BD_COffDiag,  &
           FILE = file_BD_COffDiag,  &
           ACTION = 'read')
      parms%BD_COffDiag = 0.0d0
      do i1=1,parms%nBD_COffDiag
        read(unit_BD_COffDiag,'(g25.16)') parms%BD_COffDiag(i1)
      end do
      close(unit_BD_COffDiag)

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

      ! Compute parms%BD_z  (J x BD_z_dim) matrix
      !  default = identity matrix
      if ((size(parms%BD_Z,1) .ne. parms%J) .or. (size(parms%BD_Z,2) .ne. parms%J)) then
        print *,'Parms%BD is wrong size: ',size(parms%BD_Z,1),size(parms%BD_Z,2)
        stop
      end if
      parms%BD_z = 0.0d0
      do i1=1,parms%J
        parms%BD_z(i1,i1) = 1.0d0
      end do

      ! Compute parms%BC_Z   (nBC x BC_z_dim)
      parms%BC_z = 0.0d0
      allocate(index(parms%K-1))
      do i1=2,parms%J
        ! define BC_z so that y(i1) = BC_z(i1,:)*BC_beta = BC_beta(j) if row i1 of
        ! y  corresponds to product j

        if (i1 <= parms%K) then
          index(1:i1-1) = (i1-1)*(i1-2)/2+(/1:i1-1/)
          parms%BC_z(index(1:i1-1),i1-1) = 1.0d0
        else
          index = parms%K*(parms%K-1)/2 + (parms%K-1)*(i1-parms%K-1)+(/1:parms%K-1/)
          parms%BC_z(index,i1-1) = 1.0d0
        end if
      end do
      deallocate(index)
    end if ! if (model==2) then
   
    unit_sigp = 471
    file_sigp = trim(InputDir) // '/' // trim(ParmFiles%sigp)
    open(UNIT = unit_sigp, &
         FILE = file_sigp, &
         ACTION = 'read')

    parms%sigp = 0.0d0
    do i1=1,parms%J*parms%J
      row = (i1-1)/parms%J + 1
      col = i1-parms%J*(row-1)  
      read(unit_sigp,897) parms%sigp(row,col)
    end do
    897 format(f25.0)
    close(unit_sigp)
    
  end subroutine ReadParameters

subroutine DefineIntegrationRule
  implicit none
  integer(i4b)               :: ErrFlag
  integer(i4b)               :: n,i1
  character(len=PL_VAL_LEN)  :: cTemp      ! temporary character string
  integer(i4b), allocatable  :: nQuad(:)
  integer(i4b)               :: RandomB_flag,RandomB_nall,RandomB_dim,nRead


  ! integration rule
  ! 0 = Gauss hermite            0 : not working
  ! 1 = Sparse Hermite           1 : not working
  ! 2 = Pseudo Monte Carlo       2 : not working
  ! 3 = Quasi Monte Carlo        3 : not working
  ! 6 = Gauss-Legendre on (-1,1) 6 : working
  ! 7 = pseudo-MC on (-1,1)      7 : working

  ErrFlag = GetVal(PropList,'IntegrationFlag',cTemp)
  read(cTemp,'(<parms%K>i2)') IntRule%flag
 
  ErrFlag = GetVal(PropList,'nQuadAll',cTemp)
  read(cTemp,'(<parms%K>i4)') IntRule%nAll

  allocate(nQuad(parms%K))
  ErrFlag = GetVal(PropList,'nQuad',cTemp)
  read(cTemp,'(<parms%K>i3)') nQuad

  do i1=1,parms%K
    allocate(IntRule%rule(i1)%nQuad(i1))
    IntRule%rule(i1)%nQuad = nQuad(1:i1)
    n = product(IntRule%rule(i1)%nQuad)
    if (n>IntRule%nall(i1)) then
      n = IntRule%nall(i1)
    else
      IntRule%nall(i1) = n
    end if
    allocate(IntRule%rule(i1)%nodes(n,i1))
    allocate(IntRule%rule(i1)%weights(n))
    call DefineIntegrationNodes(i1,IntRule%flag(i1),IntRule%rule(i1)%nQuad,IntRule%nAll(i1), &
                                IntRule%rule(i1)%nodes,IntRule%rule(i1)%weights)
  end do

  deallocate(nQuad)

  if (parms%model==2) then
    ! Define RandomB rule
    ErrFlag = GetVal(PropList,'RandomB_flag',cTemp)
    read(cTemp,'(i2)') RandomB_flag
  
    ErrFlag = GetVal(PropList,'RandomB_nall',cTemp)
    read(cTemp,'(i4)') RandomB_nall

    RandomB_dim = parms%dim_eta 

    allocate(RandomB%nQuad(RandomB_dim))
    RandomB%nQuad = 3
    ErrFlag = GetVal(PropList,'RandomB_nQuad',cTemp)
    nRead = min(RandomB_dim,10)
    read(cTemp,'(<nRead>i3)') RandomB%nQuad(1:nRead)
  
    RandomB%nall =product(RandomB%nQuad)
    if (RandomB%nall>RandomB_nall) then
      RandomB%nall = RandomB_nall
    else if (RandomB%nall<=RandomB_nall) then
      RandomB_nall = RandomB%nall
    end if
  
    allocate(RandomB%nodes(RandomB%nall,RandomB_dim),RandomB%weights(RandomB%nall))
    call DefineIntegrationNodes(RandomB_Dim,RandomB_flag,RandomB%nQuad, &
                                RandomB%nall, &
                                RandomB%nodes,RandomB%weights)
  end if ! (parms%model==2) then
end subroutine DefineIntegrationRule 

subroutine DefineIntegrationNodes(d,flag,n1,nAll,nodes,weights)
  use ToolsModule, only : kron1
  use nr, only : gauher,gauleg
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

  external G05KFF  ! NAG routine to initialize random number generator
  external G05SAF  ! NAG routine to generate uniform random numbers
  external G05RZF  ! NAG routine to generate normal random numbers
  nodes   = 0.0d0
  weights = 1.0d0
  if (flag==0) then
    ! 0 = Gauss-Hermite
    do i1=1,d
      allocate(x(n1(i1)))
      allocate(w(n1(i1)))
      call gauher(x,w)
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
  
  ! initialize random number generator
  genid = 3   ! Mersenne twister algorithm
  subid = 0   ! not used when genid==3
  lseed = 1   ! number of seeds needed for genid==3
  allocate(seed(lseed))
  seed = 9824
  lstate = 1260  ! min value for genid==3
  allocate(state(lstate))
  ifail = 0
  call G05KFF(genid,subid,seed,lseed,state,lstate,ifail)

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
  call G05RZF(mode,nall,d,MU,C,d,R,LR,state, &
              nodes,nAll,ifail)
  weights = 1.0d0/dble(nAll)
  deallocate(seed)
  deallocate(state)
  deallocate(C,R)
elseif (flag==6) then
  ! Gauss-Legendre on (-1,1)
    do i1=1,d
      allocate(x(n1(i1)))
      allocate(w(n1(i1)))
      call gauleg(-1.0d0,1.0d0,x,w)
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
  
  ! initialize random number generator
  genid = 3   ! Mersenne twister algorithm
  subid = 0   ! not used when genid==3
  lseed = 1   ! number of seeds needed for genid==3
  allocate(seed(lseed))
  seed = 9824
  lstate = 633  ! min value for genid==3
  allocate(state(lstate))
  ifail = 0
  call G05KFF(genid,subid,seed,lseed,state,lstate,ifail)

  ! generate random numbers
  ifail = 0
  allocate(x(nall*d))
  call G05SAF(nall*d,state,x,ifail)
  x = 2.0d0*x-1.0d0
  nodes = reshape(x,(/nall,d/))
  deallocate(x)
  deallocate(seed)
  deallocate(state)

end if

end subroutine DefineIntegrationNodes

#if USE_MPI==1
subroutine BroadcastParameters(pid)
  use mpi
  use IFPORT
  implicit none
  integer(i4b), intent(in)  :: pid
  integer(i4b)              :: ierr(30),i1,i2,nerr,ierr_barrier
  integer(i4b), allocatable :: ierr_quad(:),ierr_B(:)
  integer(i4b)              :: n,d  ! size of matrix for RandomB and RandomD

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
  if (parms%model==2) then
    call mpi_bcast(parms%BC_z_dim,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(17))
    call mpi_bcast(parms%BD_z_dim,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(18))
    call mpi_bcast(parms%dim_eta,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(19))
    call mpi_bcast(parms%BC_lo,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(20))
    call mpi_bcast(parms%BC_hi,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%BC_beta_lo,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(20))
    call mpi_bcast(parms%BC_beta_hi,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%BD_beta_lo,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(20))
    call mpi_bcast(parms%BD_beta_hi,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%InvCDiag_LO,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%InvCDiag_HI,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%InvCOffDiag_LO,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
    call mpi_bcast(parms%InvCOffDiag_HI,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(21))
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
  call mpi_bcast(IntRule%flag,parms%K,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(25))
  call mpi_bcast(IntRule%nAll,parms%K,MPI_INteger,MasterID,MPI_COMM_WORLD,ierr(26))

  nerr = 2*parms%K+parms%K*parms%K
  allocate(ierr_quad(nerr))  
  ierr_quad = 0
  do i1=1,parms%K
    ! rule(i1)%nodes
    ! rule(i1)%weights
    ! rule(i1)%nQuad
    call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
    if (pid>MasterID) then
      allocate(IntRule%rule(i1)%nQuad(i1))
      allocate(IntRule%rule(i1)%nodes(IntRule%nAll(i1),i1))
      allocate(IntRule%rule(i1)%weights(IntRule%nAll(i1)))
    end if
    call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
    call mpi_bcast(IntRule%rule(i1)%nQuad,i1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr_quad(i1))
    call mpi_bcast(IntRule%rule(i1)%weights,IntRule%nAll(i1),MPI_DOUBLE_PRECISION, &
                MasterID,MPI_COMM_WORLD,ierr_quad(parms%K+i1))
    do i2=1,i1
      call mpi_bcast(IntRule%rule(i1)%nodes(:,i2),IntRule%nAll(i1),                &
                     MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr_quad(2*parms%K+parms%K*(i1-1)+i2))
    end do            
  end do
  print 1596,'ierr_quad',pid,ierr(25:26),ierr_quad
1596 format(a10,i4,2i3,<nerr>i3)
  deallocate(ierr_quad)
  call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)

  ! broadcast information on integration rule for RandomB
  call mpi_bcast(RandomB%nall,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr(27))
  if (pid==MasterID) then
    n = size(RandomB%nodes,1)
    d = size(RandomB%nodes,2)
  end if
  call mpi_bcast(n,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(28))
  call mpi_bcast(d,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(29))
  allocate(ierr_B(d+2))
  ierr_B = 0
  if (pid>MasterID) then
    allocate(RandomB%nodes(n,d))
    allocate(RandomB%weights(n))
    allocate(RandomB%nQuad(d))
  end if
 
  call mpi_barrier(MPI_COMM_WORLD,ierr_barrier) 
  call mpi_bcast(RandomB%weights,n,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr_B(1))
  call mpi_bcast(RandomB%nQuad,d,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr(2))
  do i1=1,d
    call mpi_bcast(RandomB%nodes(:,i1),n,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr_B(2+i1))
  end do
  print 1619,'ierr_B',pid,ierr_B
1619 format(a7,i4,<d+2>i3)
  deallocate(ierr_B)
  call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
end subroutine BroadcastParameters

subroutine BroadcastParms(LocalParms,pid)
  use mpi
  implicit none
  type(ParmsStructure), intent(inout) :: LocalParms
  integer(i4b), intent(in)            :: pid
  integer(i4b)                        :: i1,nerr,ierr_barrier
  integer(i4b), allocatable           :: ierr(:)

  allocate(ierr(5*LocalParms%J+3*LocalParms%K+11))
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
  do i1=1,LocalParms%k
    call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
    call mpi_bcast(LocalParms%InvC(:,i1),LocalParms%K,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+i1))
    call mpi_bcast(LocalParms%CSig(:,i1),LocalParms%K,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+Localparms%K+i1))
    call mpi_bcast(LocalParms%Sig(:,i1),LocalParms%K,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+2*Localparms%K+i1))
  end do
  nerr = nerr+3*LocalParms%K

  if (LocalParms%model==2) then
    call mpi_barrier(MPI_COMM_WORLD,ierr_barrier)
    call mpi_bcast(LocalParms%BD_beta,LocalParms%BD_z_dim,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+1))
    call mpi_bcast(LocalParms%BD_CDiag,LocalParms%J,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+2))
    call mpi_bcast(LocalParms%BD_COffDiag,LocalParms%nBD_COffDiag,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr(nerr+3))
    nerr = nerr+3
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
  end if

  nerr = 5*LocalParms%J+3*LocalParms%K+11
  print 1677,'parms_ierr',pid,ierr
1677 format(a11,i4,<nerr>i3)
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
  call mpi_bcast(iFree%flagInvCDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagInvCOffDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)

  call mpi_bcast(iFree%flagBC_beta,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagBC_CDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagBC_COffDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagBD_beta,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagBD_CDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%flagBD_COffDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%OneAtATime,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)

  call mpi_bcast(iFree%nD,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nBC,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nInvCDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nInvCOffDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nMuE,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)

  call mpi_bcast(iFree%nBD_beta,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nBD_CDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iFree%nBD_COffDiag,1,MPI_Integer,MasterID,MPI_COMM_WORLD,ierr)

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
  implicit none
  type(ParmsStructure), intent(inout) :: LocalParms
  character(LEN=*),     intent(in)    :: LocalAction

  integer(i4b)                        :: i1
  character(len=30)                   :: TempString

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
      read(LocalParms%unit,102) LocalParms%B(i1,:)
    end do

    ! update (D,BC): convert B to spherical coordinates (D,BC)
    call MatrixToSphere(LocalParms%B,LocalParms%D,LocalParms%BC)

    ! MUE   (K x 1)
    read(LocalParms%unit,113) TempString
    print *, TempString
    read(LocalParms%unit,103) LocalParms%MUE

    ! CSIG  (K x K)  Cholesky decomposition of sig
    read(LocalParms%unit,114) TempString
    print *,TempString
    do i1=1,LocalParms%K
      read(LocalParms%unit,104) LocalParms%CSIG(i1,:)
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
    read(LocalParms%unit,107) LocalParms%BC_BETA

    ! BD_beta  (BD_Z_DIM)
    read(LocalParms%unit,118) TempString
    print *, TempString
    read(LocalParms%unit,108) LocalParms%BD_BETA

    ! BC_C  (nBC x dim_eta)
    read(LocalParms%unit,119) TempString
    print *,TempString
    do i1=1,LocalParms%nBC
      read(LocalParms%unit,109) LocalParms%BC_C(i1,:)
    end do

    ! BD_C  (J x dim_eta)
    read(LocalParms%unit,119) TempString
    print *,TempString
    do i1=1,LocalParms%J
      read(LocalParms%unit,109) LocalParms%BD_C(i1,:)
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
      write(LocalParms%unit,102) LocalParms%B(i1,:)
    end do

    ! MUE
    write(LocalParms%unit,113) 'MUE(K x 1)'
    write(LocalParms%unit,103) LocalParms%MUE

    ! CSIG  (K x K)  Cholesky decomposition of sig
    write(LocalParms%unit,114) 'CSIG(K x K)'
    do i1=1,LocalParms%K
      write(LocalParms%unit,104) LocalParms%CSIG(i1,:)
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
    write(LocalParms%unit,107) LocalParms%BC_BETA

    ! BD_beta  (BD_Z_DIM)
    write(LocalParms%unit,118) 'BD_BETA(BD_Z_DIM x 1)'
    write(LocalParms%unit,108) LocalParms%BD_BETA

    ! BC_C  (nBC x dim_eta)
    write(LocalParms%unit,119) 'BC_C(nBC x dim_eta)'
    do i1=1,LocalParms%nBC
      write(LocalParms%unit,109) LocalParms%BC_C(i1,:)
    end do

    ! BD_C  (J x dim_eta)
    write(LocalParms%unit,119) 'BD_C(J x dim_eta)'
    do i1=1,LocalParms%J
      write(LocalParms%unit,109) LocalParms%BD_C(i1,:)
    end do

  end select

  close(LocalParms%unit)

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
  102 format(<LocalParms%J-1>(f25.16,','),f25.16)        ! B
  103 format(<LocalParms%K-1>(f25.16,','),f25.16)        ! MUE  (K x 1)
  104 format(<LocalParms%K-1>(f25.16,','),f25.16)        ! CSIG
  105 format(2(i4,','),i4)                               ! (dim_eta,BC_Z_DIM,BD_Z_DIM)
  106 format(2(i4,','),f25.16,',',f25.16)                ! nBC_COffDiag,nBD_COffDiag,BC_LO,BC_HI
  107 format(<LocalParms%BC_Z_DIM-1>(f25.16,','),f25.16) ! BC_BETA
  108 format(<LocalParms%BD_Z_DIM-1>(f25.16,','),f25.16) ! BD_BETA
  109 format(<LocalParms%dim_eta-1>(f25.16,','),f25.16)  ! BC_C

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
  parms_out%sig = parms_in%sig
  parms_out%CSIG = parms_in%CSIG
  parms_out%invc = parms_in%invc
  parms_out%InvCDiag = parms_in%InvCDiag
  parms_out%InvCOffDiag = parms_in%InvCOffDiag
  
  if (parms_out%model==2) then
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
    parms_out%GradBD_C_CDiag    = parms_in%GradBD_C_CDiag
    parms_out%GradBD_C_COffDiag = parms_in%GradBD_C_COffDiag

    parms_out%sigp = parms_in%sigp

  end if ! if (parms_out%model==2) then

  parms_out%BC_lo        = parms_in%BC_lo
  parms_out%BC_hi        = parms_in%BC_hi
  parms_out%InvCDiag_lo  = parms_in%InvCDiag_lo
  parms_out%InvCDiag_hi  = parms_in%InvCDiag_hi
  parms_out%InvCOffDiag_lo  = parms_in%InvCOffDiag_lo
  parms_out%InvCOffDiag_hi  = parms_in%InvCOffDiag_hi
  parms_out%BC_beta_lo  = parms_in%BC_beta_lo
  parms_out%BC_beta_hi  = parms_in%BC_beta_hi
  parms_out%BD_beta_lo  = parms_in%BD_beta_lo
  parms_out%BD_beta_hi  = parms_in%BD_beta_hi

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
  do i1=1,parms%J
    read(TaxParmsUnit,69) tax(i1,:)
  end do
 69 format(<ntax>g13.4)  

  deallocate(itemp)

  close(TaxParmsUnit)
end subroutine LoadTaxParameters

end module GlobalModule
