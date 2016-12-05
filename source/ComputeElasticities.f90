subroutine ComputeElasticities
  use nrtype
  use GlobalModule, only : parms,HHData,ControlOptions,InputDir, &
                           DataStructure,AllocateLocalHHData,    &
                           DeallocateLocalHHData
  use nag_library, only : X04ABF,X04ACF,X04ADF,  &
                          G05KFF,G05RZF,G05SAF,  &
                          E04WBF,E04NCA,E04NDA
  implicit none

  integer(i4b)          :: i1,i2
  type(DataStructure)     :: HHData0,HHData1

  ! initialize random number generator
  integer(i4b)   ::  genid,subid,lseed,ifail,lstate
  integer(i4b), allocatable :: seed(:),state(:)

  ! variables used by E04RZF: normal random number generator
  integer(i4b)               :: mode,ldc,ldx,LR
  real(dp), allocatable      :: R(:)
  real(dp), allocatable      :: e(:,:)

  ! variables used by E04NCA: solve quadratic program
  integer(i4b)               :: LCWSAV,LLWSAV,LIWSAV,LRWSAV
  integer(i4b), allocatable  :: IWSAV(:)
  real(dp),     allocatable  :: RWSAV(:)
  logical,      allocatable  :: LWSAV(:)
  character(6)               :: RNAME
  character(80), allocatable :: CWSAV(:)
  integer(i4b)               :: M,N,NCLIN,LDA
  integer(i4b), allocatable  :: istate(:),kx(:),iwork(:)
  integer(i4b)               :: iter,liwork,lwork
  real(dp), allocatable      :: C(:,:),BL(:),BU(:),CVEC(:),x(:),A(:,:),B(:)
  real(dp), allocatable      :: clamda(:),work(:)
  real(dp)                   :: obj
  real(dp)                   :: crit

  real(dp), allocatable      :: eye(:,:),zero(:)

  ! set options for E04NCA
  integer(i4b) :: options_unit
  integer(i4b) :: err_unit
  integer(i4b) :: PriceFlag
  character(len=200) :: options_file

  ! allocate memory for HHData0
  HHData0%N = HHData%N
  call AllocateLocalHHData(HHData0)

  ! set price = average price
  HHData0%p = spread(sum(HHData%p,2)/real(HHData0%n,8),2,HHData0%N)

  ! initialize random number generator
  !  external G05KFF  ! NAG routine to set seed for random number generator
  !  external G05SAF  ! NAG routine to generate uniform random numbers

  genid = 3 ! Mersenne twister algorithm
  subid = 0 ! not used when genid = 3
  !lseed = 624 ! number of seeds needed for genid = 3
  !allocate(seed(lseed))
  lseed = 1
  allocate(seed(lseed))
  lstate = 1260  ! min value for genid=3
  allocate(state(lstate))
  ifail = -1
  seed = IMC

  call G05KFF(genid,subid,seed,lseed,state,lstate,ifail)

  ! Generate HHData%e
  mode = 2
  ifail = 0
  ldc = parms%K
  ldx = parms%K
  LR = parms%K*(parms%K+1)+1
  allocate(R(LR))
  allocate(e(HHData0%N,parms%K))
  ! generate normal random numbers
  call G05RZF(mode,HHData0%N,parms%K,parms%MuE,parms%sig,parms%K,R,LR,state,e,HHData%N,ifail)
  HHData0%e = transpose(e)
  deallocate(R,e)

  ! Random coefficients in utility
  if (parms%model==2) then
    ! Random coefficients in BD:  eta
    ifail = 0
    LR = parms%dim_eta*(parm0%dim_eta+1)+1
    allocate(zero(parms%dim_eta))
    allocate(eye(parms%dim_eta,parms%dim_eta))
    allocate(e(HHData%N,parms%dim_eta))
    allocate(R(LR))
    R    = 0.0d0
    zero = 0.0d0
    eye  = 0.0d0
    do i1=1,parms%dim_eta
      eye(i1,i1) = 1.0d0
    end do
    call G05RZF(mode,HHData%N,parms%dim_eta,zero,eye,parms%dim_eta,R,LR,state,e,HHData0%N,ifail)
    HHData0%eta = transpose(e)
    deallocate(zero,eye,e,R)
  end if

  deallocate(seed,state)

  ! Compute baseline demand
  call ComputeDemand(HHData0)

  ! copy exogenous variables from HHData0 to HHData1
  HHData1%N = HHData0%N
  call AllocateLocalHHData(HHData1)
  HHData1%


  ! Compute initial
  ! deallocate memory for HHData0
  call DeallocateLocalHHData(HHData0)
  call DeallocateLocalHHData(HHData1)

end subroutine ComputeElasticities
