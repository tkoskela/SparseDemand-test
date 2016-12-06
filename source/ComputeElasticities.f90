subroutine ComputeElasticities
  use nrtype
  use GlobalModule, only : parms,HHData,ControlOptions,InputDir, &
                           DataStructure,AllocateLocalHHData,    &
                           DeallocateLocalHHData
  use nag_library, only  : G05KFF,G05RZF,G05SAF
  use OutputModule, only : WriteElasticities,WriteDemandResults
  implicit none

  integer(i4b)              :: i1,i2,i3
  type(DataStructure)       :: HHData0,HHData1

  ! variables to compute demand
  real(dp)                  :: h
  real(dp), allocatable     :: q0(:),q1(:)
  real(dp), allocatable     :: GradQ(:,:)
  real(dp), allocatable     :: elas(:,:)
  real(dp), allocatable     :: newQ(:,:),p(:)
  integer(i4b)              :: np

  ! initialize random number generator
  integer(i4b)              :: genid,subid,lseed,ifail,lstate
  integer(i4b), allocatable :: seed(:),state(:)

  ! variables used by E04RZF: normal random number generator
  integer(i4b)              :: mode,ldc,ldx,LR
  real(dp), allocatable     :: R(:)
  real(dp), allocatable     :: e(:,:)
  real(dp), allocatable     :: eye(:,:),zero(:)

  ! allocate memory for HHData0
  HHData0%N = HHData%N
  call AllocateLocalHHData(HHData0)

  ! set price = average price
  HHData0%p = spread(sum(HHData%p,2)/real(HHData0%n,8),2,HHData0%N)

  ! initialize random number generator
  !  external G05KFF  ! NAG routine to set seed for random number generator
  !  external G05SAF  ! NAG routine to generate uniform random numbers

  genid  = 3 ! Mersenne twister algorithm
  subid  = 0 ! not used when genid = 3
  !lseed = 624 ! number of seeds needed for genid = 3
  !allocate(seed(lseed))
  lseed  = 1
  allocate(seed(lseed))
  lstate = 1260  ! min value for genid=3
  allocate(state(lstate))
  ifail  = -1
  seed   = 1

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

  ! baseline demand
  allocate(q0(Parms%J),q1(parms%J))
  q0 = 0.0d0
  q1 = 0.0d0
  call ComputeAggregateDemand(HHData0%q,HHData0%iNonZero,q0)

  ! copy exogenous variables from HHData0 to HHData1
  HHData1%N = HHData0%N
  call AllocateLocalHHData(HHData1)
  HHData1%e   = HHData0%e
  HHData1%eta = HHData0%eta
  HHData1%p   = HHData0%p

  ! Compute slope of demand
  allocate(GradQ(parms%J,parms%J))
  allocate(elas(parms%J,parms%J))
  GradQ = 0.0d0
  elas  = 0.0d0
  h = 1.0d-4
  do i1=1,parms%J
    ! size(p) = (J x N)
    HHData1%p(i1,:) = HHData0%p(i1,:) + h
    call ComputeDemand(HHData1)
    call ComputeAggregateDemand(HHData1%q,HHData1%iNonZero,q1,HHData1%nNonZero,0)
    GradQ(:,i1) = (q1-q0)/(HHData1%p(i1,1)-HHData0%p(i1,1))
    elas(:,i1)  = HHData0%p(i1,1)*GradQ(:,i1)/q0
    HHData1%p(i1,:) = HHData0%p(i1,:)
  end do

  call WriteElasticities(elas)

  ! Plot demand functions w.r.t. own price
  !   for each good save
  !   p(i1), q(i1), q1,q2,q3,q4,q5,n1,n2,n3,n4,n5
  np = 30
  allocate(newq(np,parms%k+1),p(np))
  newq = 0.0d0
  p    = HHData0%p(1,1) * (0.5d0+(2.0d0-0.5d0)*real((/0:np-1/),dp)/real(np-1,dp))
  do i1=1,parms%J
    do i2=1,np
      HHData1%p(i1,:) = p(i2)
      call ComputeDemand(HHData1)
      do i3=0,5
        ! Compute q(i1) conditional on nNonZero==n
        !    if n==0, then compute total demand
        !    save results in matrix newq
        call ComputeAggregateDemand(HHData1%q,HHData1%iNonZero,newq(i2,i3+1),HHData1%nNonZero,i3)
      end do
    end do
    call WriteDemandResults(p,NewQ,i1)
  end do

  ! deallocate memory for HHData0
  call DeallocateLocalHHData(HHData0)
  call DeallocateLocalHHData(HHData1)
  deallocate(q0,q1,GradQ)
  deallocate(newq)

end subroutine ComputeElasticities

subroutine ComputeAggregateDemand(q,iNonZero,TotalQ,nNonZero,n)
  use nrtype
  use GlobalModule, only : parms
  implicit none
  real(dp),     intent(in)  :: q(:,:)
  integer(i4b), intent(in)  :: iNonZero(:,:)
  real(dp),     intent(out) :: TotalQ(:)
  integer(i4b), intent(in)  :: nNonZero(:),n

  integer(i4b) :: i1,i2

  do i1=1,parms%J
  do i2=1,parms%K
    if
    if (n==0) then
      TotalQ(i1) = sum(q(i2,:),mask=iNonZero(i2,:)==i1)
    elseif (n>0 .and. n<=parms%K) then
      TotalQ(i1) = sum(q(i2,:),mask=(iNonZero(i2,:)==i1 .and. nNonZero==n))
    endif
  end do
  end do
end subroutine ComputeAggregateDemand
