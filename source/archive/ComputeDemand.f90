subroutine ComputeDemand(HHData1)
  use GlobalModule, only : parms,DataStructure
  use nrtype
  implicit none
  type(DataStructure), intent(inout) :: HHData1

  ! variables used by E04NCA: solve quadratic program
  integer(i4b)               :: LCWSAV,LLWSAV,LIWSAV,LRWSAV
  integer(i4b), allocatable  :: IWSAV(:)
  real(dp),     allocatable  :: RWSAV(:)
  logical,      allocatable  :: LWSAV(:)
  character(6)               :: RNAME
  character(80), allocatable :: CWSAV(:)
  integer(i4b)               :: M,N,NCLIN,LDA
  integer(i4b), allocatable  :: istate(:),kx(:),iwork(:)
  integer(i4b)               :: iter,liwork,lwork,ifail
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

  ! solve quadratic program
  ! 1) E04NCA or E04NCF  (convex)
  ! 2) E04NFA or E04NFF  (not convex)
  ! 3) E04NKA or E04NKF  (sparse)
  ! variables used by E04NCA
  RNAME = 'E04NCA'
  LCWSAV = 1
  LLWSAV = 120
  LIWSAV = 610
  LRWSAV = 475

  allocate(IWSAV(LIWSAV))
  allocate(RWSAV(LRWSAV))
  allocate(LWSAV(LLWSAV))
  allocate(CWSAV(LCWSAV))

  ! initialize workspace for E04NCA
  call E04WBF(RNAME,CWSAV,LCWSAV,LWSAV,LLWSAV,IWSAV,LIWSAV,RWSAV,LRWSAV,ifail)

  ! set advisory message unit number to standard output (i.e. unit=6)
  err_unit = 6
  CALL X04ABF(1,err_unit)

  ! Open options file for reading
  ifail = -1
  mode = 0
  ! set options for E04NCA
  options_unit=7
  options_file= trim(InputDir) // '/E04NCA_Options.txt'
  call X04ACF(options_unit,options_file,mode,ifail)
  !open(UNIT = ioptions, &
  !     FILE = 'inputs/E04NCA_Options.txt', &
  !     ACTION = 'read')
  ifail = 0
  call E04NDA(options_unit,LWSAV,IWSAV,RWSAV,ifail)
  !close(ioptions)
  ifail = -1
  call X04ADF(options_unit,ifail)

  ! max -0.5*(B*q-e)'*(B*q-e) -p'*q
  !   subject to
  !   0 <=q <= q_max
  !
  ! OR
  ! min 0.5*(B*q-e)'*(B*q-e) + p'*q
  ! min 0.5*q'*B'*B*q + (-e'*B + p')*q +0.5*e'*e
  !  A = 0.5*B'*B
  !  CVEC = -B'*e + p
  !
  ! FOC    2*A*q -B'*e + p
  M = parms%J
  N = parms%J
  NCLIN = 0
  LDC   = 1
  LDA   = M
  allocate(C(LDC,1))
  C = 0.0d0
  allocate(BL(N),BU(N))
  BL = 0.0d0
  BU = 1.0d10
  allocate(CVEC(N))
  allocate(istate(n))
  allocate(KX(N))
  allocate(x(N))
  allocate(A(parms%J,parms%J))
  allocate(B(1))
  B = 0.0d0
  allocate(clamda(N))
  LIWORK = N
  allocate(IWORK(LIWORK))
  Lwork = 10*N
  allocate(work(LWORK))
  HHData1%q = 0.0d0
  crit = 1e-4

  do i1=1,HHData1%N
    if (parms%model==2) then
      ! if model==2, each HH TempB is a random coefficient
      call ComputeCurrentB(HHData1%eta(:,i1),parms)
    end if

    CVEC = HHData1%p(:,i1) - matmul(transpose(parms%B),HHData1%e(:,i1))
    x = 0.0d0
    ifail = -1
    ! compute X to solve
    !    min 0.5 * q'*A*q + cvec'*q   subject to q>=0
    ! E04NCA transforms A: so A needs to be reset to initial value

    A = matmul(transpose(parms%B),parms%B)
    call E04NCA(M,N,NCLIN,LDC,LDA,C,BL,BU,CVEC,ISTATE,KX,X,A,B,iter,OBJ,CLAMDA, &
                IWORK,LIWORK,WORK,LWORK,LWSAV,IWSAV,RWSAV,ifail)

    HHData1%nNonZero(i1) = count(x>=crit)
    HHData1%iNonZero(1:HHData1%nNonZero(i1),i1) = pack((/1:parms%J/),x>=crit)
    HHData1%iZero(1:parms%J-HHData1%nNonZero(i1),i1) = pack((/1:parms%J/),x<crit)
    if (HHData1%nNonZero(i1)>0) then
      HHData1%q(1:HHData1%nNonZero(i1),i1) = pack(x,x>=crit)
    end if
  end do

deallocate(IWSAV,RWSAV,LWSAV,CWSAV)
deallocate(C,BL,BU,CVEC,istate,kx,x,A,B,clamda)
deallocate(WORK,IWORK)
deallocate(seed,state)

end subroutine ComputeDemand