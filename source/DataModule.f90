module DataModule
! 1) subroutine CreateData(N,M,J,K,IMC)


contains

subroutine CreateData(IMC)
  use nrtype
  use GlobalModule, only : parms0,HHData,ControlOptions,InputDir
  use nag_library, only : X04ABF,X04ACF,X04ADF,  &
                          G05KFF,G05RZF,G05SAF,  &
                          E04WBF,E04NCA,E04NDA
  use OutputModule, only : SaveData
  implicit none
  integer(i4b),         intent(in) :: IMC

  integer(i4b)          :: i1,i2
  real(dp), allocatable :: p0(:),p1(:,:)
  real(dp), allocatable :: market(:)

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

  ifail = -1
  allocate(market(HHData%N))
  ! generate uniform random numbers
  call G05SAF(HHData%N,state,market,ifail)
  HHData%market = ceiling(real(HHData%M,dp)*market)
  allocate(p1(parms0%J,HHData%M))
 
  PriceFlag = 2
  if (PriceFlag==1) then
    ! prices are drawn from uniform distribution 
    n = HHData%M * parms0%J
    allocate(p0(n))
    ifail = 0
    call G05SAF(n,state,p0,ifail)
    p1 = reshape(p0,(/parms0%J,HHData%M/))
    deallocate(p0)
  else
    ! prices equal absolute value of a vector of normal random variables
    mode = 2
    ifail = 0
    ldc = parms0%J
    ldx = parms0%J
    LR = parms0%J*(parms0%J+1)+1
    allocate(R(LR))
    allocate(zero(parms0%J))
    zero = 0.0d0
    ! generate normal random numbers
    call G05RZF(mode,HHData%M,parms0%J,zero,parms0%sigp,parms0%J,R,LR,state,p1,HHData%M,ifail)
    p1 = abs(p1)
    deallocate(R,zero)
  end if

  ! if market==i1, set price equal to p(:,i1)
  do i1=1,HHData%N
    HHData%p(:,i1) = p1(:,HHData%market(i1))
  end do
  deallocate(p1)

  do i1=1,parms0%J
    write(HHData%ColumnLabels(i1),'(a8,i4)') 'Product ',i1
  end do
  ! compute HHData%e

  mode = 2
  ifail = 0
  ldc = parms0%K
  ldx = parms0%K
  LR = parms0%K*(parms0%K+1)+1
  allocate(R(LR))
  allocate(e(HHData%N,parms0%K))
  ! generate normal random numbers
  call G05RZF(mode,HHData%N,parms0%K,parms0%MuE,parms0%sig,parms0%K,R,LR,state,e,HHData%N,ifail)
  HHData%e = transpose(e)
  deallocate(R,e)

  ! Random coefficients in utility
  if (parms0%model==2) then

    ! Random coefficients in BD:  eta
    ifail = 0
    LR = parms0%dim_eta*(parms0%dim_eta+1)+1
    allocate(zero(parms0%dim_eta))
    allocate(eye(parms0%dim_eta,parms0%dim_eta))
    allocate(e(HHData%N,parms0%dim_eta))
    allocate(R(LR))
    R    = 0.0d0
    zero = 0.0d0
    eye  = 0.0d0
    do i1=1,parms0%dim_eta
      eye(i1,i1) = 1.0d0
    end do
    call G05RZF(mode,HHData%N,parms0%dim_eta,zero,eye,parms0%dim_eta,R,LR,state,e,HHData%N,ifail)
    HHData%eta = transpose(e)
    deallocate(zero,eye,e,R)

  end if
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
  M = parms0%J
  N = parms0%J
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
  allocate(A(parms0%J,parms0%J))
  allocate(B(1))
  B = 0.0d0
  allocate(clamda(N))
  LIWORK = N
  allocate(IWORK(LIWORK))
  Lwork = 10*N 
  allocate(work(LWORK))
  HHData%q = 0.0d0
  crit = 1e-4

  do i1=1,HHData%N
    if (parms0%model==2) then
      ! if model==2, each HH TempB is a random coefficient
      call ComputeCurrentB(HHData%eta(:,i1),parms0)
    end if

    CVEC = HHData%p(:,i1) - matmul(transpose(parms0%B),HHData%e(:,i1)) 
    x = 0.0d0
    ifail = -1
    ! compute X to solve
    !    min 0.5 * q'*A*q + cvec'*q   subject to q>=0
    ! E04NCA transforms A: so A needs to be reset to initial value

    A = matmul(transpose(parms0%B),parms0%B)
    call E04NCA(M,N,NCLIN,LDC,LDA,C,BL,BU,CVEC,ISTATE,KX,X,A,B,iter,OBJ,CLAMDA, &
                IWORK,LIWORK,WORK,LWORK,LWSAV,IWSAV,RWSAV,ifail)
    HHData%nNonZero(i1) = count(x>=crit)
    HHData%iNonZero(1:HHData%nNonZero(i1),i1) = pack((/1:parms0%J/),x>=crit)
    HHData%iZero(1:parms0%J-HHData%nNonZero(i1),i1) = pack((/1:parms0%J/),x<crit)
    if (HHData%nNonZero(i1)>0) then 
      HHData%q(1:HHData%nNonZero(i1),i1) = pack(x,x>=crit)
    end if

  end do

  if (ControlOptions%SaveDataFlag==1) then
    call SaveData
  end if

  deallocate(IWSAV,RWSAV,LWSAV,CWSAV)
  deallocate(C,BL,BU,CVEC,istate,kx,x,A,B,clamda)
  deallocate(WORK,IWORK)
  deallocate(market)
  deallocate(seed,state)
end subroutine CreateData

subroutine ComputeCurrentB(eta,parms)
  use nrtype
  use GlobalModule, only : ParmsStructure   ! data type for parms
  use NewTools,     only : SphereToMatrix   ! map spherical coordinates to B
  use nag_library,  only : S15ABF           ! normal CDF
  
  implicit none
  real(dp),             intent(in)  :: eta(:)
  type(ParmsStructure), intent(inout) :: parms
  
  real(dp), allocatable :: GradBC(:,:),GradBD(:,:)
  real(dp), allocatable :: BC(:),D(:)
  integer(i4b)          :: ifail,i1
  
  allocate(GradBC(parms%K,parms%nBC))
  allocate(GradBD(parms%K,parms%J))
  allocate(BC(parms%nBC))
  allocate(D(parms%J))
    
  ! eta = random coefficients in (BC,BD)

  parms%B = 0.0d0
  D = exp(matmul(parms%BD_z,parms%BD_beta)    &
          +matmul(parms%BD_C,eta))

  ! S15ABF = normcdf(x,ifail)
  do i1=1,parms%nBC
    BC(i1) = parms%BC_lo + (parms%BC_hi - parms%BC_lo) * pi_d &
                         * S15ABF(dot_product(parms%BC_z(i1,:),parms%BC_beta) &
                                  + dot_product(parms%BC_C(i1,:),eta),ifail)
  end do
  call SphereToMatrix(BC,D,parms%K,parms%J,parms%B,GradBC,GradBD)
  if (any(isnan(parms%B))) then
    print *, "Some values of B are NAN."
    print *,"BC = ",BC
    print *,"D = ",D
    print *,"B = ",parms%B
  end if
  deallocate(GradBC,GradBD)
  deallocate(BC,D)
end subroutine ComputeCurrentB

#if USE_MPI==1
subroutine SendData(pid)
  use nrtype
  use IFPORT
  use mpi
  use GlobalModule, only : HHData,parms,AllocateHHData,MasterID,nWorkers
  implicit none
  integer(i4b), intent(in) :: pid
  integer(i4b), allocatable :: stat_array(:,:)
  integer(i4b), allocatable :: stat1(:)

  ! send subset of data to each pid
  ! HHData%N = sample size
  integer(i4b) :: N1,N2,R,iHH1,iHH2,iw,i1,ierr(6)
  real(dp),     allocatable :: temp1(:)
  integer(i4b), allocatable :: itemp1(:)
  integer(i4b), allocatable :: request(:)

  ! N1+1 = size of data sent to processors (1:R)
  ! N1   = size of data sent to processors (R+1:N)
  ! R    = number of remainder data points after subtracting numprocs * N1
  N1 = HHData%N/nWorkers     ! integer division
  R = HHData%N - nWorkers*N1
  allocate(stat_array(MPI_STATUS_SIZE,12*nworkers))
  allocate(stat1(MPI_STATUS_SIZE))
  allocate(request(12*nworkers))

  if (pid==MasterID) then
    ! send data
    ! iHH1 = index of first element in data to send
    ! iHH2 = index of last element in data to send
    print *, "Begin sending data to workers."
    do iw=1,nWorkers
      if (R==0) then
        ! Send N1 observations to each processor
        iHH1 = N1*(iw-1)+1
        iHH2 = iHH1 + N1-1
      elseif (R>0) then
        if (iw<=R) then
          ! Send (N1+1) observations to first R processors
          iHH1 = (N1+1)*(iw-1)+1
          iHH2 = iHH1 + N1
        else if (iw>R) then
          ! send N1 observations to remainder of processors
          iHH1 = (N1+1)*R + N1*(iw-R-1) + 1
          iHH2 = iHH1 + N1-1
        end if
      end if

      N2 = merge(N1,N1+1,iw>R)
      ! send q
      allocate(temp1(parms%J*N2))
      temp1(1:parms%K*N2) = reshape(HHData%q(:,iHH1:iHH2),(/parms%K*N2/))
      call mpi_send(temp1(1:parms%K*N2),parms%K*N2,MPI_DOUBLE_PRECISION,iw,1,MPI_COMM_WORLD,ierr(1))
      !call mpi_isend(temp1(1:parms%K*N2),parms%K*N2,MPI_DOUBLE_PRECISION,iw,1,MPI_COMM_WORLD, &
      !               request(6*(iw-1)+1),ierr)
   
      ! send p 
      temp1 = reshape(HHData%p(:,iHH1:iHH2),(/parms%J*N2/))
      call mpi_send(temp1,parms%J*N2,MPI_DOUBLE_PRECISION,iw,2,MPI_COMM_WORLD,ierr(2))
      !call mpi_isend(temp1,parms%J*N2,MPI_DOUBLE_PRECISION,iw,2,MPI_COMM_WORLD,               &
      !               request(6*(iw-1)+2),ierr)
      deallocate(temp1)
 
      ! send market
      allocate(itemp1(parms%J*N2))
      itemp1(1:N2) = HHData%market(iHH1:iHH2)
      call mpi_send(itemp1(1:N2),N2,MPI_INTEGER,iw,3,MPI_COMM_WORLD,ierr(3))
      !call mpi_isend(itemp1(1:N2),N2,MPI_INTEGER,iw,3,MPI_COMM_WORLD,      &
      !               request(6*(iw-1)+3),ierr)

      ! send iNonZero
      itemp1(1:parms%K*N2) = reshape(HHData%iNonZero(:,iHH1:iHH2),(/parms%K*N2/))
      call mpi_send(itemp1(1:N2*parms%K),N2*parms%K,MPI_INTEGER,iw,4,MPI_COMM_WORLD,ierr(4))
      !call mpi_isend(itemp1(1:N2*parms%K),N2*parms%K,MPI_INTEGER,iw,4,MPI_COMM_WORLD, &
      !               request(6*(iw-1)+4),ierr)

      ! end iZero
      itemp1               = reshape(HHData%iZero(:,iHH1:iHH2),(/parms%J*N2/))
      call mpi_send(itemp1,N2*parms%J,MPI_INTEGER,iw,5,MPI_COMM_WORLD,ierr(5))
      !call mpi_isend(itemp1,N2*parms%J,MPI_INTEGER,iw,5,MPI_COMM_WORLD, &
      !               request(6*(iw-1)+5),ierr)

      ! end nNonZero
      itemp1(1:N2) = HHData%nNonZero(iHH1:iHH2)
      call mpi_send(itemp1(1:N2),N2,MPI_INTEGER,iw,6,MPI_COMM_WORLD,ierr(6))
      !call mpi_isend(itemp1(1:N2),N2,MPI_INTEGER,iw,6,MPI_COMM_WORLD, &
      !               request(6*(iw-1)+6),ierr)
      deallocate(itemp1)
      print 380,'datasend',pid,ierr
380 format(a10,i4,6i3)
      !call sleep(1)
      !do i1=1,6
      !  call mpi_wait(request(6*(iw-1)+i1),stat_array(:,6*(iw-1)+i1),ierr)
      !end do
    end do
  else if (pid .ne. MasterID) then
    ! receive subset of data
    N2 = merge(N1+1,N1,pid<=R)
    HHData%N = N2
    ! allocate memory for HHData
    call AllocateHHData
    allocate(temp1(parms%J*N2))
    
    ! receive q
    call mpi_recv(temp1(1:parms%K*N2),parms%K*N2,MPI_DOUBLE_PRECISION,MasterID,1,MPI_COMM_WORLD,stat1,ierr(1))
    !call mpi_irecv(temp1(1:parms%K*N2),parms%K*N2,MPI_DOUBLE_PRECISION,MasterID,1,MPI_COMM_WORLD, &
    !               request(6*(nworkers+pid-1)+1),ierr)
    HHData%q = reshape(temp1(1:parms%K*N2),(/parms%K,n2/))
   
    ! receive p 
    call mpi_recv(temp1,parms%J*N2,MPI_DOUBLE_PRECISION,MasterID,2,MPI_COMM_WORLD,stat1,ierr(2))
    !call mpi_irecv(temp1,parms%J*N2,MPI_DOUBLE_PRECISION,MasterID,2,MPI_COMM_WORLD, &
    !               request(6*(nworkers+pid-1)+2),ierr)
    HHData%p = reshape(temp1,(/parms%J,N2/))
    deallocate(temp1)

    allocate(iTemp1(N2*parms%J))
    ! receive market
    call mpi_recv(itemp1(1:N2),N2,MPI_INTEGER,MasterID,3,MPI_COMM_WORLD,stat1,ierr(3))
    !call mpi_irecv(itemp1(1:N2),N2,MPI_INTEGER,MasterID,3,MPI_COMM_WORLD, &
    !               request(6*(nworkers+pid-1)+3),ierr)
    HHData%market = itemp1(1:N2)

    ! receive iNonZero
    call mpi_recv(itemp1(1:parms%K*N2),parms%K*N2,MPI_INTEGER,MasterID,4,MPI_COMM_WORLD,stat1,ierr(4))
    !call mpi_irecv(itemp1(1:parms%K*N2),parms%K*N2,MPI_INTEGER,MasterID,4,MPI_COMM_WORLD, &
    !               request(6*(nworkers+pid-1)+4),ierr)
    HHData%iNonZero = reshape(itemp1(1:parms%K*N2),(/parms%K,N2/))

    ! receive iZero
    call mpi_recv(itemp1,parms%J*N2,MPI_INTEGER,MasterID,5,MPI_COMM_WORLD,stat1,ierr(5))
    !call mpi_irecv(itemp1,parms%J*N2,MPI_INTEGER,MasterID,5,MPI_COMM_WORLD, &
    !               request(6*(nworkers+pid-1)+5),ierr)
    HHData%iZero = reshape(itemp1,(/parms%J,N2/))

    ! receive nNonZero
    call mpi_recv(itemp1(1:N2),N2,MPI_INTEGER,MasterID,6,MPI_COMM_WORLD,stat1,ierr(6))
    !call mpi_irecv(itemp1(1:N2),N2,MPI_INTEGER,MasterID,6,MPI_COMM_WORLD, &
    !               request(6*(nworkers+pid-1)+6),ierr)
    HHData%nNonZero = itemp1(1:N2)
    deallocate(itemp1) 
    print 433,'receive:',pid,ierr
433 format(a10,i4,6i3)
    !do i1=1,6
    !  call mpi_wait(request(6*(nworkers+pid-1)+i1),stat_array(:,6*(nworkers+pid-1)+i1),ierr)
    !end do
  end if

  !call mpi_waitall(12*nworkers,request,stat_array,ierr) 
  deallocate(stat_array,request)
end subroutine SendData
#endif

subroutine LoadData
  use nrtype
  use GlobalModule, only : HHData,parms
  implicit none
  integer(i4b)  :: DataUnit,iComma,i1
  character(400) :: AllLabels
  real(dp),     allocatable :: qp(:,:),err(:,:),tiering(:,:),expenditure(:)
  integer(i4b), allocatable :: fascia(:)

  DataUnit = 20
  open(unit = DataUnit, &
       file = HHData%RawDataFile, &
       action = 'read')
 
  ! read in variable labels
  read(DataUnit,371) AllLabels
371 format(a)

  ! copy variable labels to a vector of characters
  allocate(HHData%RawDataLabels(HHData%nRawVars))
  do i1=1,HHData%nRawVars
    ! iCOmma = end of current data field
    iComma = index(AllLabels,',')
    if (iComma>0) then
      HHData%RawDataLabels(i1) = AllLabels(1:iComma-1)
      AllLabels = AllLabels(iComma+1:400)
    else 
      HHData%RawDataLabels(i1) = AllLabels
      exit
    end if
  end do

  allocate(qp(HHData%N,2*parms%J))
  allocate(tiering(HHData%N,parms%J))
  allocate(expenditure(HHData%N))
  allocate(fascia(HHData%N))
  allocate(err(HHData%N,parms%J))

  if (HHData%RawDataFormat==1) then
    do i1=1,HHData%N
      read(DataUnit,389) HHData%HHID(i1),HHData%date(i1),HHData%shopid(i1),  &
                         qp(i1,:),HHData%nNonZero(i1),HHData%day(i1),err(i1,:)
      HHData%q(1:HHData%nNonZero(i1),i1) = pack(qp(i1,1:2*parms%J:2),qp(i1,1:2*parms%J:2)>0.0d0)
      HHData%p(:,i1) = qp(i1,2:2*parms%J:2)
      HHData%iNonZero(1:HHData%nNonZero(i1),i1) = pack((/1:parms%J/),qp(i1,1:2*parms%J:2)>0.0d0)
      HHData%iZero(1:parms%J-HHData%nNonZero(i1),i1) = pack((/1:parms%J/),qp(i1,1:2*parms%J:2)==0.0d0)
    end do
389 format(3i10,<2*parms%J>f25.0,2i10,<parms%J>f25.0)
  elseif (HHData%RawDataFormat==2) then
    do i1=1,HHData%N
      read(DataUnit,390) HHData%HHID(i1),HHData%date(i1),HHData%shopid(i1),  &
                         qp(i1,:),HHData%nNonZero(i1),         &
                         fascia(i1),expenditure(i1),HHData%day(i1),err(i1,:)
      HHData%q(1:HHData%nNonZero(i1),i1) = pack(qp(i1,1:2*parms%J:2),qp(i1,1:2*parms%J:2)>0.0d0)
      HHData%p(:,i1) = qp(i1,2:2*parms%J:2)
      HHData%iNonZero(1:HHData%nNonZero(i1),i1) = pack((/1:parms%J/),qp(i1,1:2*parms%J:2)>0.0d0)
      HHData%iZero(1:parms%J-HHData%nNonZero(i1),i1) = pack((/1:parms%J/),qp(i1,1:2*parms%J:2)==0.0d0)
    end do
390 format(3i10,<2*parms%J>f25.0,2i10,f25.0,i10,<parms%J>f25.0)
  end if 
  deallocate(qp,err)
  close(DataUnit)

end subroutine LoadData

end module DataModule
