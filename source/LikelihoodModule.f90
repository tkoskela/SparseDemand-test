module LikelihoodModule
! Revision history
! 2017SEP13 LN  add option to max one at time
! 10JUN2014 LN  continue editing random coefficients model
! 22MAR2013 LN  clean up
! 09DEC2012 LN  adapted from Like1.m
!
! Subroutines
! Line No. |  Name
!-------------------------
!
! 41       | subroutine Like2(mode,nx,x,L,GradL,nstate,iuser,ruser)
! 18       | subroutine Like1(mode,nx,x,L,GradL,nstate,iuser,ruser)                        
! 490      | subroutine UpdateParms(xFree,iFree,parms)
! 566      | subroutine LogDensityFunc(e,F,GradF_e,GradF_MuE,GradF_InvC)
! 634      | subroutine ComputeRowGroup(R,J,d,RowGroup)
! 650      | subroutine ComputeProb(p,v,w,RowGroup,R,D,Integrand,    &
! 818      | subroutine DensityFunc2(x,G0,G1,Q,MuE,S12,C22,omega11,d1,d2,i1,q1,F,GradF)
! 1816     ! subroutine Constraint(x,mode,C,GradC)
use nrtype
! variables used by monfun_E04JCF
real(dp), allocatable :: last_x(:)

contains

#if USE_MPI>0
subroutine Like2_master(mode,nx,x,L,GradL,nstate,iuser,ruser)
  use nrtype
  use mpi
  use GlobalModule, only : MasterID,nWorkers
  implicit none
  
  integer(i4b), intent(inout)        :: mode
  integer(i4b), intent(in)           :: nx
  real(dp),     intent(in)           :: x(:)
  real(dp),     intent(out)          :: L
  real(dp),     intent(inout)        :: GradL(:)
  integer(i4b), intent(in)           :: nstate
  integer(i4b), intent(in)           :: iuser(*)
  real(dp),     intent(in)           :: ruser(*)

  integer(i4b)                       :: task,ierr,received
  real(dp)                           :: L1
  real(dp), allocatable              :: result(:)
  integer(i4b)                       :: status(MPI_STATUS_SIZE)
  integer(i4b)                       :: mode_in

  mode_in=mode
  !mode = 0  ! set mode=0 since gradient is not yet defined

  L = 0.0d0
  if (mode>0) then
    !GradL = 0.0d0
    allocate(result(nx+1))
  end if

  ! task=1   compute likelihood and possibly gradient
  task = 1
  call mpi_bcast(task,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(x,nx,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(mode,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)

  received=0
  do while (received<nWorkers)
    if (mode==0) then
      ! level of likelihood only
      call mpi_recv(L1,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
      L = L+L1
    elseif (mode>0) then
      call mpi_recv(result,nx+1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
      L  = L + result(1)
    !  GradL = GradL + result(2:nx+1)
    end if
    received=received+1
  end do

  if (mode>0) then
    deallocate(result)
  end if

  !mode=mode_in

end subroutine Like2_master

! Evaluate worker tasks
! task=0:  worker receives (x,mode) and computes likelihood for subset of data points
subroutine WorkerTask(model,pid)
  use nrtype
  use mpi
  use GlobalModule, only : MasterID,nWorkers
  implicit none

  integer(i4b), intent(in)  :: model,pid

  integer(i4b)              :: task,mode,ierr
  real(dp)                  :: L
  real(dp),     allocatable :: x(:),GradL(:)
  integer(i4b)              :: stat(MPI_STATUS_SIZE)
  integer(i4b), allocatable :: iuser(:)
  real(dp),     allocatable :: ruser(:)

  integer(i4b)              :: ntemp,nstate,nx

  call mpi_bcast(nx,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ntemp,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  allocate(iuser(ntemp))
  call mpi_bcast(iuser,ntemp,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ntemp,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  allocate(ruser(ntemp))
  call mpi_bcast(ruser,ntemp,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr)

  allocate(x(nx))
  allocate(GradL(nx))

  task=1

  do while (task==1)
    ! task : dictates which job worker should do
    ! task = 1 : compute likelihood
    ! task = 0 : finish
    call mpi_bcast(task,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)

    if (task==1) then
      ! Evaluate objective function

      ! broadcast parameter vector x and mode
      call mpi_bcast(x,nx,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr)
      call mpi_bcast(mode,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)

      if (model==2) then
        call Like2(mode,nx,x,L,GradL,nstate,iuser,ruser)
        if (mode==0) then
          call mpi_send(L,1,MPI_DOUBLE_PRECISION,MasterID,mode,MPI_COMM_WORLD,ierr)
        elseif (mode>0) then
          call mpi_send((/L,GradL/),nx+1,MPI_DOUBLE_PRECISION,MasterID,mode,MPI_COMM_WORLD,ierr)
        end if
      end if
    end if
  end do

  deallocate(x,GradL)
  deallocate(iuser,ruser)
end subroutine WorkerTask
#endif
!-----------------------------------------------------------------------------
!
!  LIKE2:   Likelihood function:  Random coefficients quadratic utility
!
!   Compute likelihood and its gradient w.r.t. x where utility is  
!
!        u  = y - p'*q - 0.5 * (B*q - e)' * (B*q - e)
!
!   with (B,e) both random
!
!
!  Revision history
!   2014JUN10  LN  completed first version
!                  Tasks:
!                  1) edit UpdateParms_RC
!                  2) test
!                  3) gradient
!
!-----------------------------------------------------------------------------
subroutine Like2(mode,nx,x,L,GradL,nstate,iuser,ruser)
  use nrtype
  use GlobalModule, only : iFree,parms,HHData,RandomB,small
  use DataModule,   only : ComputeCurrentB
  implicit none
  
  integer(i4b), intent(inout)        :: mode
  integer(i4b), intent(in)           :: nx
  real(dp),     intent(in)           :: x(:)
  real(dp),     intent(out)          :: L
  real(dp),     intent(inout)        :: GradL(:)
  integer(i4b), intent(in)           :: nstate
  integer(i4b), intent(in)           :: iuser(*)
  real(dp),     intent(in)           :: ruser(*)
  
  integer(i4b)                       :: d1,d2,d3,iHH,iq
  real(dp)                           :: L0,L1
  real(dp), allocatable              :: GradL0(:),GradL1(:)

  call UpdateParms2(x,iFree,parms)
  
  ! Initial values for likelihood and gradient
  L     = 0.0d0
  !GradL = 0.0d0

  ! Loop through households
  
  allocate(GradL0(nx),GradL1(nx))
  do iHH=1,HHData%N
    L0 = 0.0d0
    GradL0 = 0.0d0

    d1 = HHData%nNonZero(iHH)
    d2 = parms%K - d1
    d3 = parms%J-d1

    do iq=1,RandomB%nall
      ! update parms%B
      call ComputeCurrentB(RandomB%nodes(iq,:),parms,HHData%month(iHH))
      call Like1_wrapper(iHH,d1,d2,d3,parms,mode,L1,GradL1)
      if (L1>0.0d0) then
        L0 = L0 + RandomB%weights(iq)*L1
        GradL0 = GradL0 + RandomB%weights(iq)*GradL1
      end if
    end do
    ! avoid log of zero
    L = L + dlog(small+L0)
    if (mode>0) then
    !  GradL = GradL + GradL0/L0
    end if
  end do

  L = -L 
  deallocate(GradL0,GradL1)
  
end subroutine Like2

!-----------------------------------------------------------------------------
!
!  LIKE2Hess:   Compute vector of log-likelihood values:
!               one element for each observation in sample
!               (used to compute information matrix approximation)
!
!-----------------------------------------------------------------------------
subroutine Like2Hess(nx,x,L)
  use nrtype
  use GlobalModule, only : iFree,parms,HHData,RandomB,small
  use DataModule,   only : ComputeCurrentB
  implicit none
  
  integer(i4b), intent(in)           :: nx
  real(dp),     intent(in)           :: x(:)
  real(dp),     intent(out)          :: L(:)

  integer(i4b)                       :: d1,d2,d3,iHH,iq
  real(dp)                           :: L1
  real(dp), allocatable              :: GradL1(:)

  integer(i4b)                       :: mode

  mode = 0 ! only compute L(iHH)

  call UpdateParms2(x,iFree,parms)
  
  ! Initial values for likelihood
  L     = 0.0d0

  ! Loop through households
  
  allocate(GradL1(nx))
  do iHH=1,HHData%N

    d1 = HHData%nNonZero(iHH)
    d2 = parms%K - d1
    d3 = parms%J-d1

    do iq=1,RandomB%nall
      ! update parms%B
      call ComputeCurrentB(RandomB%nodes(iq,:),parms,HHData%month(iHH))
      call Like1_wrapper(iHH,d1,d2,d3,parms,mode,L1,GradL1)
      if (L1>0.0d0) then
        L(iHH) = L(iHH) + RandomB%weights(iq)*L1
      end if
    end do
    ! avoid log of zero
    L(iHH) = dlog(small+L(iHH))
  end do

  deallocate(GradL1)
  
end subroutine Like2Hess

!-----------------------------------------------------------------------------
!
! Like1_wrapper:   Given (i1,B), compute likelihood for 3 cases 
!
!-----------------------------------------------------------------------------
subroutine Like1_wrapper(i1,d1,d2,d3,parms,mode,L1,GradL1)
  use nrtype
  use GlobalModule, only : HHData,ParmsStructure

  implicit none
  integer(i4b),         intent(in)    :: i1,d1,d2,d3
  type(ParmsStructure), intent(inout) :: parms
  integer(i4b),         intent(inout) :: mode
  real(dp),             intent(out)   :: L1,GradL1(:)

  if (d1==parms%K) then
    ! Case A: Simple case:   dim(q>0) = K
    call Like1A(i1,parms,mode,L1,GradL1)
  else if (d1< parms%K .and. d1>0) then
    ! Case B:  dim(q>0)<K but >0
    call Like1B(i1,d1,d2,d3,parms,mode,L1,GradL1)
  else if (d1==0) then
    ! Case C:  q = 0
    call Like1C(i1,parms,mode,L1,GradL1)
  end if
end subroutine Like1_wrapper

!-----------------------------------------------------------------------------
!
!  Like1C:  LIKELIHOOD:    CASE 3:  d1 == 0
!
!           (HHData%nNonZero(i1) == d1 and d1 == 0)
!           Mapping from 0 to e is NOT one-to-one
!           need to integrate across region of e-space
!           satisfying B'*e<=p
!
!-----------------------------------------------------------------------------
subroutine Like1C(i1,parms,mode,L,GradL)
  use nrtype
  use GlobalModule, only : ParmsStructure,HHData,IntRule
  use NewTools,     only : ComputeLQ
  implicit none
  integer(i4b),         intent(in)    :: i1
  type(ParmsStructure), intent(in)    :: parms
  integer(i4b),         intent(inout) :: mode
  real(dp),             intent(out)   :: L
  real(dp),             intent(inout) :: GradL(:)

  real(dp)                  :: M1(parms%J,parms%K)
  real(dp)                  :: DTIlde(parms%J)

  ! M1 = LL*Q
  real(dp)                  :: R(parms%J,parms%K),Q(parms%K,parms%K)
  integer(i4b)              :: RowGroup(parms%J)

  integer(i4b)              :: integrand
  integer(i4b)              :: ifail
 
  M1 = matmul(transpose(parms%B),parms%CSig)

  ! M1 = LL*Q
  !      LL = (J x K) lower triangular
  call ComputeLQ(parms%J,parms%K,M1,R,Q,ifail)
  DTilde = HHData%p(:,i1) - matmul(transpose(parms%B),parms%MuE+parms%mue_month(:,HHData%month(i1)))

  ! Prob = integral of DensityFunc over region of x satisfying R*x<=M2_tilda
  ! 
  Integrand = 3
    
  ! size(RowGroup) = d2 = parms%J
  ! RowGroup(j1) = max(i1) s.t.  R(i1,j1) .ne. 0.0d0
  call ComputeRowGroup(R,parms%J,parms%K,RowGroup)

  call ComputeProb(L,IntRule%rule(parms%K)%nodes,IntRule%rule(parms%K)%weights, &
                   RowGroup,R,DTilde,Integrand)

  if (mode>0) then
    GradL = 0.0d0
  end if
end subroutine Like1C


!--------------------------------------
!  CASE 2:  0 < d1 < K
!--------------------------------------
subroutine Like1B(iHH,d1,d2,d3,parms,mode,L,GradL)
  use nrtype
  use GlobalModule, only : ParmsStructure,HHData,IntRule
  use NewTools,     only : ComputeSVD,ComputeInverse_LU,ComputeLQ
  use nag_library,  only : F07FDF,F04BAF
  ! F07FDF = Cholesky decomposition
  ! F04BAF = matrix inversion using LU decomposition

  implicit none
  integer(i4b),         intent(in)    :: iHH,d1,d2,d3
  type(ParmsStructure), intent(in)    :: parms
  integer(i4b),         intent(inout) :: mode 
  real(dp),             intent(out)   :: L
  real(dp),             intent(inout) :: GradL(:)

  ! F07FDF,F04BAF variables
  character(1)              :: UPLO
  real(dp)                  :: rCond,ErrBound
  integer(i4b)              :: ifail
  integer(i4b), allocatable :: ipivot(:)
  real(dp),     allocatable :: temp2(:,:)

  integer(i4b)              :: i1,i2,ix
  integer(i4b), allocatable :: index1(:),index2(:)
  real(dp),     allocatable :: B1(:,:)
  real(dp),     allocatable :: U(:,:),S(:),VT(:,:)
  real(dp),     allocatable :: B2Tilde(:,:),B21Tilde(:,:),B22Tilde(:,:)
  real(dp),     allocatable :: M1(:,:),DTilde(:)
  real(dp),     allocatable :: Temp1(:)         ! used to compute DTilde and epsilon
  real(dp),     allocatable :: C(:,:),Omega(:,:)
  real(dp),     allocatable :: Omega11(:,:),Omega22(:,:),Omega12(:,:),C2(:,:)
  real(dp),     allocatable :: nu(:),epsilon1(:)
  real(dp),     allocatable :: Psi(:,:),CPsi(:,:)
  real(dp),     allocatable :: R(:,:),Q(:,:)
  integer(i4b), allocatable :: RowGroup(:)
  integer(i4b)              :: integrand

  !------------------------------------------------------------------------  
  ! Case 2:   d1 < K products were purchased
  !           (HHData%nNonZero(iHH) == d1 and d1 < K)
  !           Mapping from q(iNonZero) to e is NOT one-to-one
  !           need to integrate across region of e-space
  !           satisfying M1*z2<=DTilde
  !------------------------------------------------------------------------
  !  d1 = dim(q>0)
  !  d2 = K - d1
  !  d3 = J - d1
  allocate(index1(d1))
  allocate(index2(d2))
  
  index1 = (/(ix,ix=1,d1)/)
  index2 = (/(ix,ix=d1+1,parms%K)/)

  ! size(B1) = (K x d1)
  allocate(B1(parms%K,d1))
  B1 = parms%B(:,HHData%iNonZero(index1,iHH))

  ! Compute SVD of B1 where B1 is (K x d1)
  ! B1 = U*S*VT
  ! S is (d1 x d1) diagonal and so is stored as a vector
  allocate(U(parms%K,parms%K))
  allocate(S(d1))
  allocate(VT(d1,d1))

  call ComputeSVD(B1,U,S,VT)

  ! size(B2) = (K x d3)
  ! size(B21) = (d1 x d3)
  ! size(B22) = (d2 x d3)

  allocate(B21Tilde(d1,d3),B22Tilde(d2,d3))
  allocate(B2Tilde(parms%K,d3))

  B2Tilde  = matmul(transpose(U),parms%B(:,HHData%iZero(1:d3,iHH)))
  B21Tilde = B2Tilde(index1,:)
  B22Tilde = B2Tilde(index2,:)

  ! Compute nu = U.' * MuE
  allocate(nu(parms%K))
  nu = matmul(transpose(U),parms%MuE+parms%mue_month(:,HHData%month(ihh)))

  ! Omega = U' * C * C' * U
  allocate(Omega(parms%K,parms%K))
  Omega = matmul(parms%CSig,transpose(parms%CSig))
  Omega = matmul(Omega,U)
  Omega = matmul(transpose(U),Omega)

  allocate(Omega11(d1,d1),Omega22(d2,d2),Omega12(d1,d2),C2(d2,d2))

  Omega11 = Omega(index1,index1)
  Omega22 = Omega(index2,index2)
  Omega12 = Omega(index1,index2)

  ! Cholesky decomposition of Omega22: lower triangular form
  !      C2*C2' = Omega22
  C2 = Omega22
  UPLO = 'L'
  call F07FDF(UPLO,d2,C2,d2,ifail)

  allocate(Psi(d1,d1))

  ! Compute  Psi = Omega11 - Omega12*inv(Omega22)*Omega12'
  Psi = 0.0d0
  allocate(iPivot(d2))
  allocate(temp2(d2,d1))
  temp2 = transpose(Omega12)
  call F04BAF(d2,d1,Omega22,d2,iPivot,temp2,d2,rCond,errBound,ifail)
  !  Psi = Omega11 - Omega12*inv(Omega22)*Omega12'
  !  size( inv(Omega22)*Omega12' ) = (d2 x d1)
  Psi = Omega11 - matmul(Omega12,temp2)
  deallocate(temp2,iPivot)
   
  allocate(CPsi(d1,d1))
  ! Cholesky decomposition of Psi: lower triangular form
  !      CPsi * CPsi' = Psi
  CPsi = Psi
  UPLO = 'L'
  call F07FDF(UPLO,d1,CPsi,d1,ifail)
  if (ifail>0) then
    print *, 'Psi is not positive definite'
    print *, 'Psi'
    do i2=1,d1
      print *, Psi(i2,:)
    end do
    print *,Omega11
    do i2=1,d1
      print *,Omega11(i2,:)
    end do
    stop
  end if

  allocate(M1(d3,d2))
  allocate(DTilde(d3))
  allocate(epsilon1(d1))
  allocate(Temp1(d1))

  ! size(M1)       = (d3 x d2)
  ! M1*z2 <= DTilde
  ! size(B21Tilde) = (d1 x d3)
  ! size(S)        = (d1 x 1)
  ! size(VT)       = (d1 x d1)
  ! size(DTilde)   = (d3 x 1)
  M1    = matmul(transpose(B22Tilde),C2)
  ! extract p1 = price(iNonZero)
  Temp1 = matmul(VT,HHData%p(HHData%iNonZero(index1,iHH),iHH))

  ! Compute inv(S1)*VT*p1
  !   epsilon1 = inv(S1)*VT*p1 + S1*VT*q1
  !   DTilde   = p2 - B21Tilde.'*inv(S1)*VT*p1 - B22Tilde.'*nu2
  Temp1 = Temp1/S(index1)
  epsilon1 = Temp1
  DTilde  = HHData%p(HHData%iZero(1:d3,iHH),iHH)               &
           - matmul(transpose(B21Tilde),Temp1)

  ! extract q1 = q(index1)
  ! nonzero elements of q are stored in q(index1)
  Temp1 = matmul(VT,HHData%q(index1,iHH))
  Temp1 = Temp1*S(index1)

  epsilon1 = epsilon1 + Temp1
  DTilde = DTilde - matmul(transpose(B22Tilde),nu(index2))
  deallocate(Temp1)

  ! M1 = R*Q  : LQ Decomposition of M1
  !             R = lower triangular
  !             Q = orthogonal
  ! size(M1) = (d3 x d2)
  ! size(R)  = (d3 x d2)
  ! size(Q)  = (d2 x d2)

  allocate(R(d3,d2),Q(d2,d2))
  call ComputeLQ(d3,d2,M1,R,Q,ifail)

  ! L = integral of DensityFunc over region of x satisfying R*x<=DTilde
  !
  Integrand = 2
  
  ! size(RowGroup) = d3
  ! RowGroup(j1) = max(i1) s.t.  LL(i1,j1) .ne. 0.0d0

  allocate(RowGroup(d3))
  call ComputeRowGroup(R,d3,d2,RowGroup)

  call ComputeProb(L,IntRule%rule(d2)%nodes,IntRule%rule(d2)%weights, &
                   RowGroup,R,DTilde,Integrand,           &
                   epsilon1,nu(index1),Omega12,C2,CPsi,Q,S(index1),     &
                   d1,d2)

  deallocate(B1,index1,index2)
  deallocate(U,S,VT)
  deallocate(nu)
  deallocate(B21Tilde,B22Tilde,B2Tilde)
  deallocate(M1,DTilde)
  deallocate(Omega)
  deallocate(Omega11,Omega22,Omega12,C2)
  deallocate(Psi,CPsi)
  deallocate(R,Q)
  deallocate(RowGroup)

  if (mode>0) then
    GradL = 0.0d0
  end if

end subroutine Like1B

!--------------------------------------
! CASE 1:  d1==K
!--------------------------------------
subroutine Like1A(i1,parms,mode,L,GradL)
  use nrtype
  use GlobalModule, only : ParmsStructure,HHData,inf,ReadWriteParameters
  use nag_library, only : F04BAF,F07AAF,F07ADF,F03BAF
  ! F04BAF = matrix inversion using LU decomp
  ! F07ADF = LU decomposition
  ! F03BAF = determinant of matrix after LU decomp

  implicit none
  integer(i4b),         intent(in)    :: i1
  type(ParmsStructure), intent(inout) :: parms
  integer(i4b),         intent(inout) :: mode
  real(dp),             intent(out)   :: L
  real(dp),             intent(inout) :: GradL(:)

  real(dP)                  :: crit

  real(dp), allocatable     :: B1(:,:),B1T(:,:)
  real(dp), allocatable     :: B2(:,:)
  real(dp), allocatable     :: e(:)
  real(dp), allocatable     :: mue_month(:)
  real(dp)                  :: F
  real(dp), allocatable     :: GradF_e(:),GradF_MuE(:)
  real(dp), allocatable     :: GradF_InvC(:,:)

  ! variables used by F04BAF,F07ADF,F03BAF
  real(dp)                  :: rcond,errBound 
  integer(i4b), allocatable :: iPivot(:)
  integer(i4b)              :: ifail
  real(dp)                  :: D
  integer(i4b)              :: ID
  real(dp), allocatable     :: temp1(:)
  integer(i4b)              :: i2

  allocate(B1(parms%K,parms%K),B1T(parms%K,parms%K))
  allocate(B2(parms%K,parms%J-parms%K))
  allocate(e(parms%K))
  allocate(mue_month(parms%K))
  allocate(GradF_e(parms%K),GradF_MuE(parms%K))
  allocate(GradF_InvC(parms%K,parms%K))
  allocate(iPivot(parms%K))
  allocate(temp1(parms%K))

  B1 = parms%B(:,HHData%iNonZero(1:parms%K,i1))
  B2 = parms%B(:,HHData%iZero(1:parms%J-parms%K,i1))

  ! Compute  e = inv(B1')*p1 + B1*q1
  e = HHData%p(HHData%iNonZero(1:parms%K,i1),i1)
  B1T = transpose(B1)
  ifail = -1
  ! B1T is changed by this routine, but don't need B1T after this
  !call F04BAF(parms%K,1,B1T,parms%K,iPivot,e,parms%K,rCond,errBound,ifail)
  call F07AAF(parms%K,1,B1T,parms%K,iPivot,e,parms%K,ifail)
  if (ifail .ne. 0 .and. ifail .ne. parms%K+1) then
    print *,'ifail = ',ifail,'. B1 is singular'
    call ReadWriteParameters(parms,'write')
    do i2=1,parms%K
      print *,B1(i2,:)
    end do
    stop
  end if

  ! crit = min(  p(iZero) - B2' * (inv(B1') * p(iNonZero)
  crit = minval(HHData%p(HHData%iZero(1:parms%J-parms%K,i1),i1) &
       - matmul(transpose(B2),e))

  if (crit<0.0d0) then
    L = 0.0d0
    if (mode>0) then
      GradL = 0.0d0
    end if
  else
    e = e + matmul(B1,HHData%q(1:parms%K,i1))
   
    ! mue_month = seasonal adjustment to mue
    mue_month = parms%mue_month(:,HHData%month(i1))
    call LogDensityFunc(e,parms,mue_month, &
                        F,GradF_e,GradF_MuE,GradF_InvC)

    ! F07ADF:  compute LU factorization of B1
    ! B1 is changed
    call F07ADF(parms%K,parms%K,B1,parms%K,iPivot,ifail)

    ! F03BAF:  compute determinant of B1 after factorizing
    call F03BAF(parms%k,B1,parms%K,iPivot,D,ID,ifail)

    ! log(det(B1)) = log(D) + ID*log(2.0)
    L = exp(F + log(abs(D)) + real(ID)*log(2.0d0))

    if (mode>0) then
      GradL = 0.0d0
    end if
  end if
  deallocate(B1,B1T)
  deallocate(B2)
  deallocate(iPivot)
  deallocate(GradF_e,GradF_MuE)
  deallocate(GradF_InvC)
  deallocate(e)
  deallocate(mue_month)
  deallocate(temp1)
end subroutine Like1A


subroutine Like1(mode,nx,x,L,GradL,nstate,iuser,ruser)
  use nrtype
  use GlobalModule, only : parms,HHData,iFree
  implicit none

  integer(i4b), intent(inout)        :: mode
  integer(i4b), intent(in)           :: nx
  real(dp),     intent(in)           :: x(:)
  real(dp),     intent(out)          :: L
  real(dp),     intent(inout)        :: GradL(:)
  integer(i4b), intent(in)           :: nstate
  integer(i4b), intent(in)           :: iuser(*)
  real(dp),     intent(in)           :: ruser(*)

  real(dp)                           :: L1
  real(dp), allocatable              :: GradL1(:)

! Compute likelihood and its gradient w.r.t. x where utility is
!
!     u  = y - p'*q - 0.5 * (B*q - e)' * (B*q - e)
!
! 
! HHData%N        = (1 x 1)  number of households
!
! Revision history
! 2015JUN07 LN  move all actual work to (Like1A,Like1B,Like1C)
! 22MAR2013 LN  continue translations from matlab code
! 11DEC2012 LN  continue translation from matlab code
! 07DEC2012 LN  Start translation from matlab code: Like1.m
!
real(dp)                  :: small
integer(i4b)              :: iHH,d1,d2,d3

call UpdateParms(x,iFree,parms)

! Initial values for likelihood and gradient
L     = 0.0d0
if (mode>0) then
  GradL = 0.0d0
end if

allocate(GradL1(nx))

small=1e-6   ! avoid log of zero

! Loop through households
do iHH=1,HHData%N

  d1 = HHData%nNonZero(iHH)
  d2 = parms%K - d1
  d3 = parms%J - d1

  call Like1_wrapper(iHH,d1,d2,d3,parms,mode,L1,GradL1)
  L  = L + log(L1+small)
  if (mode>0) then
    GradL = GradL + GradL1/(L1+small)
  end if
end do        ! do iHH=1,HHData%N

L = -L
if (mode>0) then
  GradL = -GradL
end if

deallocate(GradL1)

end subroutine Like1

subroutine UpdateParms(xFree,iFree,parms)
  use nrtype
  use GlobalModule, only : SelectFreeType,ParmsStructure
  use NewTools,     only : SphereToMatrix,ComputeInverse_LU
  implicit none

! Revision history
! 09DEC2012 LN  translated from matlab file Updateparms.m
!
! Utility parameters
! B   = matrix of utility parameters
!     = SphereToMatrix(C,D,K,J)
!
! C   = K*(K-1)/2 + (K-1)*(J-K) spherical coordinate representation of B
!       stored as a vector
! D   = J x 1 scaling factors 
!
real(dp),             intent(in)    :: xFree(:)
type(SelectFreeType), intent(in)    :: iFree
type(ParmsStructure), intent(inout) :: parms

integer(i4b)                        :: ifail

! iFree%D   elements of xFree corresponding to D
! iFree%xD  elements of D that are free
if (size(iFree%D)>0) then
  parms%D(iFree%D) = xFree(iFree%xD)
end if

! iFree%BC  elements of xFree corresponnding to BC
! iFree%xBC elements of BC that are free
! BC has at most K*(K-1)/2 + (K-1)*(J-K) free elements
! parms%BC is a vector storing the normalized columns of B
!           each column stored in BC is normalized to have norm 1.
if (size(iFree%BC)>0) then
  parms%BC(iFree%BC) = xFree(iFree%xBC)
end if

! parms%B         = (K x J) B matrix
! parms%GradBC = (K x n) gradient of B w.r.t. C
! parms%GradBD   = (K x J) gradient of B w.r.t. D
call SphereToMatrix(parms%BC,parms%D,parms%K,parms%J,parms%B, &
                    parms%GradBC,parms%GradBD)

! MuE
! iFree%MuE  = elements of xFree corresponding to MuE
! iFree%xmue = elements of MuE that are free
if (size(iFree%MuE>0)) then
  parms%MuE(iFree%MuE) = xFree(iFree%xmue)    
end if

! iFree%InvCDiag  = elements of xFree corresponding to InvCDiag
! iFree%xInvCDiag = elements of InvCDiag that are free
if (size(iFree%InvCDiag)>0) then
  parms%InvCDiag(iFree%InvCDiag) = xFree(iFree%xInvCDiag)
end if

! iFree%InvCOffDiag  elements of xFree correspondning to InvCOffDiag
! iFree%xInvCOffDiag elements of InvCOffDiag that are free
! InvCOffDiag has at most K*(K-1)/2 free elements
! parms%InvCOffDiag is a vector storing the normalized columns of InvC
!           each column stored in InvC_phi is normalized to have norm 1.
if (size(iFree%InvCOffDiag)>0) then
  parms%InvCOffDiag(iFree%InvCOffDiag) = xFree(iFree%xInvCOffDiag)
end if

! parms%InvC         = (K x K) upper triangular Cholesky decomposition 
!                                 of inverse of Sig
!                                 inv(sig) = InvC'*InvC
! parms%GradInvCOffDiag = (K x n) gradient of InvC w.r.t. InvCOffDiag
! parms%GradInvCDiag    = (K x K) gradient of InvC w.r.t. InvCDiag
if (iFree%nInvCDiag>0 .or. iFree%nInvCOffDiag>0) then
  call SphereToMatrix(parms%InvCOffDiag,parms%InvCDiag,parms%K,parms%K, &
                      parms%InvC,parms%GradInvCOffDiag,parms%GradInvCDiag)
  parms%InvC = transpose(parms%InvC)
  ifail = 0
  call ComputeInverse_LU(parms%K,parms%InvC,parms%CSig,ifail)
end if

end subroutine UpdateParms

!-----------------------------------------------------------------------
!
! UpdateParms2:   update parameters for like2
!
! Revision history
!
! 2018MAY23 LN  add seasonality to model 
! 10JUN2014 LN  adapted from matlab file UpdateParms_RC.m
!
! Utility parameters
!     log(BD) = BD_z * BD_beta + BD_C * eta + BD_month(month)
!     BC      = pi * normcdf(Y)
!     Y      = BC_z * BC_beta + BC_C * eta
! 
! BD_beta = (nzD x 1)
! BD_C    = (J x dim_eta)
! BC_beta = (nzC x 1)
! BC_C    = (nBC  x dim_eta)
!
!-----------------------------------------------------------------------
subroutine UpdateParms2(x,iFree,parms)
  use nrtype
  use GlobalModule, only : SelectFreeType,ParmsStructure
  use NewTools,     only : SphereToMatrix,ComputeInverse_LU
  implicit none

  real(dp),             intent(in)    :: x(:)
  type(SelectFreeType), intent(in)    :: iFree
  type(ParmsStructure), intent(inout) :: parms

  integer(i4b)                        :: ifail
  real(dp), allocatable               :: temp1(:),temp2(:,:)

  ! log(BD)      = BD_beta * BD_z + BD_C * BD_eta
  ! size(BD_eta) = parms%dim_eta
  ! size(BD_C)   = (J x dim_eta)

  ! BD_beta
  ! iFree%BD_beta   free elements of BD_beta
  ! iFree%xBD_beta  elements of x corresponding to BD_beta
  if (iFree%nBD_beta>0) then
    parms%BD_beta(iFree%BD_beta) = x(iFree%xBD_beta)
  end if

  ! BD_CDiag : diagonal elements of BD_C
  if (iFree%nBD_CDiag>0) then
    parms%BD_CDiag(iFree%BD_CDiag) = x(iFree%xBD_CDiag)
  end if

  ! BD_COffDiag : off-diagonal elements of BD_C
  if (iFree%nBD_COffDiag>0) then
    parms%BD_COffDiag(iFree%BD_COffDiag) = x(iFree%xBD_COffDiag)
  end if

  ! BD_month: month dummy
  if (iFree%nBD_month>0) then
    allocate(temp1(11*parms%J))
    temp1 = reshape(parms%BD_month(:,2:12),(/11*parms%J/))
    temp1(ifree%bd_month) = x(ifree%xbd_month)
    parms%BD_month(:,2:12) = reshape(temp1,(/parms%J,11/))
    deallocate(temp1)
  end if

  if (iFree%nBD_CDiag>0 .or. iFree%nBD_COffDiag>0) then
    ! Compute parms%BD_C
    ! size(BD_C) = (J x dim_eta)
    ! SphereToMatrix creates upper triagular matrix of size
    !    (dim_eta x J)
    ! BD_C is lower triangular:  we need the transpose as below.
    allocate(temp2(parms%dim_eta,parms%J))
    call SphereToMatrix(parms%BD_COffDiag,parms%BD_CDiag,  &
                        parms%dim_eta,parms%J,temp2,     &
                        parms%GradBD_C_COffDiag,parms%GradBD_C_CDiag)
    parms%BD_C = transpose(temp2)
    deallocate(temp2)
  end if

  ! norminv(BC) = BC_beta * BC_z + BC_C * eta
  ! size(BC)    = (nBC x dim_eta)

  ! iFree%BC_beta   free elements of BC_beta
  ! iFree%xBC_beta  elements of x corresponding to BC_beta
  if (iFree%nBC_beta>0) then
    parms%BC_beta(iFree%BC_beta) = x(iFree%xBC_beta)
  end if

  ! BC_CDiag : diagonal elements of BC_C
  if (iFree%nBC_CDiag>0) then
    parms%BC_CDiag(iFree%BC_CDiag) = x(iFree%xBC_CDiag)
  end if

  ! BC_COffDiag : off-diagonal elements of BC_C
  if (iFree%nBC_COffDiag>0) then
    parms%BC_COffDiag(iFree%BC_COffDiag) = x(iFree%xBC_COffDiag)
  end if

  if (iFree%nBC_CDiag>0 .or. iFree%nBC_COffDiag>0) then
    ! Compute parms%BC_C
    ! size(BC_C) = (nBC x dim_eta)
    ! SphereToMatrix creates upper triangular matrix of size
    !   (dim_eta x nBC)
    ! after call, compute transpose to create lower triangular matrix
    allocate(temp2(parms%dim_eta,parms%nBC))
    call SphereToMatrix(parms%BC_COffDiag,parms%BC_CDiag,  &
                        parms%dim_eta,parms%nBC,temp2,   &
                        parms%GradBC_C_COffDiag,parms%GradBC_C_CDiag)
    parms%BC_C = transpose(temp2)
    deallocate(temp2)
  end if

  ! MuE
  ! iFree%xmue  = elements of xFree corresponding to MuE
  ! iFree%mue = elements of MuE that are free
  if (iFree%nMuE>0) then
    parms%MuE(iFree%MuE) = x(iFree%xmue)    
  end if

  ! MuE_month
  ! Month 1 is baseline
  ! iFree%mue_month = elements of mue_month that are free
  ! iFree%xmue_month  = elements of xFree corresponding to MuE_month
  if (iFree%nmue_month>0) then
    allocate(temp1(11*parms%K))
    temp1 = reshape(parms%mue_month(:,2:12),(/parms%K*11/))
    temp1(iFree%mue_month) = x(iFree%xmue_month)
    parms%mue_month(:,2:12) = reshape(temp1,(/parms%K,11/))
    deallocate(temp1)
  end if

  ! sig = covariance matrix of e
  ! sig = C * C^T
  ! inv(C) = InvC

  ! iFree%InvCDiag  = elements of InvCDiag that are free
  ! iFree%xInvCDiag = elements of x corresponding to InvCDiag
  if (iFree%nInvCDiag>0) then
    parms%InvCDiag(iFree%InvCDiag) = x(iFree%xInvCDiag)
  end if

  ! iFree%InvCOffDiag  = elements of InvCOffDiag that are free
  ! iFree%xInvCOffDiag = elements of x corresponding to InvCOffDiag
  if (iFree%nInvCOffDiag>0) then
    parms%InvCOffDiag(iFree%InvCOffDiag) = x(iFree%xInvCOffDiag)
  end if

  ! parms%InvC         = (K x K) inverse of lower triangular Cholesky decomposition
  !                              of Sig
  !                              sig = C * C'
  !                              inv(sig) = inv(C') * inv(C)
  !                                       = InvC'   * InvC
  ! parms%GradInvCOffDiag = (K x n) gradient of InvC w.r.t. InvCOffDiag
  ! parms%GradInvCDiag    = (K x K) gradient of InvC w.r.t. InvCDiag
  ! SphereToMatrix creates an upper triangular matrix
  ! must take transpose after creation
  if (iFree%nInvCDiag>0 .or. iFree%nInvCOffDiag>0) then
    call SphereToMatrix(parms%InvCOffDiag,parms%InvCDiag,parms%K,parms%K, &
                        parms%InvC,parms%GradInvCOffDiag,parms%GradInvCDiag)
    ! Make InvC a lower triangular matrix
    parms%InvC = transpose(parms%InvC)

    ! need C as well
    ifail = 0
    call ComputeInverse_LU(parms%K,parms%InvC,parms%CSig,ifail)
  end if

  !print *, 'Warning: gradients of SphereToMatrix need to be adjusted for transpose.'
end subroutine UpdateParms2

subroutine LogDensityFunc(e,parms,mue_month,F,GradF_e,GradF_MuE,GradF_InvC)
! e          = (K x 1)
! parms      = parameter structure
! mue_month  = (K x 1) seasonal adjustment
! F          = (1 x 1) log density evaluated at e
! GradF_e    = (K x 1) gradient of log density w.r.t. e
! GradF_MuE  = (K x 1) gradient of log density w.r.t. MuE
! GradF_InvC = (K x K) gradient of log density w.r.t. InvC 
!                           
! Revision history
! 09DEC2012 LN translated to Fortran from matlab LogDensity.m
  use nrtype
  use nag_library,  only : F04BAF
  use GlobalModule, only : ParmsStructure
  implicit none
  real(dp),             intent(in)  :: e(:)
  type(ParmsStructure), intent(in)  :: parms
  real(dp),             intent(in)  :: mue_month(:)
  real(dp),             intent(out) :: F
  real(dp),             intent(out) :: GradF_e(:)
  real(dp),             intent(out) :: GradF_MuE(:)
  real(dp),             intent(out) :: GradF_InvC(:,:)

  real(dp), allocatable    :: z(:),TempInvC(:,:)
  integer(i4b)             :: i1,i2
  real(dp), allocatable    :: temp1(:),temp2(:,:)
  real(dp)                 :: DetInvC

  ! workspace to compute determinant and matrix inverse
  integer(i4b)              :: ifail
  integer(i4b), allocatable :: iPivot(:)
  real(dp)                  :: rCond,errBound

  !external F04BAF  ! compute inverse of matrix using LU decomp

  allocate(iPivot(parms%K))
  allocate(z(parms%K),TempInvC(parms%K,parms%K))
  allocate(temp1(parms%K))

  z = matmul(parms%InvC,e - parms%MuE-mue_month)

  ! determinant of InvC
  DetInvC = 1.0d0
  do i1=1,parms%K
    DetInvC = DetInvC * parms%InvC(i1,i1)
  end do  

  ! Ln(L) = -0.5*z'*z -0.5*K*log(2*pi) + log(det(InvC))
  F     = -0.5d0*dot_product(z,z) - 0.5d0*real(parms%K,dp)*log(2.0d0*pi) &
          + log(DetInvC)

  GradF_e    = -matmul(transpose(parms%InvC),z)
  GradF_MuE  = -GradF_e
  GradF_InvC = 0.0d0
  allocate(temp2(parms%K,parms%K))
  do i1=1,parms%K
  do i2=1,i1
    temp2 = 0.0d0
    temp2(:,i2) = temp2(:,i2) + parms%InvC(i1,:)
    temp2(i1,:) = temp2(i1,:) + parms%InvC(i2,:)
    temp1       = matmul(temp2,e-parms%MuE-mue_month)
    GradF_InvC(i2,i1) = -0.5d0*dot_product(e-parms%MuE-mue_month,temp1)
  end do
  end do
  
  temp2 = 0.0d0
  do i1=1,parms%K
    temp2(i1,i1) = 1.0d0
  end do
  TempInvC = parms%InvC
  ! compute temp2 = inv(parms%InvC)
  call F04BAF(parms%K,parms%K,TempInvC,parms%K,iPivot,temp2, &
              parms%K,rCond,errBound,ifail)

  GradF_InvC = GradF_InvC + transpose(temp2)

  deallocate(temp2)
  deallocate(temp1)
  deallocate(z,TempInvC)
  deallocate(iPivot)
end subroutine LogDensityFunc

subroutine ComputeRowGroup(L,J,d,RowGroup)
  use nrtype
  implicit none
  integer(i4b), intent(in)  :: J,d
  real(dp), intent(in)      :: L(J,d)
  integer(i4b), intent(out) :: RowGroup(J)

  integer(i4b) :: j1,ix 
  ! for each row find the right-most column with a non-zero element
  ! RowGroup(j)==i1 indicates that L(j,i1)~=0 and L(j,i1+1:d)==0
  RowGroup = 0
  do j1=1,J
    RowGroup(j1) =  maxval(pack((/(ix,ix=1,d)/),L(j1,1:min(j1,d)) .ne. 0.0d0))
  end do
end subroutine ComputeRowGroup

subroutine ComputeProb(p,v,w,RowGroup,R,DTilde,Integrand,  &
                       epsilon1,nu1,Omega12,C2,CPsi,Q,S1,  &
                       d1,d2)
  
  use nrtype
  use nag_library, only : S15ABF
  use ToolsModule, only : InverseNormal_mkl
  use GlobalModule, only : DensityGradType   ! Gradient of F w.r.t. parameters

  implicit none
  ! Compute integral of f(x) * exp(-0.5*x'*x)/((2*pi)^d/2) over the region
  ! R*x <= DTilde.
  !
  ! v        = (n x d) integration nodes lieing in [-1,1]^d  
  ! w        = (n x 1) w(i1) integration weights for node x(i1,:)
  ! RowGroup = (J x 1) index of group that row j belongs to
  !                     if RowGroup(j)==i1, then row(j) imposes 
  !                    a constraint on x(i1) because R(j,i1)~=0
  !                    and R(j,i1+1:d)==0
  !                    a column in row(j) is i1.
  ! R        = (J x d2) lower triangular matrix of constraints
  ! DTilde   = (J x 1) right side of constraints
  ! func     = (function handle) integrand 
  real(dp),     intent(out) :: p
  real(dp),     intent(in)  :: v(:,:)
  real(dp),     intent(in)  :: w(:)
  integer(i4b), intent(in)  :: RowGroup(:)
  real(dp),     intent(in)  :: R(:,:)
  real(dp),     intent(in)  :: DTilde(:)
  integer(i4b), intent(in)  :: Integrand 
  real(dp),     optional, intent(in) :: epsilon1(:),nu1(:),Omega12(:,:),C2(:,:),CPsi(:,:)
  real(dp),     optional, intent(in) :: S1(:),Q(:,:)
  integer(i4b), optional, intent(in) :: d1,d2

  real(dp)                  :: small,big
  real(dp), allocatable     :: u(:,:),xT(:,:),J(:),F(:)
  integer(i4b)              :: n,QuadDim
  logical                   :: ZeroProb
  integer(i4b)              :: i1,i2,i3
  logical, allocatable      :: Rows_Pos(:),Rows_Neg(:)
  real(dp), allocatable     :: ULO(:),UHI(:) 
  integer(i4b)              :: n1
  real(dp), allocatable     :: t1(:),t2(:),t3(:)
  real(dp), allocatable     :: R1(:)
  real(dp)                  :: xtemp
  integer(i4b)              :: ifail

  ! gradient of F w.r.t. parameters
  type(DensityGradType) :: GradF

  small = tiny(1.0d0)   ! smallest real number
  big   = huge(1.0d0)   ! biggest real number
  n     = size(v,1)
  QuadDim     = size(v,2)
  allocate(u(n,QuadDim),xT(n,QuadDim),J(n),F(n))
  allocate(ULO(n),UHI(n))

  ! First map v into [0,1]
  u = 0.5d0*(v+1.0d0)

  ! Then map u into xT
  xT = 0.0d0

  ! J = Jacobian of mapping from v to u
  J = 1.0d0

  ! Loop through dimensions of v
  ZeroProb= .false.

  allocate(Rows_Pos(size(R,1)),Rows_Neg(size(R,1)))
  do i1=1,QuadDim
    Rows_Pos = (RowGroup .eq. i1) .and. (R(:,i1)>0)
    Rows_Neg = (RowGroup .eq. i1) .and. (R(:,i1)<0)
 
    if (i1==1) then
      if (any(Rows_Neg)) then
        xtemp = maxval(pack(DTilde,Rows_Neg)/pack(R(:,i1),Rows_Neg))
        ! S15ABF(x,ifail) = normcdf of x
        ifail = 0
        ULO = S15ABF(xtemp,ifail)
      else
        ULo = 0.0d0
      end if

      if (any(Rows_Pos)) then
        xtemp = minval(pack(DTilde,Rows_Pos)/pack(R(:,i1),Rows_Pos))
        ifail = 0
        UHI = S15ABF(xtemp,ifail)
      else
        UHi = 1.0d0
      end if
      if (all(UHi<=ULo+small)) then
        ZeroProb= .true.
        return
      end if
      u(:,1) = (UHi - ULo)*u(:,1)
      ifail = 0
      xT(:,i1) = max(InverseNormal_mkl(u(:,1),n,ifail),-big)

      J = J*(UHi-ULo)/2.0d0 
    elseif (i1>1) then 
      if (any(Rows_Neg)) then
        n1 = count(Rows_Neg)
        allocate(t1(n1),t2(n1),t3(n1),R1(n1))
        R1 = pack(R(:,i1),Rows_Neg)
        t1 = pack(DTilde,Rows_Neg)/R1
        do i3=1,i1-1
          t2 = pack(R(:,i3),Rows_Neg)/R1
          t3 = t1
          do i2=1,n 
            t3=t3-t2*xT(i2,i3)
            ! S15ABF normal CDF
            ifail = 0
            ULO(i2) = min(S15ABF(maxval(t3),ifail),1.0d0-small)
          end do
        end do
        deallocate(t1,t2,t3,R1)
      else
        ULo = 0.0d0
      end if

      if (any(Rows_Pos)) then
        ! size(t1) = nPos x n
        ! size(t2) = nPos x i1-1
        n1 = count(Rows_Pos)
        allocate(R1(n1),t1(n1),t2(n1),t3(n1))
        R1 = pack(R(:,i1),Rows_Pos)
        t1 = pack(DTilde,Rows_Pos)/R1
        do i3=1,i1-1
          t2 = pack(R(:,i3),Rows_Pos)/R1
          t3=t1
          do i2=1,n
            t3 = t3 - t2*xT(i2,i3)
            ! S15ABF = normal CDF
            ifail = 0
            UHI(i2) = max(S15ABF(maxval(t3),ifail),ULO(i2)+small)
          end do
        end do
        deallocate(R1,t1,t2,t3)
      else
        UHi = 1.0d0
      end if
      u(:,i1) = (UHi-ULo)*u(:,i1)
      ! xT(:,i1) = InverseNormal of u(:,i1)
      xT(:,i1) = max(InverseNormal_mkl(u(:,i1),n,ifail),-big)

      J = 0.5d0*J*(UHi-ULo)     
    end if ! if i1==1
  end do   ! do i1=1,d

  if (ZeroProb) then
    p=0.0d0
  else
    if (Integrand==2) then
      ! allocate memory for GradF
      call AllocateGradF(GradF,'allocate',n,d1,d2)

      ! compute F and GradF
      call DensityFunc2(xT,epsilon1,nu1,Omega12,C2,CPsi,Q,S1,d1,d2,F,GradF)

      ! deallocate memory for GradF
      call AllocateGradF(GradF,'deallocate')

    elseif (Integrand==3) then
      F = 1.0d0
    end if
    
    p = dot_product(F*J,w)
  end if
  deallocate(u,xT,J,F)
  deallocate(Rows_Pos,Rows_Neg) 
  deallocate(ULO,UHI)
end subroutine ComputeProb

subroutine AllocateGradF(GradF,action,nx,d1,d2)
use nrtype
use GlobalModule, only : DensityGradType
implicit none
type(DensityGradType), intent(inout) :: GradF
character(len=*),      intent(in)    :: action
integer(i4b), optional,intent(in)    :: nx,d1,d2

select case (action)

case ('allocate')
  allocate(GradF%S1(nx,d1))
  allocate(GradF%nu1(nx,d1))
  allocate(GradF%Omega12(nx,d1,d2))
  allocate(GradF%VT(nx,d1,d1))
  allocate(GradF%C2(nx,d2,d2))
  allocate(GradF%Q(nx,d2,d2))
  allocate(GradF%CPsi(nx,d1,d1))
  GradF%S1 = 0.0d0
  GradF%nu1 = 0.0d0
  GradF%Omega12 = 0.0d0
  GradF%VT = 0.0d0
  GradF%C2 =0.0d0
  GradF%Q = 0.0d0
  GradF%CPsi = 0.0d0
case ('deallocate')
  deallocate(GradF%S1)
  deallocate(GradF%nu1)
  deallocate(GradF%Omega12)
  deallocate(GradF%VT)
  deallocate(GradF%C2)
  deallocate(GradF%Q)
  deallocate(GradF%CPsi)

end select 
end subroutine AllocateGradF

subroutine DensityFunc2(xT,epsilon1,nu1,Omega12,C2,CPsi,Q,S1,d1,d2,F,GradF)
  use nrtype
  use GlobalModule, only : DensityGradType
  use ToolsModule, only : InvertLower
!  use NewTools, only : det   ! compute determinant of matrix
  implicit none
  ! Case 2: 0 < d1 < K
  ! density of data to be integrated
  !
  ! Compute f(q1,x) = f(epsilon1(q1),x) 
  !                   epsilon2 = C2*z2+nu2
  !                   z2 = x*Q     
  !                   v1 = nu1 + Omega12*(C2.'\z2)
  !
  ! F        = (nx x 1)   output, one value for each value of x
  ! xT       = (nx x d2)  variable of integration
  ! epsilon1 = (d1 x 1)   value of epsilon1
  !          =            inv(S1)* VT * p1 + S1 * VT * q1
  ! nu1      = (d1 x 1)   mean of epsilon1
  ! Omega12  = (d1 x d2)  covariance of (epsilon1,epsilon2)
  ! C2       = (d2 x d2)  Cholesky decomposition of Omega22, lower triangular
  ! CPsi     = (d1 x d1)  Cholesky decomposition of variance of epsilon1 conditional on epsilon2
  ! Q        = (d2 x d2)  orthogonal so that inv(Q) = Q.';
  ! S1       = (d1 x 1)   non-zero elements of S where U*S*VT = B1
  ! d1         = (1 x 1)    size of e1 and q1
  ! d2         = (1 x 1)    size of e2
  !
  ! Revision history
  ! 2015JUN07 LN major overhaul reflecting better DensityFunc.m compuations
  ! 11DEC2012 LN translated from matlab file: DensityFunc.m
  real(dp),     intent(in) :: xT(:,:)
  real(dp),     intent(in) :: epsilon1(:),nu1(:)
  real(dp),     intent(in) :: Omega12(:,:)
  real(dp),     intent(in) :: C2(:,:)
  real(dp),     intent(in) :: CPsi(:,:)
  real(dp),     intent(in) :: Q(:,:)
  real(dp),     intent(in) :: S1(:)
  integer(i4b), intent(in) :: d1,d2
  real(dp),     intent(out) :: F(:)
  type(DensityGradType), intent(inout) :: GradF

  integer(i4b)          :: i1,i2
  integer(i4b)          :: nx
  real(dp), allocatable :: e2(:,:),z1T(:,:)
  real(dp)              :: Jacobian

  nx=size(xT,1)
  ! e2 = (nx x d2)
  allocate(e2(nx,d2))
  ! e2 = x*Q*inv(C2')
  ! z2 = x*Q
  ! C2 = (d2 x d2) lower triangular
  e2 = matmul(xT,Q)

  ! these loops compute e2_new = inv(C2') * e2'
  !                            = e2 * inv(C2)
  !                              size(e2) = nx x d2
  !
  ! exploit fact that C2 is lower triangular
  !    x * C2(:,1)            = e2(:,1)
  !    x(:,2:d2) * C2(2:d2,2) = e2(:,2)
  !    x(:,d2)   * C2(d2,d2)  = e2(:,d2)
  do i1=d2,-1,1
    if (i1==d2) then
      e2(:,i1) = e2(:,i1)/C2(d2,d2)
    else if (i1<d2) then
      e2(:,i1) = (matmul(e2(:,i1+1:d2),C2(i1+1:d2,i1)) -e2(:,i1))/C2(i1,i1)
    end if
  end do

!  do i1=1,nx
!  do i2=1,d2
!    if (i2==1) then
!      e2(i1,i2) = e2(i1,i2)/C2(i2,i2)
!    else
!      e2(i1,i2) =   &
!      (e2(i1,i2)-dot_product(e2(i1,1:i2-1),C2(1:i2-1,i2)))/C2(i2,i2)
!    end if
!  end do
!  end do

  ! size(v1) = (nx x d1)
  allocate(z1T(nx,d1))
  ! z1T = epsilon1 - nu1 - xT*Q*inv(C2)*Omega12'
  z1T  = spread(epsilon1 - nu1,1,nx) - matmul(e2,transpose(Omega12))

  ! compute z1T_new = z1T * inv(CPsi)
  ! exploit fact that CPsi is lower triangular
  do i1=1,nx
  do i2=1,d1
    if (i2==1) then
      z1T(i1,i2) = z1T(i1,i2)/CPsi(i2,i2)
    else
      z1T(i1,i2) =   &
      (z1T(i1,i2)-dot_product(z1T(i1,1:i2-1),CPsi(1:i2-1,i2)))/CPsi(i2,i2)
    end if
  end do
  end do

  Jacobian = product(S1)*(2.0d0*pi)**(-real(d1,dp)/2.0d0)
  ! det(S1)   = product(diag(S1)) since S1 is diagonal
  ! det(CPsi) = product of diagonal since CPsi is lower triangular
  do i1=1,d1
    Jacobian = Jacobian/CPsi(i1,i1)
  end do

  F   = abs(Jacobian)*dexp(-0.5d0*sum(z1T*z1T,2))

  ! gradient w.r.t. nu1
  !  df/dz1 = -f * z1T_new * inv(CPsi)^T
  ! size(z1T) = (nx x d1)
  ! size(f)   = (nx x 1)
  ! InvertLower  = invert lower triangular matrix
  GradF%nu1 = - spread(f,2,d1) * matmul(z1T,transpose(InvertLower(CPsi,d1)))
  !GradF%nu1 = - spread(f,2,d1) * matmul(z1T,inv(transpose(CPsi)))
  !print *,"Warning. Need to check on CPsi vs. transpose(CPsi)."

  ! gradient of F w.r.t. S1
  !    = F./S1    (since epsilon1 = inv(S1)*VT*p1 + S1*VT*q1 is fixed
  ! size(F) = (nx x 1)
  GradF%S1 = spread(F,2,d1) / spread(S1,1,nx)

  ! Gradient w.r.t. epsilon1 = - gradient w.r.t. nu1

  ! Gradient w.r.t. Omega12 = df/dz1 * dz1/dOmega12
  !                 dz1/dOmega12 = - inv(C2^T) * Q^T *x
  !                              =  -e2
  ! size(df/dz1)        = (nx x d1)
  ! size(e2)            = (nx x d1)
  ! size(GradF%Omega12) = (nx x d1 x d1)
  GradF%Omega12 = - spread(GradF%nu1,3,d1) * spread(e2,2,d1)

  !  Gradient w.r.t. (C2,Q)
  !  dF/d(C2) = df/dz1 * dz1/d(C2)
  !             dz1/d(C2) = - 

  deallocate(e2,z1T)
end subroutine DensityFunc2

subroutine SparseGridBayes(iFree,iuser,ruser)
  use nrtype
  use GlobalModule, only : SelectFreeType,MaxOptions
  use OutputModule, only : WriteBayes
  use nag_library,  only : D01ESF,D01ZKF
  implicit none
  type(SelectFreeType), intent(in) :: iFree
  integer(i4b),         intent(inout) :: iuser(*)
  real(dp),             intent(inout) :: ruser(*)

  ! D01ESF:   multidimensional Sparse Grid quadrature
  INTEGER(i4b)              :: NI,NDIM,LIOPTS,LOPTS,ifail
  integer(i4b), allocatable :: MAXDLV(:), IVALID(:), IOPTS(:)
  REAL (dp),    allocatable :: DINEST(:), ERREST(:), OPTS(:)
  CHARACTER(len=100)        :: OPTSTR

  ! set options for D01ESF
  ifail = -1
  LIOPTS = 100
  LOPTS  = 100
  NDIM   = iFree%nall
  NI     = 1 + ndim + ndim*(ndim+1)/2

  allocate(IOPTS(LIOPTS))
  allocate(OPTS(LOPTS))
  allocate(MAXDLV(NDIM))
  allocate(IVALID(NI))
  allocate(DINEST(NI))
  allocate(ERREST(NI))

  OPTSTR = "Initialize = D01ESF"
  MAXDLV = 0  ! default

  ! Initialize options
  call D01ZKF(OPTSTR, IOPTS, LIOPTS, OPTS, LOPTS, IFAIL)

  ! set options
  if (MaxOptions%AbsTol>0.0d0) then
    write(OPTSTR,'(a21,d12.4)') "Absolute Tolerance = ",MaxOptions%AbsTol
    ifail = -1
    call D01ZKF(OPTSTR, IOPTS, LIOPTS, OPTS, LOPTS, IFAIL)
  end if

  if (MaxOptions%RelTol>0.0d0) then
    write(OPTSTR,'(a21,d12.4)') "Relative Tolerance = ",MaxOptions%RelTol
    ifail = - 1
    call D01ZKF(OPTSTR, IOPTS, LIOPTS, OPTS, LOPTS, IFAIL)
  end if

  if (MaxOptions%MaxLevel>0) then
    write(OPTSTR,'(a21,i2)') "Maximum Level = ",MaxOptions%MaxLevel
    ifail=-1
    call D01ZKF(OPTSTR, IOPTS, LIOPTS, OPTS, LOPTS, IFAIL)
  end if 
  ifail = -1
  call D01ESF(NI, NDIM, IntegrateLikeFunc, MAXDLV, DINEST, ERREST, IVALID, IOPTS, OPTS, IUSER, RUSER, IFAIL)

  call WriteBayes(DINEST,ERREST,IVALID)
  deallocate(IOPTS)
  deallocate(OPTS)
  deallocate(MAXDLV)
  deallocate(IVALID)
  deallocate(DINEST)
  deallocate(ERREST)

end subroutine SparseGridBayes

SUBROUTINE IntegrateLikeFunc(NI, NDIM, NX, XTR, NNTR, ICOLZP, IROWIX, XS, QS, FM, IFLAG, IUSER, RUSER)
  implicit none
  INTEGER(i4b), intent(in)    ::  NI,NDIM,NX,NNTR,ICOLZP(NX+1),IROWIX(NNTR),QS(NNTR)
  integer(i4b), intent(inout) :: IFLAG, IUSER(*)
  REAL(dp),     intent(in)    :: XTR, XS(NNTR)
  real(dp),     intent(out)   :: FM(NI,NX)
  real(dp),     intent(inout) :: RUSER(*)

  real(dp), allocatable :: x0(:),x1(:)
  integer(i4b)          :: i1

  allocate(x0(ndim),x1(ndim))
  x0 = XTR
  do i1=1,nx
    x0(irowix(icolzp(i1):icolzp(i1+1)-1)) = xs(icolzp(i1):icolzp(i1+1)-1)
    ! change of variable from x0 to x1
    x1 = ChangeX(x0,ndim,iuser(1))
    call LikeFunc_QuadWrapper(x1,ndim,iuser,ruser,FM(:,i1))
    x0(irowix(icolzp(i1):icolzp(i1+1)-1)) = xtr
  end do

  deallocate(x0,x1)
end subroutine IntegrateLikeFunc

function ChangeX(x0,nx,model) result(x1)
  use nrtype
  use GlobalModule, only : ifree,bayes
  implicit none
  integer(i4b), intent(in) :: nx,model
  real(dp),     intent(in) :: x0(nx)
  real(dp)                 :: x1(nx)

  if (model==1) then
    print *, "Code for Bayesian estimation of model==1 not yet completed."
    print *, "Currently, Bayesian estimation only works for model==2."
    stop
  else if (model==2) then

    ! parameters impacting BD
    if (iFree%nBD_beta>0) then
      x1(iFree%xBD_beta) = bayes%BD_beta_lo + (bayes%BD_beta_hi - bayes%BD_beta_lo)*x0(iFree%xBD_beta)
    end if
    if (iFree%nBD_CDiag>0) then
      x1(iFree%xBD_CDiag) = bayes%BD_CDiag_lo + (bayes%BD_CDiag_hi - bayes%BD_CDiag_lo)*x0(iFree%xBD_CDiag)
    end if
    if (iFree%nBD_COffDiag>0) then
      x1(iFree%xBD_COffDiag) = bayes%BD_COffDiag_lo + (bayes%BD_COffDiag_hi - bayes%BD_COffDiag_lo)*x0(iFree%xBD_COffDiag)
    end if

    ! Parameters impacting BC
    if (iFree%nBC_beta>0) then
      x1(iFree%xBC_beta) = bayes%BC_beta_lo + (bayes%BC_beta_hi - bayes%BC_beta_lo)*x0(iFree%xBC_beta)
    end if
    if (iFree%nBC_CDiag>0) then
      x1(iFree%xBC_CDiag) = bayes%BC_CDiag_lo + (bayes%BC_CDiag_hi - bayes%BC_CDiag_lo)*x0(iFree%xBC_CDiag)
    end if
    if (iFree%nBC_COffDiag>0) then
      x1(iFree%xBC_COffDiag) = bayes%BC_COffDiag_lo + (bayes%BC_COffDiag_hi - bayes%BC_COffDiag_lo)*x0(iFree%xBC_COffDiag)
    end if

    ! MUE
    if (iFree%nMUE>0) then
      x1(iFree%xMUE) = bayes%MUE_lo + (bayes%MUE_hi-bayes%MUE_lo) * x0(iFree%xMUE)
    end if

    ! MUE_month
    if (iFree%nMUE_month>0) then
      x1(iFree%xMUE_month) = bayes%MUE_month_lo + (bayes%MUE_month_hi-bayes%MUE_month_lo) * x0(iFree%xmue_month)
    end if

    ! InvCDiag
    if (iFree%nInvCDiag>0) then
      x1(iFree%xInvCDiag) = bayes%InvCDiag_lo +(bayes%InvCDiag_hi - bayes%InvCDiag_lo)* x0(iFree%xInvCDiag)
    end if

    ! InvCOffDiag
    if (iFree%nInvCOffDiag>0) then
      x1(iFree%xInvCOffDiag) = bayes%InvCOffDiag_lo +(bayes%InvCOffDiag_hi - bayes%InvCOffDiag_lo)* x0(iFree%xInvCOffDiag)
    end if
  end if
end function ChangeX

subroutine LikeFunc_QuadWrapper(x,nx,iuser,ruser,F)
  use nrtype
  implicit none
  real(dp),     intent(in)    :: x(nx)
  integer(i4b), intent(in)    :: nx
  integer(i4b), intent(inout) :: iuser(*)
  real(dp),     intent(inout) :: ruser(*)
  real(dp),     intent(out)   :: F(:)

  integer(i4b) :: mode,nstate,i1,ix
  real(dp)     :: L,GradL(nx)

  mode=0
  nstate = 0

  call LikeFunc(mode,nx,x,L,GradL,nstate,iuser,ruser)
  F(1) = L
  F(2:nx+1) = x*L
  do i1=1,nx
    F(1+nx+i1*(i1-1)/2 + (/(ix,ix=1,i1)/)) = x(i1) * F(2:i1+1)
  end do

end subroutine LikeFunc_QuadWrapper

#ifdef BAYES
subroutine SetupBayesPrior(parms,iFree)
  use nrtype
  use GlobalModule, only : ParmsStructure,SelectFreeType,Bayes
  implicit none
  type (parmsStructure), intent(in) :: parms
  type(SelectFreeTYpe),  intent(in) :: iFree

  bayes%nAll = 10000
  allocate(bayes%x(iFree%nall,bayes%nall))
  allocate(bayes%w(bayes%nall))
  allocate(bayes%prior(bayes%nall))

  if (parms%model==1) then

  else if (parms%model==2) then

    if (iFree%flagMUE==1) then
      

  end if

  deallocate(bayes%x)
  deallocate(bayes%w)
  deallocate(bayes%prior)

end subroutine SetupBayesPrior
#endif

subroutine RunMonteCarlo(IMC1,IMC2,pid)
  use nrtype
  use GlobalModule, only : parms,parms0,iFree,            &
                           MaxOptions,                    &
                           ResultStructure,               &
                           ReadWriteParameters,CopyParameters, &
                           ControlOptions,MasterID
#if USE_MPI==0
  use DataModule, only   : CreateData,LoadData
#else
  use DataModule, only   : CreateData,SendData,LoadData
  use GlobalModule, only : BroadcastIFree,BroadcastParms
#endif
  use OutputModule, only : SaveMCOutputs
  implicit none
  integer(i4b), intent(in) :: IMC1,IMC2
  integer(i4b), intent(in) :: pid
  integer(i4b)             :: IMC
  real(dp), allocatable    :: x(:),MCX(:,:),MCLambda(:)
  real(dp), allocatable    :: Grad(:),Hess(:,:)
  real(dp)                 :: LVALUE
  type(ResultStructure)    :: stats
  integer(i4b)             :: NMC
  integer(i4b)             :: ifail
  integer(i4b)             :: DateTime(8)
  logical                  :: ExistFlag

  ! Copy parameters from parms to parms0.
  ! Then read parameters from disk
  inquire(file=parms%file,exist=ExistFlag)
  if (pid==MasterID .and. ControlOptions%HotStart==1 .and. ExistFlag) then
    call CopyParameters(parms,parms0)
    call ReadWriteParameters(parms,'read')
  end if

  ! Set initial value of iFree and broadcast iFree
  if (pid==MasterID) then
    call SelectFreeParameters(parms,iFree)
    !if (MaxOptions%Algorithm==6) then
    !  call SetupBayesPrior(parms,iFree)
    !end if
  end if

#if USE_MPI==1
  call BroadcastParms(parms,pid)
  print *,'pid = ',pid,'Broadcastparms complete.'
  call BroadcastIFree(pid)
  print *,'pid = ',pid,'BroadcastIFree complete.'
#endif

  ! copy true parameters from parms to parms0
  call CopyParameters(parms,parms0)

  NMC = IMC2-IMC1+1
  do iMC=IMC1,IMC2
    ! reset parms to be equal to parms0
    if (iMC>IMC1) then
      call CopyParameters(parms0,parms)
    end if

    if (pid==MasterID) then
      if (ControlOptions%SimulateData==1) then
        call CreateData(iMC)
      else
        call date_and_time(values=DateTime)
        print *,"Begin load data. (day,hour,min) = ",DateTime(3),DateTime(5),DateTime(6)
        call LoadData
        call date_and_time(values=DateTime)
        print *,"Completed load data. (day,hour,min) = ",DateTime(3),DateTime(5),DateTime(6)
      end if
    end if

#if USE_MPI==1
    call SendData(pid)
#endif

    if (iFree%OneAtATime==1) then
       call MaxOneAtTime(pid)
    else
      if (pid==MasterID) then
        allocate(x(iFree%NALL))
        allocate(MCX(iFree%NALL,NMC+1))
        allocate(MCLambda(NMC))
        allocate(Grad(iFree%NALL),Hess(iFree%NALL,iFree%NALL))
        call ComputeInitialGuess(parms,iFree,x)
        MCX(:,NMC+1) = x
        call MaximizeLikelihood(x,LValue,Grad,Hess,ifail)
        MCX(:,IMC-IMC1+1) = x
        ! update parms
        call UpdateParms2(x,iFree,parms)
        ! b) save parms to disk
        call ReadWriteParameters(parms,'write')
        deallocate(x,Grad,Hess)
        deallocate(MCX,MCLambda)
      elseif (pid .ne. MasterID) then
#if USE_MPI>0
        call WorkerTask(parms%model,pid)
#endif
      end if ! if (pid==MasterID) then
    end if   ! if (iFree%OneAtATime==1) then
  end do     ! do IMC=IMC1,IMC2
end subroutine RunMonteCarlo

! loop over elements of x
! in each iteration, max w.r.t. current element of x
! one element of x at a time
subroutine MaxOneAtTime(pid)
  use nrtype
#if USE_MPI=1
  use mpi
#endif
  use GlobalModule, only : MasterID,iFree,SelectFreeType,parms, &
                           ReadWriteParameters,DeallocateIFree

  implicit none
  integer(i4b), intent(in) :: pid
  type(SelectFreeType)     :: iFree0
  integer(i4b)             :: i1,ifail,ierr
  real(dp), allocatable    :: x(:),Grad(:),Hess(:,:)
  real(dp)                 :: LValue

  call CopyIFree(iFree,iFree0)

  do i1=1,iFree0%NALL
    call DeallocateIFree(iFree)
    call UpdateIFree(i1,iFree0,iFree)
    allocate(x(iFree%NALL))
    allocate(Grad(iFree%NALL),Hess(iFree%NALL,iFree%NALL))

    if (pid==MasterID) then
      call ComputeInitialGuess(parms,iFree,x)
      call MaximizeLikelihood(x,LValue,Grad,Hess,ifail)
    elseif (pid .ne. MasterID) then
#if USE_MPI>0
      call WorkerTask(parms%model,pid)
#endif
    end if

#if USE_MPI>0
    call mpi_bcast(x,size(x,1),MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr)
#endif

    ! update parms
    call UpdateParms2(x,iFree,parms)
    ! b) save parms to disk
    call ReadWriteParameters(parms,'write')
    deallocate(x,Grad,Hess)
  end do

  call DeallocateIFree(iFree)
  call CopyIFree(iFree0,iFree)
  call DeallocateIFree(iFree0)

end subroutine MaxOneAtTime

! Copy iFree to iFreeCopy
subroutine CopyIFree(iFree,iFreeCopy)
  use GlobalModule, only : SelectFreeType
  implicit none
  type(SelectFreeType), intent(in) :: iFree
  type(SelectFreeType), intent(inout) :: iFreeCopy

  iFreeCopy%ND           = iFree%ND
  iFreeCopy%NBC          = iFree%NBC
  iFreeCopy%NMUE         = iFree%NMUE
  iFreeCopy%NMUE_month   = iFree%NMUE_month
  iFreeCopy%NINVCDIAG    = iFree%NINVCDIAG
  iFreeCopy%NINVCOFFDIAG = iFree%NINVCOffDiag

  iFreeCopy%NBD_beta     = iFree%NBD_beta
  iFreeCopy%NBD_CDIAG    = iFree%NBD_CDIAG
  iFreeCopy%NBD_COFFDIAG = iFree%NBD_COFFDIAG
  iFreeCopy%nBD_month    = iFree%nBD_month

  iFreeCopy%NBC_beta     = iFree%NBC_beta
  iFreeCopy%NBC_CDIAG    = iFree%NBC_CDIAG
  iFreeCopy%NBC_COFFDIAG = iFree%NBC_COFFDIAG

  iFreeCopy%NALL         = iFREE%NALL

  iFreeCopy%flagD           = iFree%flagD
  iFreeCopy%flagBC          = iFree%flagBC
  iFreeCopy%flagMUE         = iFree%flagMUE
  iFreeCopy%flagMUE_month   = iFree%flagMUE_month
  iFreeCopy%flagInvCDiag    = iFree%flagInvCDiag
  iFreeCopy%flagInvCOffDiag = iFree%flagInvCOffDiag

  iFreeCopy%flagBC_beta     = iFree%flagBC_beta
  iFreeCopy%flagBC_CDiag    = iFree%flagBC_CDiag
  iFreeCopy%flagBC_COffDiag = iFree%flagBC_COffDiag
  iFreeCopy%flagBD_month    = iFree%flagBD_month

  iFreeCopy%flagBD_beta     = iFree%flagBD_beta
  iFreeCopy%flagBD_CDiag    = iFree%flagBD_CDiag
  iFreeCopy%flagBD_COffDiag = iFree%flagBD_COffDiag
  iFreeCopy%OneAtATime      = iFree%OneAtATime

  if (iFree%flagD>0) then
    allocate(iFreeCopy%D(iFree%ND))
    allocate(iFreeCopy%xD(iFree%ND))
    iFreeCopy%D = iFree%D
    iFreeCopy%xD = iFree%xD
  end if

  if (iFree%flagBC>0) then
    allocate(iFreeCopy%BC(iFree%NBC))
    allocate(iFreeCopy%xBC(iFree%NBC))
    iFreeCopy%BC = iFree%BC
    iFreeCopy%xBC = iFree%xBC
  end if

  if (iFree%flagMUE>0) then
    allocate(iFreeCopy%MUE(iFree%NMUE))
    allocate(iFreeCopy%xMUE(iFree%NMUE))
    iFreeCopy%MUE = iFree%MUE
    iFreeCopy%xMUE = iFree%xMUE
  end if

  if (iFree%flagMUE_month>0) then
    allocate(iFreeCopy%MUE_month(iFree%NMUE_month))
    allocate(iFreeCopy%xMUE_month(iFree%NMUE_month))
    iFreeCopy%MUE_month = iFree%MUE_month
    iFreeCopy%xMUE_month = iFree%xMUE_month
  end if

  if (iFree%flagInvCDiag>0) then
    allocate(iFreeCopy%InvCDiag(iFree%NInvCDiag))
    allocate(iFreeCopy%xInvCDiag(iFree%NInvCDiag))
    iFreeCopy%InvCDiag = iFree%InvCDiag
    iFreeCopy%xInvCDiag = iFree%xInvCDiag
  end if

  if (iFree%flagInvCOffDiag>0) then
    allocate(iFreeCopy%InvCOffDiag(iFree%NInvCOffDiag))
    allocate(iFreeCopy%xInvCOffDiag(iFree%NInvCOffDiag))
    iFreeCopy%InvCOffDiag = iFree%InvCOffDiag
    iFreeCopy%xInvCOffDiag = iFree%xInvCOffDiag
  end if

  if (iFree%flagBC_beta>0) then
    allocate(iFreeCopy%BC_beta(iFree%NBC_beta))
    allocate(iFreeCopy%xBC_beta(iFree%NBC_beta))
    iFreeCopy%BC_beta = iFree%BC_beta
    iFreeCopy%xBC_beta = iFree%xBC_beta
  end if

  if (iFree%flagBC_CDiag>0) then
    allocate(iFreeCopy%BC_CDiag(iFree%NBC_CDiag))
    allocate(iFreeCopy%xBC_CDiag(iFree%NBC_CDiag))
    iFreeCopy%BC_CDiag = iFree%BC_CDiag
    iFreeCopy%xBC_CDiag = iFree%xBC_CDiag
  end if

  if (iFree%flagBC_COffDiag>0) then
    allocate(iFreeCopy%BC_COffDiag(iFree%NBC_COffDiag))
    allocate(iFreeCopy%xBC_COffDiag(iFree%NBC_COffDiag))
    iFreeCopy%BC_COffDiag = iFree%BC_COffDiag
    iFreeCopy%xBC_COffDiag = iFree%xBC_COffDiag
  end if

  if (iFree%flagBD_beta>0) then
    allocate(iFreeCopy%BD_beta(iFree%NBD_beta))
    allocate(iFreeCopy%xBD_beta(iFree%NBD_beta))
    iFreeCopy%BD_beta = iFree%BD_beta
    iFreeCopy%xBD_beta = iFree%xBD_beta
  end if
  
  if (iFree%flagBD_month>0) then
    allocate(iFreeCopy%BD_month(iFree%NBD_month))
    allocate(iFreeCopy%xBD_month(iFree%NBD_month))
    iFreeCopy%BD_month = iFree%BD_month
    iFreeCopy%xBD_month = iFree%xBD_month
  end if

  if (iFree%flagBD_CDiag>0) then
    allocate(iFreeCopy%BD_CDiag(iFree%NBD_CDiag))
    allocate(iFreeCopy%xBD_CDiag(iFree%NBD_CDiag))
    iFreeCopy%BD_CDiag = iFree%BD_CDiag
    iFreeCopy%xBD_CDiag = iFree%xBD_CDiag
  end if

  if (iFree%flagBD_COffDiag>0) then
    allocate(iFreeCopy%BD_COffDiag(iFree%NBD_COffDiag))
    allocate(iFreeCopy%xBD_COffDiag(iFree%NBD_COffDiag))
    iFreeCopy%BD_COffDiag = iFree%BD_COffDiag
    iFreeCopy%xBD_COffDiag = iFree%xBD_COffDiag
  end if

end subroutine CopyIFree

! Set iFreeNew equal to element i1 of iFree0
subroutine UpdateIFree(i1,iFree0,iFreeNew)
  use GlobalModule, only : SelectFreeType
  implicit none
  integer(i4b), intent(in) :: i1
  type(SelectFreeType), intent(in)    :: iFree0
  type(SelectFreeType), intent(inout) :: iFreeNew

  iFreeNew%nall         = 1
  iFreeNew%nD           = 0
  iFreeNew%nBC          = 0
  iFreeNew%nMUE         = 0
  iFreeNew%nMUE_month   = 0
  iFreeNew%nInvCDiag    = 0
  iFreeNew%nInvCOffDiag = 0
  iFreeNew%nBC_beta     = 0
  iFreeNew%nBC_CDiag    = 0
  iFreeNew%nBC_COffDiag = 0
  iFreeNew%nBD_beta     = 0
  iFreeNew%nBD_CDiag    = 0
  iFreeNew%nBD_COffDiag = 0
  iFreeNew%nBD_month    = 0
  
  iFreeNew%flagD           = 0
  iFreeNew%flagBC          = 0
  iFreeNew%flagMUE         = 0
  iFreeNew%flagMUE_month   = 0
  iFreeNew%flagInvCDiag    = 0
  iFreeNew%flagInvCOffDiag = 0
  iFreeNew%flagBC_beta     = 0
  iFreeNew%flagBC_CDiag    = 0
  iFreeNew%flagBC_COffDiag = 0
  iFreeNew%flagBD_beta     = 0
  iFreeNew%flagBD_CDiag    = 0
  iFreeNew%flagBD_COffDiag = 0
  iFreeNew%flagBD_month    = 0
 
  allocate(iFreeNew%xlabels(1))
  write(iFreeNew%xlabels(1),'(i3.0)') i1

  if (iFree0%nD>0) then
    if (any(iFree0%xD==i1)) then
      allocate(iFreeNew%D(1))
      allocate(IFreeNew%xD(1))
      iFreeNew%xD = 1
      iFreeNew%D  = pack(iFree0%D,iFree0%xD==i1)
      iFreeNew%nD = 1
      iFreeNew%flagD = 1
    end if
  end if

  if (iFree0%nBC>0) then
    if (any(iFree0%xBC==i1)) then
      allocate(iFreeNew%BC(1))
      allocate(IFreeNew%xBC(1))
      iFreeNew%xBC = 1
      iFreeNew%BC  = pack(iFree0%BC,iFree0%xBC==i1)
      iFreeNew%nBC = 1
      iFreeNew%flagBC = 1
    end if
  end if

  if (iFree0%nMUE>0) then
    if (any(iFree0%xMUE==i1)) then
      allocate(iFreeNew%MUE(1))
      allocate(IFreeNew%xMUE(1))
      iFreeNew%xMUE = 1
      iFreeNew%MUE  = pack(iFree0%MUE,iFree0%xMUE==i1)
      iFreeNew%nMUE = 1
      iFreeNew%flagMUE = 1
    end if
  end if

  if (iFree0%nMUE_month>0) then
    if (any(iFree0%xMUE_month==i1)) then
      allocate(iFreeNew%MUE_month(1))
      allocate(IFreeNew%xMUE_month(1))
      iFreeNew%xMUE_month = 1
      iFreeNew%MUE_month  = pack(iFree0%MUE_month,iFree0%xMUE_month==i1)
      iFreeNew%nMUE_month = 1
      iFreeNew%flagMUE_month = 1
    end if
  end if

  if (iFree0%nInvCDiag>0) then
    if (any(iFree0%xInvCDiag==i1)) then
      allocate(iFreeNew%InvCDiag(1))
      allocate(IFreeNew%xInvCDiag(1))
      iFreeNew%xInvCDiag = 1
      iFreeNew%InvCDiag  = pack(iFree0%InvCDiag,iFree0%xInvCDiag==i1)
      iFreeNew%nInvCDiag = 1
      iFreeNew%flagInvCDiag = 1
    end if
  end if

  if (iFree0%nInvCOffDiag>0) then
    if (any(iFree0%xInvCOffDiag==i1)) then
      allocate(iFreeNew%InvCOffDiag(1))
      allocate(IFreeNew%xInvCOffDiag(1))
      iFreeNew%xInvCOffDiag = 1
      iFreeNew%InvCOffDiag  = pack(iFree0%InvCOffDiag,iFree0%xInvCOffDiag==i1)
      iFreeNew%nInvCOffDiag = 1
      iFreeNew%flagInvCOffDiag = 1
    end if
  end if

  if (iFree0%nBC_beta>0) then
    if (any(iFree0%xBC_beta==i1)) then
      allocate(iFreeNew%BC_beta(1))
      allocate(IFreeNew%xBC_beta(1))
      iFreeNew%xBC_beta = 1
      iFreeNew%BC_beta  = pack(iFree0%BC_beta,iFree0%xBC_beta==i1)
      iFreeNew%nBC_beta = 1
      iFreeNew%flagBC_beta = 1
    end if
  end if

  if (iFree0%nBD_beta>0) then
    if (any(iFree0%xBD_beta==i1)) then
      allocate(iFreeNew%BD_beta(1))
      allocate(IFreeNew%xBD_beta(1))
      iFreeNew%xBD_beta = 1
      iFreeNew%BD_beta  = pack(iFree0%BD_beta,iFree0%xBD_beta==i1)
      iFreeNew%nBD_beta = 1
      iFreeNew%flagBD_beta = 1
    end if
  end if

  if (iFree0%nBD_month>0) then
    if (any(iFree0%xBD_month==i1)) then
      allocate(iFreeNew%BD_month(1))
      allocate(IFreeNew%xBD_month(1))
      iFreeNew%xBD_month = 1
      iFreeNew%BD_month  = pack(iFree0%BD_month,iFree0%xBD_month==i1)
      iFreeNew%nBD_month = 1
      iFreeNew%flagBD_month = 1
    end if
  end if

  if (iFree0%nBC_CDiag>0) then
    if (any(iFree0%xBC_CDiag==i1)) then
      allocate(iFreeNew%BC_CDiag(1))
      allocate(IFreeNew%xBC_CDiag(1))
      iFreeNew%xBC_CDiag = 1
      iFreeNew%BC_CDiag  = pack(iFree0%BC_CDiag,iFree0%xBC_CDiag==i1)
      iFreeNew%nBC_CDiag = 1
      iFreeNew%flagBC_CDiag = 1
    end if
  end if

  if (iFree0%nBD_CDiag>0) then
    if (any(iFree0%xBD_CDiag==i1)) then
      allocate(iFreeNew%BD_CDiag(1))
      allocate(IFreeNew%xBD_CDiag(1))
      iFreeNew%xBD_CDiag = 1
      iFreeNew%BD_CDiag  = pack(iFree0%BD_CDiag,iFree0%xBD_CDiag==i1)
      iFreeNew%nBD_CDiag = 1
      iFreeNew%flagBD_CDiag = 1
    end if
  end if

  if (iFree0%nBC_COffDiag>0) then
    if (any(iFree0%xBC_COffDiag==i1)) then
      allocate(iFreeNew%BC_COffDiag(1))
      allocate(IFreeNew%xBC_COffDiag(1))
      iFreeNew%xBC_COffDiag = 1
      iFreeNew%BC_COffDiag  = pack(iFree0%BC_COffDiag,iFree0%xBC_COffDiag==i1)
      iFreeNew%nBC_COffDiag = 1
      iFreeNew%flagBC_COffDiag = 1
    end if
  end if

  if (iFree0%nBD_COffDiag>0) then
    if (any(iFree0%xBD_COffDiag==i1)) then
      allocate(iFreeNew%BD_COffDiag(1))
      allocate(IFreeNew%xBD_COffDiag(1))
      iFreeNew%xBD_COffDiag = 1
      iFreeNew%BD_COffDiag  = pack(iFree0%BD_COffDiag,iFree0%xBD_COffDiag==i1)
      iFreeNew%nBD_COffDiag = 1
      iFreeNew%flagBD_COffDiag = 1
    end if
  end if

end subroutine UpdateIFree


!!   BEGIN ANALYSIS of DEMAND
subroutine AnalyseResults(IMC)
  use nrtype
  use GlobalModule, only : ControlOptions,parms,ReadWriteParameters
  use DataModule, only   : CreateData,LoadData
  implicit none
  integer(i4b), intent(in) :: IMC
  integer(i4b)             :: DateTime(8)

  if (ControlOptions%HotStart==1) then
    call ReadWriteParameters(parms,'read')
  end if

  if (ControlOptions%SimulateData==1) then
    call CreateData(iMC)
  else
    call date_and_time(values=DateTime)
    print *,"Begin load data. (day,hour,min) = ",DateTime(3),DateTime(5),DateTime(6)
    call LoadData
    call date_and_time(values=DateTime)
    print *,"Completed load data. (day,hour,min) = ",DateTime(3),DateTime(5),DateTime(6)
  end if

  ! reset prices = average price
  call ComputeElasticities

end subroutine AnalyseResults

subroutine ComputeElasticities
#ifndef __GFORTRAN__
  use IFPORT
#endif
  use nrtype
  use GlobalModule, only : parms,HHData,ControlOptions,InputDir, &
                           DataStructure,AllocateLocalHHData,    &
                           DeallocateLocalHHData,LoadBasePrice,  &
                           LoadTaxParameters,OutDir,ParmFiles
  use nag_library, only  : G05KFF,G05RZF,G05SAF
  use OutputModule, only : WriteElasticities,WriteDemandResults1,WriteDemandResults2, &
                           WriteTaxResults1,WriteTaxResults2, &
                           WritePrediction
  use DataModule,  only  : DrawRandomCoefficients
  use ToolsModule, only  : InverseNormal_mkl
  implicit none

  integer(i4b)              :: i1,i2,i3
  type(DataStructure)       :: HHDataSim1,HHDataSim2,SimData1,HHFit

  ! variables to compute demand
  real(dp)                  :: h
  real(dp), allocatable     :: q0(:),q1(:),qdata(:),qhat(:),ExcessQ(:)
  real(dp), allocatable     :: GradQ(:,:)
  real(dp), allocatable     :: elas(:,:)
  real(dp), allocatable     :: newQ(:,:),p(:),p0(:)
  integer(i4b)              :: np,nProb

  ! variables to store results from tax experiments
  integer(i4b)                    :: ntax
  integer(i4b),       allocatable :: taxid(:)
  character(len=100), allocatable :: taxlabel(:)
  character(len=10),  allocatable :: taxtype(:)
  real(dp),           allocatable :: tax(:,:)
  real(dp),           allocatable :: qtax(:,:),ptax(:,:)
  real(dp),           allocatable :: HHSpending(:,:),HHUtility(:,:)
  character(len=1000)             :: copyfilecommand

  ! initialize random number generator
  integer(i4b)              :: genid,subid,lseed,ifail,lstate
  integer(i4b), allocatable :: localseed(:),state(:)

  ! variables used by E04RZF: normal random number generator
!  integer(i4b)              :: mode,ldc,ldx,LR
!  real(dp), allocatable     :: R(:)
!  real(dp), allocatable     :: e(:,:)
!  real(dp), allocatable     :: eye(:,:),zero(:)

  ! SimData1 will simulate results for a small number of individual people
  !   nprob = number of values of eta in each dimension
  !   N     = total number of individuals in simulation
  !   prob  = quantiles of eta distribution at which results will be computed
  !   eta   = corresponding points in eta distribution
  real(dp), allocatable     :: eta(:),prob(:)
  real(dp), allocatable     :: e2(:)

  nprob      = 3
  SimData1%N = 9  ! simulate and graph for 9 people with varying levels of 
                  ! (eta1,eta2)
  call AllocateLocalHHData(SimData1)
  ! Define month
  SimData1%month = 5  ! all are May
  allocate(eta(nprob),prob(nprob))
  prob = (/0.25d0,0.5d0,0.75d0/)
  ! SimData1% e = MuE for all people
  eta = InverseNormal_mkl(prob,nprob,ifail)
  SimData1%e        = spread(Parms%MuE,2,SimData1%N)
  do i1=1,12
    SimData1%e = SimData1%e &
       + merge(spread(parms%mue_month(:,i1),2,SimData1%N),0.0d0,spread(SimData1%month,1,12)==i1)
  end do
  SimData1%eta(1,:) = reshape(spread(eta,2,nprob),(/nprob*nprob/))
  SimData1%eta(2,:) = reshape(spread(eta,1,nprob),(/nprob*nprob/))

  ! allocate memory for HHData0
  HHDataSim1%N = HHData%NSim   ! NSIM might be different from N
  call AllocateLocalHHData(HHDataSim1)

  ! initialize random number generator
  !  external G05KFF  ! NAG routine to set seed for random number generator
  !  external G05SAF  ! NAG routine to generate uniform random numbers

  genid  = 3 ! Mersenne twister algorithm
  subid  = 0 ! not used when genid = 3
  !lseed = 624 ! number of seeds needed for genid = 3
  !allocate(seed(lseed))
  lseed  = 1
  allocate(localseed(lseed))
  lstate = 1260  ! min value for genid=3
  allocate(state(lstate))
  ifail  = -1
  localseed   = HHData%SimulationSeed

  call G05KFF(genid,subid,localseed,lseed,state,lstate,ifail)

  ! set price and month
  if (HHDataSim1%N==HHData%N) then
    HHDataSim1%p = HHData%p
    HHDataSim1%month = HHData%month
  else
    print *,'HHDataSim1%N = ',HHDataSim1%N,'.'
    print *,'Sample prices with replacement from raw data.'
    allocate(e2(HHDataSim1%N))
    call G05SAF(HHDataSim1%N,state,e2,ifail)
    do i1=1,HHDataSim1%N
      i2 =ceiling(HHData%n * e2(i1))
      HHDataSim1%p(:,i1) = HHData%p(:,i2)
      HHDataSim1%month(i1) = HHData%month(i2)
    end do
    deallocate(e2)
  end if

  ! Generate (e,eta)
  call DrawRandomCoefficients(HHDataSim1,state)

  ! Copy HHData to HHDataFit
  HHFit%N     = HHData%N
  call AllocateLocalHHData(HHFit)
  HHFit%p     = HHData%p
  HHFit%month = HHData%month
  call DrawRandomCoefficients(HHFit,state)
  ! Compute predictions for fitted data
  call ComputeDemand(HHFit)

  call WritePrediction(HHData,HHFit)

  ! Aggregate demand from data
  ! baseline demand
  allocate(q0(Parms%J),q1(parms%J),qdata(parms%J),qhat(parms%j))
  allocate(ExcessQ(parms%J))
  allocate(p0(parms%J))
  q0 = 0.0d0
  q1 = 0.0d0
  qdata=0.0d0
  qhat =0.0d0
 
  ! aggregate demand in raw data
  call ComputeAggregateDemand(HHData%q,HHData%iNonZero,qdata,HHData%nNonZero,0)

  ! predicted demand and welfare at prices in data
  call ComputeDemand(HHDataSim1)
  call ComputeAggregateDemand(HHDataSim1%q,HHDataSim1%iNonZero,qhat,HHDataSim1%nNonZero,0)

  ! set price = average price
  call LoadBasePrice(p0)
  p0 = sum(HHDataSim1%p,2)/real(HHDataSim1%n,dp) + p0
  HHDataSim1%p = spread(p0,2,HHDataSim1%N)
  SimData1%p = spread(sum(HHDataSim1%p,2)/real(HHDataSim1%n,dp),2,SimData1%N)

  ! Compute baseline demand and welfare
  call ComputeDemand(HHDataSim1)
  call ComputeDemand(SimData1)
  call ComputeAggregateDemand(HHDataSim1%q,HHDataSim1%iNonZero,q0,HHDataSim1%nNonZero,0)

  ExcessQ = q0-qdata
  call WriteDemandResults1(qdata,qhat,q0,p0)

  ! copy exogenous variables from HHDataSim1 to HHDataSim2
  HHDataSim2%N = HHDataSim1%N
  call AllocateLocalHHData(HHDataSim2)
  HHDataSim2%e   = HHDataSim1%e
  HHDataSim2%eta = HHDataSim1%eta
  HHDataSim2%p   = HHDataSim1%p
  HHDataSim2%month = HHDataSim1%month

  ! Compute slope of demand
  allocate(GradQ(parms%J,parms%J))
  allocate(elas(parms%J,parms%J))
  GradQ = 0.0d0
  elas  = 0.0d0
  h = 1.0d-4
  do i1=1,parms%J
    ! size(p) = (J x N)
    HHDataSim2%p(i1,:) = HHDataSim1%p(i1,:) + h
    call ComputeDemand(HHDataSim2)
    call ComputeAggregateDemand(HHDataSim2%q,HHDataSim2%iNonZero,q1,HHDataSim2%nNonZero,0)
    GradQ(:,i1) = (q1-q0)/(HHDataSim2%p(i1,1)-HHDataSim1%p(i1,1))
    elas(:,i1)  = HHDataSim1%p(i1,1)*GradQ(:,i1)/q0
    HHDataSim2%p(i1,:) = HHDataSim1%p(i1,:)
    print *,'Computing elasticities: ',i1,' of ',parms%J
  end do

  call WriteElasticities(elas)

  ! Plot demand functions w.r.t. own price
  !   for each good save
  !   p(i1), q(i1), q1,q2,q3,q4,q5,n1,n2,n3,n4,n5
  np = parms%nPrices_plot
  allocate(newq(np,parms%k+1),p(np))
  newq = 0.0d0
  do i1=1,parms%J
    p    = HHDataSim1%p(i1,1) * (0.5d0+(2.0d0-0.5d0)*real((/0:np-1/),dp)/real(np-1,dp))
    do i2=1,np
      HHDataSim2%p(i1,:) = p(i2)
      call ComputeDemand(HHDataSim2)
      do i3=0,parms%K
        ! Compute q(i1) conditional on nNonZero==n
        !    if n==0, then compute total demand
        !    save results in matrix newq
       call ComputeAggregateDemand(HHDataSim2%q,HHDataSim2%iNonZero,q1,HHDataSim2%nNonZero,i3)
       newq(i2,i3+1) = q1(i1)
      end do
    end do
    call WriteDemandResults2(p,NewQ,i1)
  end do

  ! Compute demands for 9 individuals
  ! write output to file
  call ComputeIndividualDemand(SimData1)

  ! Compute tax and other counterfactuals
  !  1) raise price of 1 good to a very large number
  !  2) raise price of foreign goods due to Brexit
  !  3) subsidy to fruit consumption
  ! set price = average price

  call LoadTaxParameters(taxid,taxlabel,taxtype,tax)
  ntax = size(taxid,1)
  allocate(qtax(parms%J,ntax),ptax(parms%J,ntax))
  allocate(HHSpending(HHDataSim2%n,ntax),HHUtility(HHDataSim2%n,ntax))

  qtax       = 0.0d0
  ptax       = 0.0d0  
  HHSpending = 0.0d0
  HHUtility  = 0.0d0
  do i1=1,ntax
    if (taxtype(i1)=="ad valorem") then
      ptax(:,i1)   = (sum(HHDataSim1%p,2)/real(HHDataSim1%n)) * tax(:,i1)
    else if (taxtype(i1)=="excise") then
      ptax(:,i1)   = (sum(HHDataSim1%p,2)/real(HHDataSim1%n)) + tax(:,i1)
    end if
    HHDataSim2%p = spread(ptax(:,i1),2,HHDataSim2%N)
    call ComputeDemand(HHDataSim2)
    HHSpending(:,i1) = HHDataSim2%expenditure
    HHUtility(:,i1)  = HHDataSim2%utility
    call ComputeAggregateDemand(HHDataSim2%q,HHDataSim2%iNonZero,qtax(:,i1),HHDataSim2%nNonZero,0)
  end do

  ! Write description of tax experiments
  copyfilecommand = "cp " // trim(ParmFiles%TaxparmsFile) // " " // trim(OutDir) // "/taxparms.csv"
  ifail = system(trim(copyfilecommand))
  
  ! Write aggregate results:  baseline (q0,p0) and alternative (qtax,ptax) (J x ...)
  call WriteTaxResults1(q0,p0,qtax,ptax,filename="taxresults_aggregate.csv")

  ! Write HH results on (expenditure,utility)
  !    HHDataSim1%expenditure,HHDataSim1%utility = baseline results
  !    (HHSpending,HHUtility) = alternative results
  call WriteTaxResults2(HHDataSim1%expenditure,HHDataSim1%utility, &
                        HHSpending,HHUtility,filename="taxresults_hh.csv")

  ! deallocate memory for HHData0
  call DeallocateLocalHHData(HHDataSim1)
  call DeallocateLocalHHData(HHDataSim2)
  call DeallocateLocalHHData(HHFit)
  deallocate(q0,q1,p0,GradQ,ExcessQ)
  deallocate(qdata,qhat)
  deallocate(newq)
  deallocate(localseed,state)
  deallocate(eta,prob)
  deallocate(tax,qtax,ptax)
  deallocate(taxlabel)
  deallocate(HHSpending,HHUtility)

end subroutine ComputeElasticities

subroutine ComputeIndividualDemand(SimData1)
  use nrtype
  use GlobalModule, only : DataStructure,parms
  use OutputModule, only : MakeFullFileName
  implicit none
  type(DataStructure), intent(inout) :: SimData1

  integer(i4b)              :: nHH,np,SimUnit
  integer(i4b)              :: i1,i2,ihh
  real(dp), allocatable     :: p(:)
  character(len=99)         :: SimFile,fmt1

  nHH = SimData1%n
  np  = 30
  allocate(p(np))
  p = 0.0d0

  SimUnit = 45
  SimFile = MakeFullFileName('IndividualDemand.txt')
  open(Unit= SimUNit, &
       File= SimFile, &
       action = 'write')
  
  write(fmt1,'(a11,i1,a18,i2,a9,i2,a12,i3,a14)') &
             '(3(i2,","),',parms%dim_eta,'(g12.4,","),i2,","', &
             parms%k,'(i3,","),',parms%k,'(g12.4,","),',       &
             parms%J,'(g12.4,:,","))'

  do i1=1,parms%J
    p = SimData1%p(i1,1) * (0.5d0+(2.0d0-0.5d0)*real((/0:np-1/),dp)/real(np-1,dp))
    do i2=1,np
      SimData1%p(i1,:) = p(i2)
      call ComputeDemand(SimData1)

      do ihh=1,nHH
        ! iHH  i1 i2  eta(1)  eta(2)  nNonZero i1-i5 q1-q5 p1-pJ
        write(SimUnit,fmt1) ihh,i1,i2,SimData1%eta(:,ihh),SimData1%nNonZero(ihh), &
                            SimData1%iNonZero(:,ihh),SimData1%q(:,ihh), &
                            SimData1%p(:,ihh)
      end do
    end do
  end do
  close(SimUnit)

  deallocate(p)
end subroutine ComputeIndividualDemand

subroutine ComputeAggregateDemand(q,iNonZero,TotalQ,nNonZero,n)
use nrtype
use GlobalModule, only : parms
implicit none
real(dp),     intent(in)  :: q(:,:)
integer(i4b), intent(in)  :: iNonZero(:,:)
real(dp),     intent(out) :: TotalQ(:)
integer(i4b), intent(in)  :: nNonZero(:),n
! size(q)        = (K x N)
! size(inonzero) = (K x N)
! size(TotalQ)   = (J x 1)
! size(nnonzero) = (N x 1)
! n              = (1 x 1)  used to compute demand conditional on n
!                           for n=1,...,K
!                           n=0 computes total demand
integer(i4b) :: i1,i2

TotalQ = 0.0d0
do i1=1,parms%J
do i2=1,parms%K
  if (n==0) then
    TotalQ(i1) = TotalQ(i1)+sum(q(i2,:),mask=iNonZero(i2,:)==i1)
  elseif (n>0 .and. n<=parms%K) then
    TotalQ(i1) = TotalQ(i1)+sum(q(i2,:),mask=(iNonZero(i2,:)==i1 .and. nNonZero==n))
  endif
end do
end do

! Per capita demand
TotalQ = TotalQ / real(size(q,2),dp)
end subroutine ComputeAggregateDemand

subroutine ComputeDemand(HHData1)
  use GlobalModule, only : parms,DataStructure,InputDir
  use DataModule,   only : ComputeCurrentB
  use nag_library,  only : X04ABF,X04ACF,X04ADF, &
                           E04WBF,E04NCA,E04NDA
  use nrtype
  implicit none
  type(DataStructure), intent(inout) :: HHData1

  ! variables used by E04NCA: solve quadratic program
  integer(i4b)               :: LCWSAV,LLWSAV,LIWSAV,LRWSAV,LDC
  integer(i4b), allocatable  :: IWSAV(:)
  real(dp),     allocatable  :: RWSAV(:)
  logical,      allocatable  :: LWSAV(:)
  character(6)               :: RNAME
  character(80), allocatable :: CWSAV(:)
  integer(i4b)               :: M,N,NCLIN,LDA
  integer(i4b), allocatable  :: istate(:),kx(:),iwork(:)
  integer(i4b)               :: iter,liwork,lwork,ifail,mode,i1
  real(dp), allocatable      :: C(:,:),BL(:),BU(:),CVEC(:),q(:),A(:,:),B(:)
  real(dp), allocatable      :: clamda(:),work(:)
  real(dp)                   :: obj
  real(dp)                   :: crit

  real(dp), allocatable      :: eye(:,:),zero(:)

  ! set options for E04NCA
  integer(i4b)       :: options_unit
  integer(i4b)       :: err_unit
  integer(i4b)       :: PriceFlag
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
  allocate(q(N))
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
      call ComputeCurrentB(HHData1%eta(:,i1),parms,HHData1%month(i1))
    end if

    CVEC = HHData1%p(:,i1) - matmul(transpose(parms%B),HHData1%e(:,i1))
    q = 0.0d0
    ifail = -1
    ! Compute solution to 
    !    min 0.5 * q'*A*q + cvec'*q   subject to q>=0
    !        0.5 * q' * (B'*B) * q + (p - B'*e)' * q
    !  = max  - 0.5 * (Bq-e)' * (Bq-e) - p'*q + 0.5 * e'*e
    ! E04NCA transforms A: so A needs to be reset to initial value

    A = matmul(transpose(parms%B),parms%B)
    ! minimise quadratic program
    call E04NCA(M,N,NCLIN,LDC,LDA,C,BL,BU,CVEC,ISTATE,KX,q,A,B,iter,OBJ,CLAMDA, &
                IWORK,LIWORK,WORK,LWORK,LWSAV,IWSAV,RWSAV,ifail)

    ! utility = -OBJ + constant
    HHData1%utility(i1)     = -OBJ
    HHData1%expenditure(i1) = dot_product(HHData1%p(:,i1),q)
    HHData1%nNonZero(i1)    = count(q>=crit)
    HHData1%iNonZero(1:HHData1%nNonZero(i1),i1) = pack((/1:parms%J/),q>=crit)
    HHData1%iZero(1:parms%J-HHData1%nNonZero(i1),i1) = pack((/1:parms%J/),q<crit)
    if (HHData1%nNonZero(i1)>0) then
      HHData1%q(1:HHData1%nNonZero(i1),i1) = pack(q,q>=crit)
    end if
  end do

  deallocate(IWSAV,RWSAV,LWSAV,CWSAV)
  deallocate(C,BL,BU,CVEC,istate,kx,q,A,B,clamda)
  deallocate(WORK,IWORK)

end subroutine ComputeDemand

!!   END OF ANALYSIS of DEMAND

subroutine ComputeInitialGuess(parms,iFree,x)
  use GlobalModule, only : ParmsStructure,SelectFreeType
  use nrtype
  implicit none
  type(ParmsStructure), intent(in) :: parms
  type(SelectFreeType), intent(in) :: iFree
  real(dp),             intent(out) :: x(:)
  real(dp), allocatable             :: temp1(:)

  x = 0.0d0
  if (parms%model==1) then
    if (allocated(iFree%D)) then
      x(iFree%xD) = parms%D(iFree%D)
    end if
    if (allocated(iFree%BC)) then
      x(iFree%xBC) = parms%BC(iFree%BC)
    end if
  end if

  if (allocated(iFree%MuE)) then
    x(iFree%xMUE) = parms%MuE(iFree%MUE)
  end if

  if (allocated(iFree%MuE_month)) then
    allocate(temp1(12*parms%K))
    temp1 = reshape(parms%mue_month,(/12*parms%K/))
    x(iFree%xMUE_month) = temp1(iFree%MUE_month)
    deallocate(temp1)
  end if
  if (allocated(iFree%InvCDiag)) then
    x(iFree%xInvCDiag) = parms%InvCDiag(iFree%InvCDiag)
  end if
  if (allocated(iFree%InvCOffDiag)) then
    x(iFree%xInvCOffDiag) = parms%InvCOffDiag(iFree%InvCOffDiag)
  end if

  if (parms%model==2) then
    if (allocated(iFree%BC_beta)) then
      x(iFree%xBC_beta) = parms%BC_beta(iFree%BC_beta)
    end if
    if (allocated(iFree%BC_CDiag)) then
      x(iFree%xBC_CDiag) = parms%BC_CDiag(iFree%BC_CDiag)
    end if
    if (allocated(iFree%BC_COffDiag)) then
      x(iFree%xBC_COffDiag) = parms%BC_COffDiag(iFree%BC_COffDiag)
    end if

    if (allocated(iFree%BD_beta)) then
      x(iFree%xBD_beta) = parms%BD_beta(iFree%BD_beta)
    end if
    if (allocated(iFree%BD_month)) then
      allocate(temp1(11*parms%J))
      temp1 = reshape(parms%BD_month(:,2:12),(/parms%J*11/))
      x(iFree%xBD_month) = temp1(iFree%BD_month)
    end if
    if (allocated(iFree%BD_CDiag)) then
      x(iFree%xBD_CDiag) = parms%BD_CDiag(iFree%BD_CDiag)
    end if
    if (allocated(iFree%BD_COffDiag)) then
      x(iFree%xBD_COffDiag) = parms%BD_COffDiag(iFree%BD_COffDiag)
    end if
  end if
 
end subroutine ComputeInitialGuess

subroutine SelectFreeParameters(parms,iFree)
  use nrtype
  use GlobalModule, only : ParmsStructure,SelectFreeType
  implicit none

  type(ParmsStructure), intent(in)    :: parms
  type(SelectFreeType), intent(inout) :: iFree

  integer(i4b) :: i1,fruit,month,ix

  !  FreeFlags.D           = 1;
  !  FreeFlags.C           = 1;
  !  FreeFlags.MuE         = 0;
  !  FreeFlags.InvCDiag    = 0;
  !  FreeFlags.InvCOffDiag = 0;

  ! parms%D  = (J x 1) scaling factors in B
  ! iFree.D  = elements of D that are free
  ! iFree%xD = elements of xFree corresponding to D
  iFree%nall = 0

  if (parms%model==2) then
    iFree%flagD=0
    iFree%flagBC=0
  end if

  if (iFree%flagD==1) then
    allocate(iFree%D(parms%J))
    allocate(iFree%xD(parms%J))
    iFree%D  = (/1:parms%J/)
    iFree%xD = iFree%D
    iFree%nD = parms%J
  else
    iFree%nD = 0
  end if
  iFree%nall = iFree%nD

  ! parms%BC  = (nBC x 1) spherical representation of normalized B
  ! iFree%xBC = elements of xFree corresponding to BC
  ! iFree%BC  = elements of BC that are free
  if (iFree%flagBC==1) then
    allocate(iFree%BC(parms%nBC))
    allocate(iFree%xBC(parms%nBC))
    iFree%xBC = iFree%nall + (/1:parms%nBC/)
    iFree%BC = (/1:parms%nBC/)
    iFree%nBC = parms%nBC
  else
    iFree%nBC = 0
  end if
  iFree%nall = iFree%nall + iFree%nBC

  ! parms%MuE  = (K x 1) mean of e
  ! iFree%xMuE = elements of xFree corresponding to MuE
  ! iFree.MuE  = elements of MuE that are free
  if (iFree%flagMUE==1) then
    allocate(iFree%MuE(parms%K-ifree%mue1+1))
    allocate(iFree%xMuE(parms%K-ifree%mue1+1))
    iFree%MuE = (/ifree%mue1:parms%K/)
    iFree%xMuE = iFree%nall + (/1:size(iFree%MuE)/)
    iFree%nMuE = size(iFree%MuE)
  else
    iFree%nMUE = 0
  end if
  iFree%nall = iFree%nall + iFree%nMUE

  ! parms%mue_month  = (K x 12) mean of e.   mue_month(:,1)=0.0d
  ! iFree%xMuE_month = elements of xFree corresponding to MuE_month
  ! iFree%MuE_month  = elements of MuE_month that are free
  if (iFree%flagMUE_month==1) then
    allocate(iFree%MuE_month(11*parms%K - ifree%mue_month1+1))
    allocate(iFree%xMuE_month(11*parms%K - ifree%mue_month1+1))
    iFree%nmue_month = size(ifree%mue_month)
    iFree%MuE_month  = (/(ix,ix=1,11*parms%K)/)
    iFree%xMuE_month = iFree%nall + (/(ix,ix=1,iFree%nmue_month)/)
  else
    iFree%nMUE_month = 0
  end if
  iFree%nall = iFree%nall + iFree%nMUE_month

  ! parms%InvCDiag  = (K x 1)  inverse of standard deviation Sig
  ! iFree%xInvCDiag = elements of xFree corresponding to InvCDiag
  ! iFree%InvCDiag  = elements of InvCDiag that are free
  if (iFree%flagInvCDiag==1) then
    allocate(iFree%InvCDiag(parms%K - ifree%invcdiag1+1))
    allocate(iFree%xInvCDiag(parms%K - ifree%invcdiag1+1))
    iFree%InvCDiag = (/(ix,ix=ifree%invcdiag1,parms%K)/)
    iFree%nInvCDiag = size(iFree%InvCDiag)
    iFree%xInvCDiag = iFree%nall + (/(ix,ix=1,iFree%nInvCDiag)/)
  else
    iFree%nInvCDiag = 0
  end if
  iFree%nall=iFree%nall+iFree%nInvCDiag

  ! parms%InvCOffDiag  = (K*(K-1)/2 x 1)  angles from spherical representation of 
  !                                    Cholesky decomposition of inverse of Sig 
  ! iFree%xInvCOffDiag = elements of xFree corresponding to InvCOffDiag
  ! iFree%InvCOffDiag  = elements of InvCOffDiag that are free
  if (iFree%flagInvCOffDiag==1) then
    iFree%nInvCOffDiag = parms%K*(parms%K-1)/2 - ifree%invcoffdiag1 + 1
    allocate(iFree%InvCOffDiag(iFree%nInvCOffDiag))
    allocate(iFree%xInvCOffDiag(iFree%nInvCOffDiag))
    iFree%InvCOffDiag  = (/(ix,ix=ifree%invcoffdiag1,(parms%k*(parms%k-1)/2))/)
    iFree%xInvCOffDiag = iFree%nall + (/(ix,ix=1,iFree%nInvCOffDiag)/)
  else
    iFree%nInvCOffDiag = 0
  end if
  iFree%nall = iFree%nall + iFree%nInvCOffDiag

  if (parms%model==2) then
    iFree%nBD_beta     = 0
    iFree%nBD_CDiag    = 0
    iFree%nBD_COffDiag = 0
    iFree%nBC_beta     = 0
    iFree%nBC_CDiag    = 0
    iFree%nBC_COffDiag = 0
    iFree%nBD_month    = 0 

    if (iFree%flagBD_beta==1) then
      iFree%nBD_beta = parms%BD_z_dim - ifree%bd_beta1 + 1
      allocate(iFree%BD_beta(iFree%nBD_beta))
      allocate(iFree%xBD_beta(iFree%nBD_beta))
      iFree%BD_beta = (/(ix,ix=ifree%bd_beta1,parms%BD_z_dim)/)
      iFree%xBD_beta = iFree%nAll + (/(ix,ix=1,ifree%nbd_beta)/)
      iFree%nall = iFree%nall + iFree%nBD_beta
    end if

    if (iFree%flagBC_beta==1) then
      iFree%nBC_beta = parms%BC_z_dim - ifree%bc_beta1 + 1
      allocate(iFree%BC_beta(iFree%nBC_beta))
      allocate(iFree%xBC_beta(iFree%nBC_beta))
      iFree%BC_beta = (/(ix,ix=ifree%bc_beta1,parms%BC_z_dim)/)
      iFree%xBC_beta = iFree%nAll + (/(ix,ix=1,ifree%nbc_beta)/)
      iFree%nall = iFree%nall + iFree%nBC_beta
    end if

    if (iFree%flagBD_CDiag==1) then
      iFree%nBD_CDiag = parms%J - ifree%bd_cdiag1 + 1
      allocate(iFree%BD_CDiag(iFree%nBD_CDiag))
      allocate(iFree%xBD_CDiag(iFree%nBD_CDiag))
      iFree%BD_CDiag = (/(ix,ix=ifree%bd_cdiag1,parms%J)/)
      iFree%xBD_CDiag = iFree%nAll + (/(ix,ix=1,ifree%nbd_cdiag)/)
      iFree%nall = iFree%nall + iFree%nBD_CDiag
    end if

    if (iFree%flagBD_COffDiag==1) then
      iFree%nBD_COffDiag = size(parms%BD_COffDiag,1) - ifree%bd_coffdiag1 + 1
      allocate(iFree%BD_COffDiag(iFree%nBD_COffDiag))
      allocate(iFree%xBD_COffDiag(iFree%nBD_COffDiag))
      iFree%BD_COffDiag = (/(ix,ix=ifree%bd_coffdiag1,size(parms%bd_coffdiag))/)
      iFree%xBD_COffDiag = iFree%nAll + (/(ix,ix=1,ifree%nbd_coffdiag)/)
      iFree%nall = iFree%nall + iFree%nBD_COffDiag
    end if

    if (iFree%flagBC_CDiag==1) then
      iFree%nBC_CDiag = parms%nBC - ifree%bc_cdiag1 + 1
      allocate(iFree%BC_CDiag(iFree%nBC_CDiag))
      allocate(iFree%xBC_CDiag(iFree%nBC_CDiag))
      iFree%BC_CDiag = (/(ix,ix=ifree%bc_cdiag1,parms%nbc)/)
      iFree%xBC_CDiag = iFree%nAll + (/(ix,ix=1,ifree%nbc_cdiag)/)
      iFree%nall = iFree%nall + iFree%nBC_CDiag
    end if

    if (iFree%flagBC_COffDiag==1) then
      iFree%nBC_COffDiag = size(parms%BC_COffDiag,1) - ifree%bc_coffdiag1 + 1
      allocate(iFree%BC_COffDiag(iFree%nBC_COffDiag))
      allocate(iFree%xBC_COffDiag(iFree%nBC_COffDiag))
      iFree%BC_COffDiag = (/(ix,ix=ifree%bc_coffdiag1,size(parms%bc_coffdiag))/)
      iFree%xBC_COffDiag = iFree%nAll + (/(ix,ix=1,size(parms%bc_coffdiag))/)
      iFree%nall = iFree%nall + iFree%nBC_COffDiag
    end if

    ! MONTH 1 is ALWAYS BASE CATEGORY
    if (iFree%flagBD_month==1) then
      iFree%nBD_month = 11* parms%J - ifree%bd_month1 + 1
      allocate(iFree%BD_month(iFree%nBD_month))
      allocate(iFree%xBD_month(iFree%nBD_month))
      iFree%BD_month  = (/(ix,ix=ifree%bd_month1,11*parms%j)/)
      iFree%xBD_month = iFree%nAll + (/(ix,ix=1,ifree%nbd_month)/)
      iFree%nall = iFree%nall + iFree%nBD_month
    end if
  end if ! if (parms%model==2)

  !  create labels for x variables
  allocate(iFree%xlabels(iFree%nall))
  do i1=1,iFree%nD
    write(iFree%xlabels(iFree%xD(i1)),'(a2,i2.2,a1)') 'D(', iFree%D(i1),')'
  end do

  do i1=1,iFree%nBC
    write(iFree%xlabels(iFree%xBC(i1)),'(a3,i2.2,a1)') 'BC(', iFree%BC(i1),')'
  end do

  do i1=1,iFree%nInvCDiag
    write(iFree%xlabels(iFree%xInvCDiag(i1)),'(a9,i2.2,a1)') 'InvCDiag(', iFree%InvCDiag(i1), ')'
  end do

  do i1=1,iFree%nInvCOffDiag
    write(iFree%xlabels(iFree%xInvCOffDiag(i1)),'(a12,i2.2,a1)') 'InvCOffDiag(', iFree%InvCOffDiag(i1), ')'
  end do

  do i1=1,iFree%nMuE
    write(iFree%xlabels(iFree%xMUE(i1)),'(a4,i2.2,a1)') 'MUE(', iFree%MUE(i1), ')'
  end do

  do i1=1,iFree%nBD_BETA
    write(iFree%xlabels(iFree%xBD_BETA(i1)),'(a8,i2.2,a1)') 'BD_beta(', iFree%BD_BETA(i1), ')'
  end do

  do i1=1,iFree%nBD_CDiag
    write(iFree%xlabels(iFree%xBD_CDiag(i1)),'(a9,i2.2,a1)') 'BD_CDiag(', iFree%BD_CDiag(i1), ')'
  end do

  do i1=1,iFree%nBD_COffDiag
    write(iFree%xlabels(iFree%xBD_COffDiag(i1)),'(a12,i2.2,a1)') 'BD_COffDiag(', iFree%BD_COffDiag(i1), ')'
  end do

  do i1=1,iFree%nBC_BETA
    write(iFree%xlabels(iFree%xBC_BETA(i1)),'(a8,i2.2,a1)') 'BC_beta(', iFree%BC_BETA(i1), ')'
  end do

  do i1=1,iFree%nBC_CDiag
    write(iFree%xlabels(iFree%xBC_CDiag(i1)),'(a9,i2.2,a1)') 'BC_CDiag(', iFree%BC_CDiag(i1), ')'
  end do

  do i1=1,iFree%nBC_COffDiag
    write(iFree%xlabels(iFree%xBC_COffDiag(i1)),'(a12,i2.2,a1)') 'BC_COffDiag(', iFree%BC_COffDiag(i1), ')'
  end do

  do i1=1,iFree%nmue_month
    !                                             'mue_month(fruit,month)'
    fruit = (i1-1)/11 + 1
    month = i1 - 11* ((i1-1)/11)
    write(iFree%xLabels(iFree%xmue_month(i1)),'(a10,i2.2,",",i2.2,a1)') 'mue_month(',fruit,month,')'
  end do

  do i1=1,iFree%nBD_month
    !                                             'BD_month(fruit,month)'
    fruit = (i1-1)/11 + 1
    month = i1 - 11* ((i1-1)/11)
    write(iFree%xLabels(iFree%xBD_month(i1)),'(a9,i2.2,",",i2.2,a1)') 'BD_month(',fruit,month,')'
  end do
end subroutine SelectFreeParameters

subroutine MaximizeLikelihood(x,LValue,Grad,Hess,ierr)
  use nrtype
  use GlobalModule, only : MaxOptions,Penalty
  implicit none
  real(dp),     intent(inout) :: x(:)
  real(dp),     intent(out)   :: LValue
  real(dp),     intent(out)   :: Grad(:)
  real(dp),     intent(out)   :: Hess(:,:)
  integer(i4b), intent(out)   :: ierr

  if (MaxOptions%Algorithm==1) then
    ! E04WDF:  Dense optimization: constrained non-linear optimization for dense problem
    if (Penalty%method==0) then
      call MaximizeLikelihood1(x,LValue,Grad,Hess,ierr)
    else if (Penalty%method==1) then
      call MaximizePenalizedLikelihood(x,LValue,Grad,Hess,ierr)
    end if
  elseif (MaxOptions%Algorithm==2) then
    ! E04VHF:  Sparse optimization: constrained non-linear optimization for sparse problem
!    call MaximizeLikelihood2(x,LValue,Grad,Hess,ierr)
  elseif (MaxOptions%Algorithm==3) then
     call Max_E04JCF(x,LValue,Grad,Hess,ierr)
  elseif (MaxOptions%Algorithm==6) then
    ! Bayesian estimation
    call BayesLikelihood1(x,LValue,Grad,Hess,ierr)
    !call ComputeBayes(x,LValue,Grad,Hess,ierr)
  end if
end subroutine MaximizeLikelihood

subroutine ComputeBayes(x,L,Grad,Hess,ierr)
  use nrtype
  use GlobalModule, only : bayes

  implicit none
  real(dp), intent(in) :: x(:)
  real(dp), intent(out) :: L,Grad(:),Hess(:,:)
  integer(i4b), intent(out) :: ierr

  integer(i4b)  :: mode,nx
  integer(i4b)  :: nstate,i1,i2
  real(dp), allocatable :: GradL(:)
  integer(i4b) :: iuser(2)
  real(dp)     :: ruser(1)

  ! moments of x: zero, one and two
  real(dp)              :: m0
  real(dp), allocatable :: m1(:),m2(:,:)
  
  mode=0
  nx = size(x)
  allocate(m1(nx))
  allocate(m2(nx,nx))

  m0 = 0.0d0
  m1 = 0.0d0
  m2 = 0.0d0

  do i1=1,bayes%nAll
    call LikeFunc(mode,nx,bayes%x(:,i1),L,GradL,nstate,iuser,ruser)
    m0 = m0 + bayes%w(i1) * bayes%prior(i1) * L
    m1 = m1 + bayes%w(i1) * bayes%prior(i1) * L * bayes%x(:,i1)
    do i2=1,nx
      m2(:,i2) = m2(:,i2) + bayes%w(i1) * bayes%prior(i1) * L  &
                          * bayes%x(:,i1) * bayes%x(i2,i1)
    end do
  end do

  L    = m0
  Grad = m1
  Hess = m2
  deallocate(m1,m2)
  deallocate(bayes%x,bayes%w,bayes%prior)
  ierr = 0
end subroutine ComputeBayes

! Maximise using Constrained maximization: E04WDF
subroutine MaximizeLikelihood1(x,LValue,Grad,Hess,ierr)
  use nrtype
#if USE_MPI==1
  use mpi
#endif
  use GlobalModule, only : ControlOptions,parms,InputDir,OutDir,MasterID,HHData,iFree, &
                           MaxOptions,ResultStructure,ReadWriteParameters
  use nag_library,  only : E04WCF,E04WDF,E04WDP,E04WEF,E04WFF, &
                           X04AAF,X04ABF,X04ACF,X04ADF
  use OutputModule, only : ComputeStats, SaveOutputs

  implicit none
  real(dp),     intent(inout) :: x(:)
  real(dp),     intent(out)   :: LValue
  real(dp),     intent(out)   :: Grad(:)
  real(dp),     intent(out)   :: Hess(:,:)
  integer(i4b), intent(out)   :: ierr

  ! Default values of parameters required for NAG routine e04WDF
  integer(i4b), parameter     :: leniw=600,lenrw=600
  integer(i4b), parameter     :: eflag=1
  integer(i4b)                :: eunit

  ! used to transmit completion message to worker in MPI version
  integer(i4b)                :: task

  ! Dimensions of optimization problem
  integer(i4b)                :: nx,nc_lin,nc_nonlin,ldcj,ldh,lda,nctotal

  ! Loop counter
  integer(i4b)                :: i1,ix

  ! initial value of x and likelihood function
  real(dp), allocatable       :: x0(:)
  real(dp)                    :: LValue0

  ! statistics from MLE
  type(ResultStructure)      :: LikeStats

  ! linear constraints on optimization problem
  real(dp), allocatable       :: A(:,:)

  ! LOWER AND UPPER BOUNDS for optimisation
  real(dp), allocatable       :: BL(:),BU(:)            

  ! MULTIPLIERS ON CONSTRAINTS
  REAL(DP), allocatable       :: CLAMBDA(:)                          

  ! PARAMETERS FOR E04WCF and E04WDF (constrained optimisation subroutines from NAG)
  INTEGER(I4B)                :: IFAIL,iter,nstate
  integer(i4b), allocatable   :: ISTATE(:)
  real(dp), allocatable       :: ccon(:)
  real(dp), allocatable       :: cjac(:,:)
  integer(i4b)                :: E04Unit   ! unit number for E04WDF options file
  character(99)               :: E04File   ! file name for E04WDF options file
  character(99)               :: OutFile   ! file name for E04WDF options file
  character(20)               :: SaveFreqString
 
! Workspace for E04WCF and E04WDF (constrained optimisation subroutines from NAG)
  real(dp)                    :: RW(lenrw)
  integer(i4b)                :: IW(leniw)
  
  ! E04WDF outputs
  !real(DP)                  :: OBJF

  ! User INPUTS TO OBJECTIVE FUNCTION
  REAL(DP)                    :: RUSER(1)
  INTEGER(I4B)                :: iuser(2)

  ! used to test NAGConstraintWrapper
  integer(i4b)                :: mode_constraint
  integer(i4b), allocatable   :: needc(:)
    
  ! Dimensions of maximization problem
  nx=size(x,1)    
  allocate(x0(nx))
  x0 = x

  iuser(1) = parms%model
  iuser(2) = HHData%N
  ruser    = 0.0d0

#if USE_MPI==1
  ! broadcast (nx,iuser,ruser) 
  call mpi_bcast(nx,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  ! broadcast length of iuser
  call mpi_bcast(2,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iuser,2,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  ! broadcast length of ruser
  call mpi_bcast(1,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ruser,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr)
#endif

  nc_lin=0                 ! number of linear constraints
  if (ControlOptions%TestLikeFlag<4) then
    ! normal non-linear constraints
    call ComputeNC(nc_nonlin)   ! number of nonlinear constraints
  else if (ControlOptions%TestLikeFlag>=4) then
    ! no non-linear constraints
    nc_nonlin=0
  end if
  lda=max(1,nc_lin)
  ldcj=max(1,nc_nonlin)
  ldh=size(HESS,1)
  nctotal=nx+nc_lin+nc_nonlin

  ! Allocate memory for linear constraints, bound constraints, and some aux. outputs
  allocate(A(lda,1))
  allocate(BL(nctotal),BU(nctotal))
  allocate(CLAMBDA(nctotal),ISTATE(nctotal))
  if (nc_nonlin==0) then
    allocate(CJAC(LDCJ,1))
  else
    allocate(CJAC(LDCJ,nx))
  end if
  allocate(CCON(max(1,nc_nonlin)))

  ! linear constraints:   BL <= A*x <= BU
  !     NONE
  !  A is not referenced when nclin==0
  ! lower and upper bounds on xFree 
  BL = -5.0d0
  !-huge(1.0d0)
  BU = 5.0d0
  !huge(1.0d0)
  call SetBounds(x,BL(1:nx),BU(1:nx))

  ! Linear constraints
  !   none

  ! Nonlinear constraints
  if (nc_nonlin>0) then
    BU(nx+nc_lin+(/(ix,ix=1,nc_nonlin)/)) = 0.0d0
  end if

  nstate = 1
  if (ControlOptions%TestLikeFlag==1) then
    ! test gradient of likelihood
    call TestGrad(LikeFunc,nx,x,nstate,iuser,ruser)
  else if (ControlOptions%TestLikeFlag==2) then
    ! plot likelihood function 
    call PlotLike(LikeFunc,nx,x,iFree%xlabels,BL(1:nx),BU(1:nx),nstate,iuser,ruser)
  else if (ControlOptions%TestLikeFlag==0 .or. ControlOptions%TestLikeFlag>2) then 
    ! Initialize E04WDF by calling E04WCF
    ifail=-1
    call e04wcf(IW,LENIW,RW,LENRW,ifail)

    ! Suppress NAG error messages
    eunit = -1
    call x04aaf(eflag,eunit)
    ! set standard output for advisory message
    eunit = 6
    CALL X04ABF(1,eunit)
    ! set options for E04WDF
    E04Unit = 50
    E04File = MaxOptions%OptionsFile
    open(UNIT = E04Unit,  &
         FILE = E04File, &
         ACTION = 'read')
    ifail=-1
    call E04WEF(E04Unit,IW,RW,ifail)
    close(E04Unit)

    ! open files for Backup Basis File and New Basis File
    if (MaxOptions%SaveBasis==1) then
        open(unit = 101, &
             file = MaxOptions%BasisFile, &
             action='WRITE')
        open(unit=103, &
             file=MaxOptions%BackupBasisFile, &
             action='WRITE')
        call E04WFF('New Basis File = 101',iw,rw,ifail)
        call E04WFF('Backup Basis File = 103',iw,rw,ifail)
        write(SaveFreqString,2484) 'Save Frequency = ',MaxOptions%SaveBasisFreq
2484 format(a17,i3)
        call E04WFF(SaveFreqString,iw,rw,ifail)
    end if
    ! open OldBasisFile to load information on basis
    if (MaxOptions%LoadBasis==1) then
      open(unit=105, &
           file = MaxOptions%OldBasisFile, &
           action='READ')
      call E04WFF('Old Basis File = 105',iw,rw,ifail)
    end if

    LValue = 0.0D0
    Grad   = 0.0d0
    Hess   = 0.0d0
    do iter=1,nx
      Hess(iter,iter)=1.0d0
    end do

    if (ControlOptions%TestLikeFlag==0) then
       ! Maximise likelihood 
       print *,'Begin maximization.'
       if (ControlOptions%OutputFlag >0) then
          call ComputeHess(x0,LValue0,GRAD,Hess,iuser,ruser)
       end if
       print *,'LValue0',LValue0
       ifail = -1
       call E04WDF(nx,nc_lin,nc_nonlin,LDA,LDCJ,LDH,A,BL,BU,                     &
                   NAGConstraintWrapper,LikeFunc,iter,ISTATE,CCON,CJAC,CLAMBDA,  &
                   LValue,GRAD,HESS,x,IW,LENIW,RW,LENRW,iuser,RUSER,ifail)
       ! update parms
       call UpdateParms2(x,iFree,parms)
       ! b) save parms to disk
       call ReadWriteParameters(parms,'write')
       print *,'E04WDF complete.'

       if (ControlOptions%OutputFlag ==2) then
         call ComputeHess(x0,LValue0,GRAD,Hess,iuser,ruser)
         call ComputeHess(x,LValue,GRAD,Hess,iuser,ruser)
       else if (ControlOptions%OutputFlag==1) then
         call ComputeHess2(x,LValue,Grad,Hess,ControlOptions%ComputeHessFlag)
       end if
    else if (ControlOptions%TestLikeFlag==3) then
      ! test non-linear contraint
      mode_constraint=0
      allocate(needc(nc_nonlin))
      needc = (/(ix,ix=1,nc_nonlin)/)
      call NAGConstraintWrapper(mode_constraint,nc_nonlin,nx,LDCJ,NEEDC,X,CCON,CJAC,nstate,iuser,ruser)
      deallocate(needc)
    else if (ControlOptions%TestLikeFlag==4) then
       ! maximise likelihood with no non-linear constraints
       print *,'Begin maximization with no non-linear constraints.'
       if (ControlOptions%OutputFlag > 0) then
         call ComputeHess(x0,LValue0,GRAD,Hess,iuser,ruser)
       end if
       print *,'LValue0',LValue0
       ifail = -1
       call E04WDF(nx,nc_lin,nc_nonlin,LDA,LDCJ,LDH,A,BL,BU,                     &
                   E04WDP,LikeFunc,iter,ISTATE,CCON,CJAC,CLAMBDA,                &
                   LValue,GRAD,HESS,x,IW,LENIW,RW,LENRW,iuser,RUSER,ifail)
       ! update parms
       call UpdateParms2(x,iFree,parms)
       ! b) save parms to disk
       call ReadWriteParameters(parms,'write')
       print *,'E04WDF complete.'
       if (ControlOptions%OutputFlag ==2) then
         call ComputeHess(x0,LValue0,GRAD,Hess,iuser,ruser)
         call ComputeHess(x,LValue,GRAD,Hess,iuser,ruser)
       else if (ControlOptions%OutputFlag==1) then
         call ComputeHess2(x,LValue,Grad,Hess,ControlOptions%ComputeHessFlag)
       end if
    else if (ControlOptions%TestLikeFlag==5) then
      ! Compute Hess only
      call ComputeHess2(x,LValue,Grad,Hess,ControlOptions%ComputeHessFlag)
    end if

    if (ControlOptions%OutputFlag==1) then
      call ComputeStats(x,LValue,Grad,Hess,HHData%N,LikeStats)
      call SaveOutputs(x,LValue,Grad,Hess,LikeStats)
    else 
      OutFile = trim(OutDir) // '/results.txt'
      open(unit = 130,file = OutFile,action = 'write')
      write(130,'(a20,3a25)') 'Var. Name','x0','x','Gradient'
      do i1=1,nx
        write(130,'(a20,3d25.12)') trim(iFree%xlabels(i1)), x0(i1),x(i1), Grad(i1)
      end do
      write(130,'(a25,d25.12)') 'LValue0',LValue0
      write(130,'(a25,d25.12)') 'LValue',LValue
      close(130)
    end if
    if (MaxOptions%SaveBasis==1) then
      close(101)
      close(103)
    end if
    if (MaxOptions%LoadBasis==1) then
      close(105)
    end if
 end if  ! if (ControlOptions%TestLikeFlag==1) then

  deallocate(A,BL,BU,CLAMBDA,ISTATE,CJAC)
  deallocate(CCON)

#if USE_MPI==1
  ! broadcast signal indicating that all mpi tasks are complete
  task=0
  call mpi_bcast(task,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
#endif
  ierr = 0   ! no error in subroutine
end subroutine MaximizeLikelihood1

! compute Bayes estimator using D01ESF
subroutine BayesLikelihood1(x,LValue,Grad,Hess,ierr)
  use nrtype
#if USE_MPI==1
  use mpi
#endif
  use GlobalModule, only : parms,HHData,iFree,MasterID

  implicit none
  real(dp),     intent(inout) :: x(:)
  real(dp),     intent(out)   :: LValue
  real(dp),     intent(out)   :: Grad(:)
  real(dp),     intent(out)   :: Hess(:,:)
  integer(i4b), intent(out)   :: ierr

  ! Dimensions of optimization problem
  integer(i4b)                :: nx,iter,task

  ! LOWER AND UPPER BOUNDS for optimisation
  real(dp), allocatable       :: BL(:),BU(:)            

  ! User INPUTS TO OBJECTIVE FUNCTION
  REAL(DP)                    :: RUSER(1)
  INTEGER(I4B)                :: iuser(2)

  ! Dimensions of maximization problem
  nx=size(x,1)    

  allocate(BL(nx),BU(nx))

  BL = -5.0d0
  BU = 5.0d0
  call SetBounds(x,BL(1:nx),BU(1:nx))

  iuser(1) = parms%model
  iuser(2) = HHData%N
  ruser    = 0.0d0

  LValue = 0.0D0
  Grad   = 0.0d0
  Hess   = 0.0d0
  do iter=1,nx
    Hess(iter,iter)=1.0d0
  end do

#if USE_MPI==1
  ! broadcast (nx,iuser,ruser)
  call mpi_bcast(nx,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(2,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iuser,2,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(1,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ruser,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr)
#endif
  call SparseGridBayes(iFree,iuser,ruser)

  deallocate(BL,BU)

#if USE_MPI==1
  task=0
  call mpi_bcast(task,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
#endif
  ierr = 0   ! no error in subroutine
end subroutine BayesLikelihood1

! Maximise using Constrained maximization: E04WDF
subroutine MaximizePenalizedLikelihood(x,LValue,Grad,Hess,ierr)
  use nrtype
  use GlobalModule, only : Penalty,ControlOptions,MaxOptions,InputDir
  use nag_library, only : E04WCF,E04WDF,E04WDP,E04WEF,E04WFF, &
                          X04AAF,X04ABF,X04ACF,X04ADF

  implicit none
  real(dp), intent(inout)  :: x(:)
  real(dp), intent(out)    :: LValue
  real(dp), intent(out)    :: Grad(:)
  real(dp), intent(out)    :: Hess(:,:)
  integer(i4b), intent(out) :: ierr

  ! xP  = [x;xPlus;xMinus]
  real(dp), allocatable :: xP(:)
  real(dp), allocatable :: GradP(:)
  real(dp), allocatable :: HessP(:,:)

  ! Default values of parameters required for NAG routine e04WDF
  integer(i4b), parameter   :: leniw=600,lenrw=600
  integer(i4b), parameter   :: eflag=1
  integer(i4b)              :: eunit
  integer(i4b)              :: E04Unit 
  character(len=30)         :: E04File

  ! used to transmit completion message to worker in MPI version
  integer(i4b)              :: mode

  integer(i4b)              :: i1
  ! Dimensions of optimization problem
  integer(i4b)              :: nxP,nclin,ncnln,ldcj,ldh,lda,nctotal

  ! linear constraints on optimization problem
  real(dp), allocatable     :: A(:,:)

  ! LOWER AND UPPER BOUNDS for optimisation
  real(dp), allocatable     :: BL(:),BU(:)            

  ! MULTIPLIERS ON CONSTRAINTS
  REAL(DP), allocatable     :: CLAMBDA(:)                          

  ! PARAMETERS FOR E04WCF and E04WDF (constrained optimisation subroutines from NAG)
  INTEGER(I4B)              :: IFAIL,iter
  integer(i4b), allocatable :: ISTATE(:)
  real(dp), allocatable     :: ccon(:)
  real(dp), allocatable     :: cjac(:,:)
  
  ! Workspace for E04WCF and E04WDF (constrained optimisation subroutines from NAG)
  real(dp)                  :: RW(lenrw)
  integer(i4b)              :: IW(leniw)
  
  ! E04WDF outputs
  real(DP)                  :: OBJF

  ! User INPUTS TO OBJECTIVE FUNCTION
  REAL(DP)                  :: RUSER(1)
  INTEGER(I4B)              :: iuser(1)

  ! Control parameter for saving intermediate results
  !   0 DO NOT SAVE intermediate results
  !   1 SAVE intermediate results
  integer(i4b), parameter   :: NAG_SAVE=0    
  
  ! Dimensions of maximization problem
  nxP   = Penalty%nxP
  nclin = Penalty%nx2
  ncnln = 0
  lda   = max(1,NCLIN)
  ldcj  = max(1,NCNLN)
  ldh   = nxP
  nctotal=nxP+nclin+ncnln

  ! Allocate memory for linear constraints, bound constraints, and some aux. outputs
  allocate(A(lda,nxP))
  allocate(BL(nctotal),BU(nctotal))
  allocate(CLAMBDA(nctotal),ISTATE(nctotal))
  allocate(CJAC(LDCJ,1))
  allocate(CCON(max(1,NCNLN)))
  allocate(xP(nxP))
  allocate(GradP(nxP))
  allocate(HessP(nxP,nxP))

  ! Initialize E04WDF by calling E04WCF
  ifail=-1
  call e04wcf(IW,LENIW,RW,LENRW,ifail)

  ! Suppress NAG error messages
  eunit=-1
  call x04aaf(eflag,eunit)

  ! constraints:   BL <= A*xP <= BU
  !  A is not referenced when nclin==0
  !    xP = [x1;x2;xPlus;xMinus]
  !    0 <= x2 - xPlus+xMinus <=0
  ! size(A) = nx2 x nxP
  A = 0.0d0
  do i1=1,nclin
    A(i1,Penalty%xP_index2(i1)) = 1.0d0
    A(i1,Penalty%xP_index_xPlus(i1)) = -1.0d0
    A(i1,Penalty%xP_index_xMinus(i1)) = 1.0d0
  end do  

  ! lower and upper bounds on xP
  ! x1L  <= x1              <= x1H
  ! x2L  <= x2              <= x2H
  ! 0    <= xPlus           <= x2H
  ! 0    <= xMinus          <= -x2L
  ! 0    <= x2-xPlus+xMinus <= 0
  BL = 0.0d0
  BU = 0.0d0

  ! set bounds on x = (x1,x2)
  call SetBounds(xP(1:Penalty%nx),BL(1:Penalty%nx),BU(1:Penalty%nx))

  ! set bounds on (xPlus,xMinus)
  BU(Penalty%xP_index_xPlus)  = BU(Penalty%xP_index2)
  BU(Penalty%xP_index_xMinus) = -BL(Penalty%xP_index2)
 
  ! set options for E04WDF
  E04Unit = 50
  E04File = trim(InputDir) // '/E04WDF.opt'
  open(UNIT = E04Unit,  &
       FILE = E04File, &
       ACTION = 'read')
  ifail=-1
  call E04WEF(E04Unit,IW,RW,ifail)
  close(E04Unit)

  if (NAG_SAVE==1) then  
      open(101,FILE = 'output/NewBasis1.txt',ACTION='WRITE')
      open(103,FILE = 'output/BackupBasis1.txt',ACTION='WRITE')
      call E04WFF('New Basis File = 101',iw,rw,ifail)   ! unit number to save intermediate results
      call E04WFF('Backup Basis File = 103',iw,rw,ifail)   ! unit number to save intermediate results
      call E04WFF('Save Frequency = 100',iw,rw,ifail)
  end if

  ifail=-1
  iuser(1) = 1
  ruser    = 0.0d0
 
  LValue = 0.0D0
  Grad   = 0.0d0
  Hess   = 0.0d0
  GradP  = 0.0d0
  HessP  = 0.0d0
  do iter=1,nxP
    HessP(iter,iter)=1.0d0
  end do
  Hess = HessP(1:Penalty%nx,1:Penalty%nx)

  ! Initial guess for xP
  xP=0.0d0
  if (Penalty%nx1>0) then
    xP(Penalty%xP_index1) = x(Penalty%x_index1)
  end if
  xP(Penalty%xP_index2) = x(Penalty%x_index2)
  do i1=1,Penalty%nx2
    xP(Penalty%xP_index_xPlus(i1)) = max(x(Penalty%x_index2(i1)),0.0d0)
    xP(Penalty%xP_index_xMinus(i1)) = max(-x(Penalty%x_index2(i1)),0.0d0)
  end do

 ! Maximise likelihood 
  if (ControlOptions%TestLikeFlag==1) then
    call TestGrad(LikeFunc,nxP,xP,1,iuser,ruser)
  else
    print *,'Begin maximization.'
    call E04WDF(nxP,NCLIN,NCNLN,LDA,LDCJ,LDH,A,BL,BU,                &
                E04WDP,LikeFunc,iter,ISTATE,CCON,CJAC,CLAMBDA,     &
                LValue,GRADP,HESSP,xP,IW,LENIW,RW,LENRW,iuser,RUSER, &
                IFAIL)
    !call ComputeHess(xP(1:Penalty%nx),LValue,GRAD,Hess,iuser,ruser)
    if (Penalty%nx1>0) then
      x(Penalty%x_index1)    = xP(Penalty%xP_index1)
      Grad(Penalty%x_index1) = GradP(Penalty%xP_index1)
      Hess(Penalty%x_index1,Penalty%x_index1) = HessP(Penalty%xP_index1,Penalty%xP_index1)
      Hess(Penalty%x_index1,Penalty%x_index2) = HessP(Penalty%xP_index1,Penalty%xP_index2)
      Hess(Penalty%x_index2,Penalty%x_index1) = HessP(Penalty%xP_index2,Penalty%xP_index1)
    end if
    x(Penalty%x_index2) = xP(Penalty%xP_index2)
    Grad(Penalty%x_index2) = GradP(Penalty%xP_index2)
    Hess(Penalty%x_index2,Penalty%x_index2) = HessP(Penalty%xP_index2,Penalty%xP_index2)

  end if

  deallocate(A,BL,BU,CLAMBDA,ISTATE,CJAC)
  deallocate(CCON)
  deallocate(xP)
  if (NAG_SAVE==1) then
      close(101)
      close(103)
  end if

  ierr = 0   ! no error in subroutine
end subroutine MaximizePenalizedLikelihood

!------------------------------------------------------------------------------
! subroutine LikeFunc
!
!------------------------------------------------------------------------------
subroutine LikeFunc(mode,nx,x,L,GradL,nstate,iuser,ruser)
  use nrtype
  use GlobalModule, only : Penalty
  implicit none
  integer(i4b), intent(inout) :: mode 
  integer(i4b), intent(in)    :: nx,nstate
  integer(i4b), intent(inout) :: iuser(*)
  real(dp), intent(in)        :: x(nx)
  real(dp), intent(inout)     :: ruser(*)
  real(dp), intent(out)       :: L
  real(dp), intent(inout)     :: GradL(nx)
  
  integer(i4b)                :: model,N

  model = iuser(1)
  N     = iuser(2)  ! sample size
#if USE_MPI==0
  if (model==1) then
    ! Logit model
    if (Penalty%method==0) then
      call Like1(mode,nx,x,L,GradL,nstate,iuser,ruser)
    else
      ! Penalized likelihood
      call PenalizedLikelihood(mode,nx,x,L,GradL,nstate,iuser,ruser)
    end if
  elseif (model==2) then
    call Like2(mode,nx,x,L,GradL,nstate,iuser,ruser)
    L = L/real(N,dp)
    if (mode>0) then
    !  GradL = GradL/real(N,dp)
    end if
  end if
#elif USE_MPI==1
  if (model==1) then
 !   call Like1_master(mode,nx,x,L,GradL,nstate,iuser,ruser)
  elseif (model==2) then
    call Like2_master(mode,nx,x,L,GradL,nstate,iuser,ruser)
    L = L/real(N,dp)
    if (mode>0) then
      ! GradL = GradL/real(N,dp)
    end if  
  end if
#endif
end subroutine LikeFunc
 
! PenalizedLikelihood:
!  compute objective function for penalized likelihood for model == iuser(1)
subroutine PenalizedLikelihood(mode,nxP,xP,LP,GradLP,nstate,iuser,ruser)
  use nrtype
  use GlobalModule, only : Penalty,HHData
  implicit none
! xP       =  (nxP x 1) free parameters
!
! Revision history
! 19sep2012 LN   created subroutine 

  integer(i4b), intent(inout)        :: mode
  integer(i4b), intent(in)           :: nxP
  real(dp),     intent(in)           :: xP(:)
  real(dp),     intent(out)          :: LP
  real(dp),     intent(inout)        :: GradLP(:)
  integer(i4b), intent(in)           :: nstate
  integer(i4b), intent(in)           :: iuser(*)
  real(dp),     intent(in)           :: ruser(*)

  ! Method 1 : max L(x1,x2) - P(x2)    subject to      x2 = x2Plus - x2Minus
  !                                                    x2Plus >= 0
  !                                                    x2Minus >= 0
  !                                                    P(x2) = LASSO penalty
  real(dp)                           :: x(Penalty%nx)
  real(dp)                           :: xPlus(Penalty%nx2)
  real(dp)                           :: xMinus(Penalty%nx2)
  real(dp)                           :: L,GradL(Penalty%nx)
  real(dp)                           :: P,GradP(Penalty%nxP)
  integer(i4b)                       :: model
  model = iuser(1)

  x = 0.0d0
  if (Penalty%nx1>0) then
    x(Penalty%x_index1) = xP(Penalty%xP_index1)
  end if

  if (Penalty%nx2>0) then
    x(Penalty%x_index2) = xP(Penalty%xP_index2)
  end if
  
  if (model==1) then 
    ! Logit model
    call Like1(mode,Penalty%nx,x,L,GradL,nstate,iuser,ruser)
  end if

  xPlus  = xP(Penalty%xP_index_xPlus)
  xMinus = xP(Penalty%xP_index_xMinus)

  call PenaltyFunction(xPlus,xMinus,P,GradP)

  LP = L + real(HHData%N)*P
  if (mode>0) then
    GradLP = 0.d0
    ! elements of GradLP corresponding to x1
    if (Penalty%nx1>0) then
      GradLP(Penalty%xP_index1) = GradL(Penalty%x_index1)
    end if

    ! elements of GradLP corresponding to x2
    GradLP(Penalty%xP_index2) = GradL(Penalty%x_index2)

    ! elements of GradLP corresponding to xPlus
    ! elements of GradLP corresponding to xMinus
    GradLP(Penalty%xP_index_xPlus)  = real(HHData%N)*GradP(Penalty%xPlus_index)
    GradLP(Penalty%xP_index_xMinus) = real(HHData%N)*GradP(Penalty%xMinus_index)
  end if
end subroutine PenalizedLikelihood

subroutine PenaltyFunction(xPlus,xMinus,P,GradP)
  use nrtype
  use GlobalModule, only : Penalty
  implicit none
  real(dp), intent(in) :: xPlus(:)
  real(dp), intent(in) :: xMinus(:)
  real(dp), intent(out) :: P
  real(dp), intent(out) :: GradP(:)

  if (Penalty%method==1) then
    ! LASSO penalty
    !P = Penalty%lambda*(sum(abs(xPlus)) + sum(abs(xMinus)))
    !GradP(Penalty%xPlus_index) = Penalty%lambda*sign(1.0d0,xPlus)
    !GradP(Penalty%xMinus_index) = Penalty%lambda*sign(1.0d0,xMinus)
    P = Penalty%lambda*(sum(xPlus) + sum(xMinus))
    GradP(Penalty%xPlus_index) = Penalty%lambda
    GradP(Penalty%xMinus_index) = Penalty%lambda
  end if

end subroutine PenaltyFunction

subroutine TestGrad(LikeFunc,nx,x,nstate,iuser,ruser)
  use nrtype
  use OutputModule, only : MakeFullFileName
  implicit none
  integer(i4b), intent(in) :: nx,nstate
  real(dp),     intent(in) :: x(:)
  integer(i4b), intent(inout) :: iuser(:)
  real(dp),     intent(inout) :: ruser(:)

  interface
    subroutine LikeFunc(mode,n,x,L,GradL,nstate,iuser,ruser)
      use nrtype
      use GlobalModule, only : Penalty
      integer(i4b), intent(inout) :: mode
      integer(i4b), intent(in)    :: n,nstate
      integer(i4b), intent(inout) :: iuser(*)
      real(dp),     intent(in)    :: x(n)
      real(dp),     intent(inout) :: ruser(*)
      real(dp),     intent(out)   :: L
      real(dp),     intent(inout) :: GradL(n)
    end subroutine LikeFunc
  end interface

  integer(i4b)  :: mode
  real(dp)      :: L0,GradL0(nx)
  real(dp)      :: L1,L2,GradL1(nx),DummyGrad(nx)
  integer(i4b)  :: i1
  real(dp)      :: x1(nx),x2(nx)

  mode=2
  call LikeFunc(mode,nx,x,L0,GradL0,nstate,iuser,ruser)
  
  open(UNIT   = 1033, &
       File   = MakeFullFileName('TestGrad.txt'), &
       action = 'write')
  GradL1 = 0.0d0
  do i1=1,nx
    x1=x
    x2=x
    x1(i1) = x(i1) + 1.0d-6
    x2(i1) = 2.0*x(i1) - x1(i1)
    mode=0
    call LikeFunc(mode,nx,x1,L1,DummyGrad,nstate,iuser,ruser)
    call LikeFunc(mode,nx,x2,L2,DummyGrad,nstate,iuser,ruser)
    GradL1(i1) = (L1-L2)/(x1(i1)-x2(i1))

    write(1033,1033) i1,x(i1),GradL0(i1),GradL1(i1)
    1033 format(i5,3d25.16)
  end do
  close(1033)
end subroutine TestGrad

!------------------------------------------------------------------------------
! subroutine ComputeHess(x,L,GradL,Hess,iuser,ruser)
!------------------------------------------------------------------------------
subroutine ComputeHess(x,L,GradL,Hess,iuser,ruser)
  use nrtype
  implicit none
  real(dp),     intent(in)  :: x(:)
  real(dp),     intent(out) :: L
  real(dp),     intent(out) :: GradL(:)
  real(dp),     intent(out) :: Hess(:,:)
  integer(i4b), intent(inout)  :: iuser(*)
  real(dp),     intent(inout)  :: ruser(*)

  integer(i4b)              :: i1,i2
  integer(i4b)              :: mode,nstate
  real(dp), allocatable     :: x1(:),x2(:)
  real(dp)                  :: L1,L2
  real(dp), allocatable     :: GradL1(:),GradL2(:)
  real(dp)                  :: h
  integer(i4b)              :: n

  n = size(x)
  allocate(x1(n),x2(n),GradL1(n),GradL2(n))
  Hess=0.0d0
  mode=1
  nstate=1
  h=1.0d-6
  call LikeFunc(mode,n,x,L,GradL,nstate,iuser,ruser)
  do i1=1,n
    x1=x
    x2=x
    x1(i1) = x(i1) + h
    x2(i1) = 2.0d0*x(i1) - x1(i1)
    call LikeFunc(mode,n,x1,L1,GradL1,nstate,iuser,ruser)
    call LikeFunc(mode,n,x2,L2,GradL2,nstate,iuser,ruser)
    Hess(:,i1) = (GradL1 - GradL2)/(x1(i1)-x2(i1))
  end do
  Hess = 0.5d0 * (Hess+transpose(Hess))
  deallocate(x1,x2,GradL1,GradL2)
end subroutine ComputeHess

#if OLDHESS==1
! Compute hessian = variance of score of likelihood
subroutine ComputeHess2(x,L,Grad,Hess)
  use GlobalModule, only : HHData,ifree
  use OutputModule, only : WriteHess
  use IFPORT, only : time
  use nrtype
  implicit none
  real(dp), intent(in)  :: x(:)
  real(dp), intent(out) :: L,Grad(:),Hess(:,:)

  integer(i4b)              :: i0,i1,i2,nx
  integer(i4b)              :: nmax,niter
  integer(i4b), allocatable :: blockindex(:)
  real(dp), allocatable     :: LHH0(:),LHH1(:),x1(:)
  real(dp), allocatable     :: GradLHH(:,:),GradLHH1(:)
  real(dp)                  :: h
  integer(i4b)              :: TotalTime
  integer(i4b), allocatable :: SubTime(:)
  character(len=10)         :: StartTime 

  ! when nx >50, break up computation into iterations
  ! compute in max block size of 50 x 50
  nx = size(x,1)
  nmax = iFree%HessNMax
  niter = nx/nmax+1  ! number of iterations
  allocate(blockIndex(nmax))

  allocate(LHH0(HHData%n),LHH1(HHData%n))
  allocate(GradLHH(HHData%n,nmax))
  allocate(GradLHH1(HHData%n))
  allocate(x1(nx))
  allocate(SubTime(nx)) 

  L       = 0.0d0
  Grad    = 0.0d0
  Hess    = 0.0d0
  GradLHH = 0.0d0
  LHH0    = 0.0d0
  LHH1    = 0.0d0
  h       = 1.0e-4

  call date_and_time(time=StartTime)
  print *,'Begin ComputeHess2. Start time is (hhmmss.sss):',StartTime
  if (iFree%flagmue>0) then
    print *,'ifree%xmue:',iFree%xmue(1),maxval(iFree%xmue)
  end if
  if (iFree%flagInvCDiag>0) then
    print *,'ifree%xInvCDiag:',iFree%xInvCDiag(1),maxval(iFree%xInvCDiag)
  end if
  if (iFree%flagInvCOffDiag>0) then
    print *,'ifree%xInvCOffDiag:',iFree%xInvCOffDiag(1),maxval(iFree%xInvCOffDiag)
  end if
  if (iFree%flagBD_beta>0) then
    print *,'ifree%xBD_beta:',iFree%xBD_beta(1),maxval(iFree%xBD_beta)
  end if
  if (iFree%flagBC_beta>0) then
    print *,'ifree%xBC_beta:',iFree%xBC_beta(1),maxval(iFree%xBC_beta)
  end if
  if (iFree%flagBD_CDiag>0) then
    print *,'ifree%xBD_CDiag:',iFree%xBD_CDiag(1),maxval(iFree%xBD_CDiag)
  end if
  if (iFree%flagBD_COffDiag>0) then
    print *,'ifree%xBD_COffDiag:',iFree%xBD_COffDiag(1),maxval(iFree%xBD_COffDiag)
  end if
  if (iFree%flagBC_CDiag>0) then
    print *,'ifree%xBC_CDiag:',iFree%xBC_CDiag(1),maxval(iFree%xBC_CDiag)
  end if
  if (iFree%flagBC_COffDiag>0) then
    print *,'ifree%xBC_COffDiag:',iFree%xBC_COffDiag(1),maxval(iFree%xBC_COffDiag)
  end if

  call Like2Hess(nx,x,LHH0)
  L  = sum(LHH0)/real(HHData%n,dp)
  TotalTime = time()
  SubTime   = 0.0d0

  do i0=iFree%HessIter0,niter
    blockindex = (i0-1)*nmax+(/1:nmax/)
    GradLHH = 0.0d0

    ! fill in diagonal blocks of hessian of size (50,50)
    ! hess(blockindex,blockindex)
    do i1=1,nmax
      if (blockIndex(i1)>nx) then
        exit
      end if
      SubTime(i1) = time()
      x1=x
      x1(blockindex(i1)) = x(blockindex(i1))+h
      call Like2Hess(nx,x1,LHH1)
      GradLHH(:,i1) = (LHH1-LHH0)/(x1(blockIndex(i1))-x(blockIndex(i1)))
      Grad(blockIndex(i1)) = sum(GradLHH(:,i1))/real(HHData%N,dp)
      do i2=1,i1
        Hess(blockIndex(i2),blockIndex(i1)) = &
          sum((GradLHH(:,i1)-Grad(blockIndex(i1))) &
              * (GradLHH(:,i2)-Grad(blockIndex(i2))))/real(HHData%N,dp)
        Hess(blockIndex(i1),blockIndex(i2)) &
           = Hess(BlockIndex(i2),blockIndex(i1))
      end do
      SubTime(i1) = time() - SubTime(i1) ! elapsed time in seconds
      write(6,'(2a4,2a11)') 'nx','i1','time','TotalTime'
      write(6,'(2i4,2i11)') nx,blockindex(i1),SubTime(i1),time()-TotalTime
      Call WriteHess(Hess)
    end do

    ! fill in off-diagonal blocks
    if (blockindex(nmax)<nx) then
      do i1=blockindex(nmax)+1,nx
        SubTime(i1) = time()
        x1=x
        x1(i1) = x(i1)+h
        call Like2Hess(nx,x1,LHH1)
        GradLHH1 = (LHH1-LHH0)/(x1(i1)-x(i1))
        Grad(i1) = sum(GradLHH1)/real(HHData%N,dp)
        do i2=1,nmax
          Hess(blockindex(i2),i1) = &
            sum((GradLHH1-Grad(i1)) &
                * (GradLHH(:,blockindex(i2))-Grad(blockindex(i2))))/real(HHData%N,dp)
          Hess(i1,blockindex(i2)) = Hess(blockindex(i2),i1)
        end do
        SubTime(i1) = time() - SubTime(i1) ! elapsed time in seconds
        write(6,'(2a4,a6,2a11)') 'nx','block','i1','time','TotalTime'
        write(6,'(2i4,i6,2i11)') nx,blockindex(nmax),i1,SubTime(i1),time()-TotalTime
        call WriteHess(hess)
      end do
    end if
  end do

  deallocate(LHH0,LHH1,x1)
  deallocate(GradLHH,GradLHH1)
  deallocate(SubTime)
  deallocate(blockIndex)

end subroutine ComputeHess2
#endif

! Compute hessian = variance of score of likelihood
subroutine ComputeHess2(x,L,Grad,Hess,ComputeHessFlag)
  use GlobalModule, only : HHData,ifree
  use OutputModule, only : WriteHess,ReadWriteGradLHH
#ifndef __GFORTRAN__
  use IFPORT, only : time
#endif
  use nrtype
  implicit none
  real(dp), intent(in)     :: x(:)
  real(dp), intent(out)    :: L,Grad(:),Hess(:,:)
  integer(i4b), intent(in) :: ComputeHessFlag

  integer(i4b)              :: i1,i2,nx,i1min,i1max
  real(dp), allocatable     :: LHH0(:),LHH1(:),x1(:)
  real(dp), allocatable     :: GradLHH(:,:)
  real(dp)                  :: h
  integer(i4b)              :: TotalTime
  integer(i4b), allocatable :: SubTime(:)
  character(len=10)         :: StartTime 

  nx = size(x,1)

  allocate(LHH0(HHData%n),LHH1(HHData%n))
  allocate(GradLHH(HHData%n,nx))
  allocate(x1(nx))
  allocate(SubTime(nx)) 

  L       = 0.0d0
  Grad    = 0.0d0
  Hess    = 0.0d0
  GradLHH = 0.0d0
  LHH0    = 0.0d0
  LHH1    = 0.0d0
  h       = 1.0e-4

  call date_and_time(time=StartTime)
  print *,'Begin ComputeHess2. Start time is (hhmmss.sss):',StartTime
  if (iFree%flagmue>0) then
    print *,'ifree%xmue:',iFree%xmue(1),maxval(iFree%xmue)
  end if
  if (iFree%flagInvCDiag>0) then
    print *,'ifree%xInvCDiag:',iFree%xInvCDiag(1),maxval(iFree%xInvCDiag)
  end if
  if (iFree%flagInvCOffDiag>0) then
    print *,'ifree%xInvCOffDiag:',iFree%xInvCOffDiag(1),maxval(iFree%xInvCOffDiag)
  end if
  if (iFree%flagBD_beta>0) then
    print *,'ifree%xBD_beta:',iFree%xBD_beta(1),maxval(iFree%xBD_beta)
  end if
  if (iFree%flagBC_beta>0) then
    print *,'ifree%xBC_beta:',iFree%xBC_beta(1),maxval(iFree%xBC_beta)
  end if
  if (iFree%flagBD_CDiag>0) then
    print *,'ifree%xBD_CDiag:',iFree%xBD_CDiag(1),maxval(iFree%xBD_CDiag)
  end if
  if (iFree%flagBD_COffDiag>0) then
    print *,'ifree%xBD_COffDiag:',iFree%xBD_COffDiag(1),maxval(iFree%xBD_COffDiag)
  end if
  if (iFree%flagBC_CDiag>0) then
    print *,'ifree%xBC_CDiag:',iFree%xBC_CDiag(1),maxval(iFree%xBC_CDiag)
  end if
  if (iFree%flagBC_COffDiag>0) then
    print *,'ifree%xBC_COffDiag:',iFree%xBC_COffDiag(1),maxval(iFree%xBC_COffDiag)
  end if

  call Like2Hess(nx,x,LHH0)
  L  = sum(LHH0)/real(HHData%n,dp)
  
  if (ComputeHessFlag==1) then
    TotalTime = time()
    SubTime   = 0.0d0

    i1min = iFree%nHess0
    i1max = merge(nx,min(nx,iFree%nHess1),iFree%nHess1==0)
    if (i1min>1) then
      call ReadWriteGradLHH(GradLHH,'read')
    end if
    do i1=i1min,i1max
      SubTime(i1) = time()
      x1=x
      x1(i1) = x(i1)+h
      call Like2Hess(nx,x1,LHH1)
      GradLHH(:,i1) = (LHH1-LHH0)/(x1(i1)-x(i1))
      Grad(i1) = sum(GradLHH(:,i1))/real(HHData%N,dp)
      do i2=1,i1
        Hess(i2,i1) = &
          sum((GradLHH(:,i1)-Grad(i1)) &
              * (GradLHH(:,i2)-Grad(i2)))/real(HHData%N,dp)
        Hess(i1,i2) &
           = Hess(i2,i1)
      end do
      SubTime(i1) = time() - SubTime(i1) ! elapsed time in seconds
      write(6,'(2a4,2a11)') 'nx','i1','time','TotalTime'
      write(6,'(2i4,2i11)') nx,i1,SubTime(i1),time()-TotalTime
      call ReadWriteGradLHH(GradLHH,'write')
      Call WriteHess(Hess)
    end do
  else if (ComputeHessFlag==2) then
    call ReadWriteGradLHH(GradLHH,'read')
    do i1=1,nx
      Grad(i1) = sum(GradLHH(:,i1))/real(HHData%N,dp)
      do i2=1,i1
        Hess(i2,i1) = &
          sum((GradLHH(:,i1)-Grad(i1)) &
              * (GradLHH(:,i2)-Grad(i2)))/real(HHData%N,dp)
        Hess(i1,i2) &
           = Hess(i2,i1)
      end do
    end do
  end if
  Hess = Hess*real(HHData%N,dp)
  Call WriteHess(Hess)
  deallocate(LHH0,LHH1,x1)
  deallocate(GradLHH)
  deallocate(SubTime)

end subroutine ComputeHess2

subroutine SetBounds(x,BL,BU)
  use nrtype
  use GlobalModule, only : iFree,MaxOptions,bayes,parms
  implicit none
  real(dp), intent(in)  :: x(:)
  real(dp), intent(out) :: BL(:)
  real(dp), intent(out) :: BU(:)
  integer(i4b)          :: lastindex,index1(1)
  ! transpose(BC_C) = (KC x JC)
  ! transpose(BD_C) = (KC x JD)  
  integer(i4b) :: KC,KD,JC,JD,i1

  BL = x-10.0d0
  BU = x+10.0d0

  if (iFree%nD>0) then
    BL(iFree%xD) = 1e-6
    BU(iFree%xD) = x(iFree%xD) + 10.0d0
  end if

  if (iFree%nBC>0) then
    BL(iFree%xBC) = 0.1d0 * pi_d
    BU(iFree%xBC) = 0.9d0 * pi_d
    ! BC(:,1:K) column i1 has i1-1 elements in BC
    !           upper bound = 2*pi_d for element i1-1
    !           upper bound = pi_d otherwise
    do i1=2,parms%K
      ! Find elements for which upper bound is 2*pi
      lastindex = (i1*(i1-1))/2 
      if (any(iFree%BC==lastindex)) then
        index1 = pack(iFree%xBC,iFree%BC==lastindex)
        BU(index1) = 1.9d0 * pi_d 
      end if
    end do
    do i1=parms%K+1,parms%J
      lastindex = (parms%K*(parms%K-1))/2 + (i1-parms%K)*(parms%K-1)
      if (any(iFree%BC==lastindex)) then
        index1 = pack(iFree%xBC,iFree%BC==lastindex)
      end if
      BU(index1) = 1.9d0 * pi_d
    end do
  end if

  if (iFree%nMUE>0) then
    BL(iFree%xMUE) = x(iFree%xmue)-10.0d0
    BU(iFree%xmue) = x(iFree%xmue)+10.0d0
  end if

  if (iFree%nMUE_month>0) then
    BL(iFree%xMUE_month) = x(iFree%xmue_month)-10.0d0
    BU(iFree%xmue_month) = x(iFree%xmue_month)+10.0d0
  end if

  if (iFree%nInvCDiag>0) then
    BL(iFree%xInvCDiag) = parms%InvCDiag_LO
    BU(iFree%xInvCDiag) = x(iFree%xInvCDiag)+parms%InvCDiag_HI
  end if

  if (iFree%nInvCOffDiag>0) then
    BL(iFree%xInvCOffDiag) = parms%InvCOffDiag_LO * pi_d
    BU(iFree%xInvCOffDiag) = parms%InvCOffDiag_HI * pi_d
    do i1=2,parms%K
      lastindex = (i1*(i1-1))/2
      if (any(iFree%InvCOffDiag==lastindex)) then
        index1 = pack(iFree%xInvCOffDiag,iFree%InvCOffDiag==lastindex)
        BU(index1) = parms%InvCOffDiag_HI * 2.0d0*pi_d
      end if
    end do 
  end if
  
  if (iFree%nBD_beta>0) then
    BL(iFree%xBD_beta) = max(x(iFree%xBD_beta)-2.0d0,parms%BD_beta_lo)
    BU(iFree%xBD_beta) = min(x(iFree%xBD_beta)+1.0d0,parms%BD_beta_hi)
  end if

  if (iFree%nBD_month>0) then
    BL(iFree%xBD_month) = max(x(iFree%xBD_month)-2.0d0,parms%BD_month_lo)
    BU(iFree%xBD_month) = max(x(iFree%xBD_month)+1.0d0,parms%BD_month_hi)
  end if

  if (iFree%nBC_beta>0) then
    BL(iFree%xBC_beta) = max(x(iFree%xBC_beta)-2.0d0,parms%BC_beta_lo)
    BU(iFree%xBC_beta) = min(x(iFree%xBC_beta)+3.0d0,parms%BC_beta_hi)
  end if
 
  if (iFree%nBC_CDiag>0) then
    ! BC   = norminv( BC_beta * z + BC_C * eta
    !        want diag of BC_C to be between [0,0.8]
    !        harmless to allow for BC_CDiag  in [-1,1]
    !        need to use care when BC_CDiag <=0
    !        also note: Like(BC_CDiag= -x) = Like(BC_CDiag = x)
    BL(iFree%xBC_CDiag) = max(x(iFree%xBC_CDiag)-1.0d0,-1.0d0)
    BU(iFree%xBC_CDiag) = min(x(iFree%xBC_CDiag)+1.0d0,1.0d0)
  end if 
  if (iFree%nBD_CDiag>0) then
    ! log(BD) = BD_Beta * z + BD_C * eta
    !           want diag of BD_C between 0.0  and  1.0
    !        harmless to allow for BD_CDiag  in [-1,1]
    !        need to use care when BD_CDiag <=0
    !        also note: Like(BD_CDiag= -x) = Like(BD_CDiag = x)
    BL(iFree%xBD_CDiag) = max(x(iFree%xBD_CDiag)-1.0d0,-1.0d0)
    BU(iFree%xBD_CDiag) = min(x(iFree%xBD_CDiag)+1.0d0,1.0d0)
  end if 

  if (iFree%nBC_COffDiag>0) then
    BL(iFree%xBC_COffDiag) = 0.10d0*pi_d
    BU(iFree%xBC_COffDiag) = 0.90d0*pi_d
    JC = size(parms%BC_C,1)
    KC = size(parms%BC_C,2)
    do i1=2,KC
      lastindex = (i1*(i1-1))/2
      if (any(iFree%BC_COffDiag==lastindex)) then
        index1 = pack(iFree%xBC_COffDiag,iFree%BC_COffDiag==lastindex)
        BU(index1) = 1.9d0 * pi_d
      end if
    end do
    do i1=KC+1,JC
      lastindex = (KC*(KC-1))/2 +(i1-KC)*(KC-1)
      if (any(iFree%xBC_COffDiag==lastindex)) then
        index1 = pack(iFree%xBC_COffDiag,iFree%BC_COffDiag==lastindex)
        BU(index1) = 1.9d0 * pi_d
      end if
    end do
  end if
  
  if (iFree%nBD_COffDiag>0) then
    BL(iFree%xBD_COffDiag) = 0.10d0*pi_d
    BU(iFree%xBD_COffDiag) = 0.90d0*pi_d
    JD = size(parms%BD_C,1)
    KD = size(parms%BD_C,2)
    do i1=2,KD
      lastindex = (i1*(i1-1))/2
      if (any(iFree%BD_COffDiag==lastindex)) then
        index1 = pack(iFree%xBD_COffDiag,iFree%BD_COffDiag==lastindex)
        BU(index1) = 1.9d0 * pi_d
      end if
    end do
    do i1=KD+1,JD
      lastindex = (KD*(KD-1))/2 +(i1-KD)*(KD-1)
      if (any(iFree%xBD_COffDiag==lastindex)) then
        index1 = pack(iFree%xBD_COffDiag,iFree%BD_COffDiag==lastindex)
        BU(index1) = 1.9d0 * pi_d
      end if
    end do
  end if

  if (MaxOptions%Algorithm==6) then
    bayes%BD_beta_lo = -2.0d0
    bayes%BD_beta_hi = 3.0d0
    bayes%BD_month_lo = -2.0d0
    bayes%BD_month_hi = 3.0d0
    bayes%BD_CDiag_lo = -1.0d0
    bayes%BD_CDiag_hi = 1.0d0
    bayes%BD_COffDiag_lo = 0.10d0*pi
    bayes%BD_COffDiag_hi = 0.90d0 * pi

    bayes%BC_beta_lo = -2.0d0
    bayes%BC_beta_hi = 3.0d0
    bayes%BC_CDiag_lo = -1.0d0
    bayes%BC_CDiag_hi = 1.0d0
    bayes%BC_COffDiag_lo = -.90d0 * pi
    bayes%BC_COffDiag_hi = 0.90d0 * pi
    bayes%MUE_month_lo  = 10.0d0
    bayes%MUE_month_hi  = 30.0d0
    bayes%MUE_lo      = 10.0d0
    bayes%MUE_hi      = 30.0d0
    bayes%InvCDiag_lo = 0.1d0
    bayes%InvCDiag_hi = 1.0d0 / 0.1d0
    bayes%InvCOffDiag_lo = 0.10d0 * pi
    bayes%InvCOffDiag_hi = 0.90d0 * pi

    bayes%nall = 100
    allocate(bayes%x(iFree%nall,bayes%nall))
    allocate(bayes%w(bayes%nall),bayes%prior(bayes%nall))
  end if

end subroutine SetBounds

subroutine PlotLike(LikeFunc,nx,x,xlabels,xlo,xhi,nstate,iuser,ruser)
  use nrtype
  use ToolsModule, only : linspace
  use GlobalModule, only : OutDir
  implicit none
  integer(i4b), intent(in)     :: nx
  real(dp),     intent(in)     :: x(:),xlo(:),xhi(:)
  character(len=*), intent(in) :: xlabels(:)
  integer(i4b), intent(in)     :: nstate
  integer(i4b), intent(inout)  :: iuser(:)
  real(dp),     intent(inout)  :: ruser(:)

  interface
    subroutine LikeFunc(mode,n,x,L,GradL,nstate,iuser,ruser)
      use nrtype
      use GlobalModule, only : Penalty
      integer(i4b), intent(inout) :: mode
      integer(i4b), intent(in)    :: n,nstate
      integer(i4b), intent(inout) :: iuser(*)
      real(dp),     intent(in)    :: x(n)
      real(dp),     intent(inout) :: ruser(*)
      real(dp),     intent(out)   :: L
      real(dp),     intent(inout) :: GradL(n)
    end subroutine LikeFunc
  end interface
  integer(i4b)          :: ix,i1,n,mode
  real(dp), allocatable :: x1(:),xplot(:,:),L(:,:),GradL(:)
  integer(i4b)          :: plot_unit
  character(200)        :: plot_file
  ! 
  plot_unit = 50
  plot_file = trim(OutDir) // '/LikeData.txt' 
  open(UNIT = plot_unit,  &
       FILE = plot_file,  &
       ACTION = 'write')
  mode = 0 
  n=30
  allocate(x1(nx),xplot(n,nx),L(n,nx),GradL(nx))
  x1 = x
  xplot = 0.0d0
  do ix=1,nx
    call linspace(xlo(ix),xhi(ix),n,xplot(:,ix))
    do i1=1,n
    !  xplot(i1,ix) = xlo(ix) + (xhi(ix)-xlo(ix))*real(i1-1,dp)/real(n-1,dp)
      x1(ix) = xplot(i1,ix)
      call LikeFunc(mode,nx,x1,L(i1,ix),GradL,nstate,iuser,ruser)
      write(plot_unit,1762) ix,trim(xlabels(ix)),i1,xplot(i1,ix),L(i1,ix)
1762 format(i4,',',a20,',',i4,',',d25.16,',',d25.16)
      x1(ix) = x(ix)
    end do
  end do
  close(plot_unit)
  deallocate(x1,xplot,L,GradL)
end subroutine PlotLike

subroutine Constraint(x,mode,F,GradF)
  use nrtype
  use GlobalModule, only : parms,iFree,HHData,IntRule
  use nag_library,  only : F04BAF,F07FDF
  use NewTools,     only : ComputeSVD,ComputeInverse_LU,ComputeLQ
  implicit none
  real(dp),     intent(in)    :: x(:)
  integer(i4b), intent(in)    :: mode
  real(dp),     intent(inout) :: F(:)
  real(dp),     intent(inout) :: GradF(:,:)

! Compute constraints:
!
!  for those households that choose nNonZero==K, there are J - nNonZero
!  constraints:
!      p2>= B2'*inv(B1')*p1
!
! for households that choose nNonZero<K, the constraint is that 
!
!    Pr(  p2 >= B2'*(B1*q1 - e) ) > 0
!
! For households that choose nNonZero==K,
!
!    Pr( B'*e <=p) >0
!
!
! F               = (NC x 1) NC = number of constraints
! HHData%N        = (1 x 1)  number of households
! HHData%nNonZero = (N x 1)  number of distinct products purchased
! HHData%iNonZero = (N x 1)  cell array containing indexes of products purchased
!                          HHData%iNonZero{i1} list of products purchased
!                          by household i1
! HHData%p        = (N x J)  prices for each household
! HHData%q        = (N x J)  quantities purchased by household i1
! parms%J         = (1 x 1)  Number of products
! parms%K         = (1 x 1)  Dimension of e
! parms%B         = (K x J)  Matrix in utility function
! parms%GradB_phi = (K x n) gradient of B w.r.t. phi
!                           GradB_phi(:,i1) = gradient of column j(i1) of B
!                           w.r.t. phi(i1)
! parms%GradB_D   = (K x J) gradient of B w.r.t. D
!                           GradB_D(:,i1) = gradient of column i1 of B 
!                           w.r.t. D(i1)
! Revision history
! 2015JUN07 LN  major overhaul to clean up and make like TechnicalAppendix.tex
! 16apr2013 LN  adapted from Constraint.m
! 15nov2012 LN  adapted from Like1.m
real(dp)                  :: small
integer(i4b)              :: n
integer(i4b)              :: JK
real(dp),     allocatable :: c1(:),c2(:)
integer(i4b)              :: n1,n2
integer(i4b)              :: i1,i2

integer(i4b)              :: d1,d2,d3
real(dp),     allocatable :: B1(:,:),B2(:,:)
real(dp),     allocatable :: p1(:),p2(:)

! workspace for F04BAF: NAG matrix inversion using LU decomp
real(dp), allocatable     :: temp1(:)
integer(i4b)              :: ifail
integer(i4b), allocatable :: iPivot(:)
real(dp)                  :: rCond,ErrBound

integer(i4b), allocatable :: index(:)
integer(i4b), allocatable :: index1(:)
integer(i4b), allocatable :: index2(:)
real(dp),     allocatable :: S(:),U(:,:),VT(:,:)
real(dp),     allocatable :: B2Tilde(:,:)
real(dp),     allocatable :: B21Tilde(:,:),B22Tilde(:,:)
real(dp),     allocatable :: M1(:,:),DTilde(:)

real(dp),     allocatable :: C(:,:),Omega(:,:)
real(dp),     allocatable :: Omega11(:,:),Omega22(:,:),Omega12(:,:),C22(:,:)
character(1)              :: UPLO
real(dp),     allocatable :: nu(:),epsilon1(:)
real(dp),     allocatable :: Psi(:,:),CPsi(:,:)
real(dp),     allocatable :: temp2(:,:)
real(dp),     allocatable :: R(:,:),Q(:,:)
integer(i4b)              :: Integrand
integer(i4b), allocatable :: RowGroup(:)
real(dp)                  :: prob

small = 1e-6
call UpdateParms(x,iFree,parms)

! Initial values for likelihood and gradient
n  = count(HHData%nNonZero==parms%K)
JK = parms%J-parms%K;
allocate(c1(n*JK),c2(HHData%N-n))

c1 = 0.0d0
c2 = 0.0d0

! Loop through households
n1=0
n2=0

do i1=1,HHData%N
  d1 = HHData%nNonZero(i1)  ! number of items purchased 
  d2 = parms%K-d1           ! K - d1
  d3 = parms%J-d1           ! number of items NOT purchased
  allocate(B1(parms%K,d1),B2(parms%K,d3))
  allocate(p1(d1),p2(d3))

  B1 = parms%B(:,HHData%iNonZero(1:d1,i1))
  B2 = parms%B(:,HHData%iZero(1:d3,i1))
  p1 = HHData%p(HHData%iNonZero(1:d1,i1),i1)
  p2 = HHData%p(HHData%iZero(1:d3,i1),i1)
   
  ! c(index) = constraints corresponding to household i1
  if (d1==parms%K) then
    allocate(index(JK))
    index = n1 + (/1:JK/)
    n1=n1+JK
    !------------------------------------------------------------------------  
    ! Case 1:   K products were purchased
    !           (HHData%nNonZero == K)
    !           -p2 + B2.'* (B1.')\p1 <=0
    !------------------------------------------------------------------------
    allocate(iPivot(d1))
    allocate(temp1(d1))
    temp1 = p1
    ifail=-1
    B1=transpose(B1)
    call F04BAF(d1,1,B1,d1,iPivot,temp1,d1,rCond,errBound,ifail)
    c1(index) = - p2 + matmul(transpose(B2),temp1)
    deallocate(temp1,iPivot)
    deallocate(index) 
  elseif (d1<parms%K .and. d1>0) then
    allocate(index(1))
    index = n2+1
    n2=n2+1
    !------------------------------------------------------------------------  
    ! Case 2:   d1 < K products were purchased
    !           (HHData%nNonZero == d1 and d1 < K)
    !
    !           -p2 -B2.'*(B1*q1-e)<=0
    !------------------------------------------------------------------------

    allocate(index1(d1),index2(parms%K))
    index1 = (/1:d1/)
    index2 = d1+(/1:parms%K/)

    ! Compute SVD of B1
    ! size(B1) = (K x d1)
    !   B1 = U*S*VT
    allocate(S(d1))
    allocate(U(parms%K,parms%K))
    allocate(VT(d1,d1))
    call ComputeSVD(B1,U,S,VT)

    ! size(B2Tilde) = (K x J-d1)
    allocate(B2Tilde(parms%K,d3))
    allocate(B21Tilde(d1,d3),B22Tilde(d2,d3))
    B2Tilde = matmul(transpose(U),parms%B(:,HHData%iZero(1:d3,i1)))
    B21Tilde = B2Tilde(index1,:)
    B22Tilde = B2Tilde(index2,:)

    ! mean of epsilon
    allocate(nu(parms%K))
    nu = matmul(transpose(U),parms%MuE)

    ! C = inv(InvC)
    ! Sig = C*C.'
    allocate(C(parms%K,parms%K),Omega(parms%K,parms%K))
    call ComputeInverse_LU(parms%K,parms%InvC,C,ifail)
    Omega = matmul(C,transpose(C))
    Omega = matmul(Omega,U)
    Omega = matmul(transpose(U),Omega)

    allocate(Omega11(d1,d1),Omega22(d2,d2),Omega12(d1,d2),C22(d2,d2))

    Omega11 = Omega(index1,index1)
    Omega22 = Omega(index2,index2)
    Omega12 = Omega(index1,index2)

    ! Cholesky decomposition of Omega22: lower triangular form
    !      C22*C22' = Omega22
    C22 = Omega22
    UPLO = 'L'
    call F07FDF(UPLO,d2,C22,d2,ifail)

    allocate(Psi(d1,d1))

    ! Compute  Psi = Omega11 - Omega12*inv(Omega22)*Omega12'
    Psi = 0.0d0
    allocate(iPivot(d2))
    allocate(temp2(d2,d1))
    temp2 = transpose(Omega12)
    call F04BAF(d2,d1,Omega22,d2,iPivot,temp2,d2,rCond,errBound,ifail)
    !  Psi = Omega11 - Omega12*inv(Omega22)*Omega12'
    !  size( inv(Omega22)*Omega12' ) = (d2 x d1)
    Psi = Omega11 - matmul(Omega12,temp2)
    deallocate(temp2,iPivot)
   
    allocate(CPsi(d1,d1))
    ! Cholesky decomposition of omega11: lower triangular form
    !      CPsi * CPsi' = Psi
    CPsi = Psi
    UPLO = 'L'
    call F07FDF(UPLO,d1,CPsi,d1,ifail)
    if (ifail>0) then
      print *, 'Psi is not positive definite'
      print *, 'Psi'
      do i2=1,d1
        print *, Psi(i2,:)
      end do
      print *,Omega11
      do i2=1,d1
        print *,Omega11(i2,:)
      end do
      stop
    end if

    allocate(M1(d3,d2),DTilde(d3))
    ! M1*epsilon2 <= DTilde
    M1 = matmul(transpose(B22Tilde),C22)

    allocate(temp1(d1))
    allocate(epsilon1(d1))
    temp1 = matmul(VT,HHData%p(1:d1,i1))

    ! Compute inv(S1)*VT*p1
    !   epsilon1 = inv(S1)*VT*p1 + S1*VT*q1
    !   DTilde   = p2 - B21Tilde.'*inv(S1)*VT*p1 - B22Tilde.'*S1*VT*q1
    Temp1 = Temp1/S(1:d1)
    
    epsilon1 = Temp1
    DTilde  = HHData%p(HHData%iZero(1:d3,i1),i1)               &
             - matmul(transpose(B21Tilde),Temp1)

    Temp1 = matmul(VT,HHData%q(1:d1,i1))
    Temp1 = Temp1*S(1:d1)
    
    epsilon1 = epsilon1 + Temp1
    DTilde = DTilde - matmul(transpose(B2Tilde),Temp1)
    deallocate(Temp1)

    ! M1 = R*Q  : LQ Decomposition of M1
    !             R = lower triangular
    !             Q = orthogonal
    ! size(M1) = (d3 x d2)
    ! size(R)  = (d3 x d2)
    ! size(Q)  = (d2 x d2)

    allocate(R(d3,d2),Q(d2,d2))
    call ComputeLQ(d3,d2,M1,R,Q,ifail)

    ! Prob = integral of DensityFunc over region of x satisfying R*x<=M2_tilda
    ! 
    Integrand = 2

    allocate(RowGroup(d3))
    call ComputeRowGroup(R,d3,d2,RowGroup)

    call ComputeProb(Prob,                                            &
                     IntRule%rule(d2)%nodes,IntRule%rule(d2)%weights, &
                     RowGroup,R,DTilde,Integrand,                     &
                     epsilon1,nu(index1),Omega12,C22,CPsi,Q(1:d2,:),S(index1),         &
                     d1,d2)

    ! Prob>=small
    ! small - Prob<=0
    c2(index) = small - Prob
   
    deallocate(index)
    deallocate(index1,index2)
    deallocate(U,S,VT)
    deallocate(B2Tilde,B21Tilde,B22Tilde)
    deallocate(M1,DTilde)
    deallocate(C,Omega)
    deallocate(Omega12,Omega11,Omega22,C22)
    deallocate(epsilon1,nu)
    deallocate(Psi,CPsi)
    deallocate(iPivot,temp2)
    deallocate(R,Q)
    deallocate(RowGroup)
    !------------------------------------------------------------------------  
    ! Case 3:   d == 0 products were purchased
    !           (HHData%nNonZero == d and d == 0)
    !           Mapping from 0 to e is NOT one-to-one
    !           need to integrate across region of e-space
    !           satisfying B'*e<=p
    !------------------------------------------------------------------------
  elseif (d1==0) then
    allocate(index(1))
    index = n2+1
    n2    = n2+1

    ! need inverse of InvC
    allocate(C(parms%K,parms%K))
    call ComputeInverse_LU(parms%K,parms%InvC,C,ifail)

    allocate(M1(d3,d2),DTilde(d3))
    M1 = matmul(transpose(parms%B),C)

    allocate(R(d3,d2),Q(d3,d2))
    ! M1 = R*Q
    !    R = (d3 x d2) lower triangular
    call ComputeLQ(d3,d2,M1,R,Q,ifail)
    DTilde = p2 - matmul(transpose(parms%B),parms%MuE)
      
    ! Prob = integral of DensityFunc over region of x satisfying R*x<=M2_tilda
    ! 
    Integrand = 3
    
    allocate(RowGroup(d3))
    call ComputeRowGroup(R,d3,d2,RowGroup)
   
    call ComputeProb(Prob,IntRule%rule(d2)%nodes,IntRule%rule(d2)%weights, &
                     RowGroup,R,DTilde,Integrand)
    c2(index) = small - Prob
    deallocate(index)
    deallocate(C)
    deallocate(M1,DTilde)
    deallocate(R ,Q)
    deallocate(RowGroup)
  end if       ! if HHData%nNonZero(i1)==parms%k
  
   deallocate(B1,B2)
   deallocate(p1,p2)
end do         ! do i1=1,HHData%N

F(1:n*JK) = c1
F(1+n*JK:n*JK+HHData%N-n)= c2

deallocate(c1,c2)
end subroutine Constraint

subroutine NAGConstraintWrapper(mode,nc_nonlin,N,LDCJ,NEEDC,X,CCON,CJAC,nstate,iuser,ruser)
  use nrtype
  implicit none
  integer(i4b), intent(inout) :: mode
  integer(i4b), intent(in)    :: nc_nonlin,n,ldcj
  integer(i4b), intent(in)    :: needc(nc_nonlin)
  real(dp),     intent(in)    :: x(n)
  real(dp),     intent(out)   :: ccon(max(1,nc_nonlin))
  real(dp),     intent(inout) :: cjac(ldcj,n)
  integer(i4b), intent(in)    :: nstate
  integer(i4b), intent(inout) :: iuser(*)
  real(dp),     intent(inout) :: ruser(*)
  real(dp), allocatable       :: F(:),GradF(:,:)

  allocate(F(nc_nonlin),GradF(ldcj,n))
  call Constraint(x,mode,F,GradF)
  if (mode==0) then
    ! compute ccon(needc)
    ccon(needc) = F(needc)
  else if (mode==1) then
    ! compute cjac(:,needc)
    cjac(needc,:) = GradF(needc,:)
  else if (mode==2) then
    ! compute ccon(needc) and cjac(:,needc)
    ccon(needc)   = F(needc)
    cjac(needc,:) = GradF(needc,:)
  end if
  deallocate(F,GradF)
end subroutine NAGConstraintWrapper

! Compute number of nonlinear constraints
subroutine ComputeNC(NC)
  use nrtype
  use GlobalModule, only : HHData,parms
  implicit none
  integer(i4b), intent(out) :: NC
  integer(i4b)              :: n,JK
  n = count(HHData%nNonZero==parms%K)
  JK = parms%J-parms%K
  NC = n*JK + HHData%N-n
end subroutine ComputeNC


subroutine Max_E04JCF(x,LValue,Grad,Hess,ierr)
  use nrtype
#if USE_MPI==1
  use mpi
#endif
  USE GLOBALMODULE, ONLY : PARMS,MASTERID,HHDATA,IFREE,InputDir,OutDir
  use nag_library, only : E04JCF,E04JCP
  real(dp),     intent(inout) :: x(:)
  real(dp),     intent(out)   :: LValue,Grad(:),Hess(:,:)
  integer(i4b), intent(out)   :: ierr

  INTEGER(i4b) 	              :: nx, NPT, MAXCAL, NF, IUSER(2),i1
  real(dp)                    :: RHOBEG, RHOEND, RUSER(1)
  real(dp), allocatable       :: BL(:),BU(:)

  integer(i4b)                :: npt_lo,npt_hi
  integer(i4b)                :: E04UNIT
  character(len=128)          :: E04FILE
  character(len=50)           :: TempString
  character(len=100)          :: ResultFile

  ! variables to compute gradient of OBJFUN
  real(dp) :: h,L1,L2
  real(dp), allocatable :: x0(:),x1(:),x2(:)
  integer(i4b)          :: inform
  integer(i4b)          :: ntest

  nx = size(x,1)

  allocate(BL(NX),BU(nx))
  allocate(x0(nx))
  x0 = x
  BL = -5.0d0
  BU = 5.0d0
  call SetBounds(x,BL,BU)

  iuser(1) = parms%model
  iuser(2) = HHData%N
  ruser    = 1.0d0   ! set = 1.0 for first pass on function evaluation

#if USE_MPI==1
  ! broadcast (nx,iuser,ruser)
  call mpi_bcast(nx,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  ! broadcast length of iuser
  call mpi_bcast(2,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iuser,2,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  ! broadcast length of ruser
  call mpi_bcast(1,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ruser,1,MPI_DOUBLE_PRECISION,MasterID,MPI_COMM_WORLD,ierr)
#endif

  LVALUE = 0.0d0
  Grad = 0.0d0
  Hess = 0.0d0
  do i1=1,nx
    Hess(i1,i1) = 1.0d0
  end do

  ! set options for E04JCF
  E04Unit = 50
  E04File = trim(InputDir) // '/e04jcf_options.txt'
  open(UNIT = E04Unit,  &
       FILE = E04File, &
       ACTION = 'read')
  ! TempString = maxcal,rhobeg,rhoend,npt
  read(50,*) TempString
  read(50,'(i6)') maxcal      ! maximum nuumber of function calls
  read(50,'(d12.5)') rhobeg   ! initial trust region radius
  read(50,'(d12.5)') rhoend   ! final trust region radius
  read(50,'(i4)')    npt      ! number of interpolation points

  ! number of interpolation conditions
  npt_lo = 2*nx+1
  npt_hi = (nx+1)*(nx+2)/2
  if (npt<npt_lo) then
    npt = npt_lo
  end if
  if (npt>npt_hi) then
    npt = npt_hi
  end if
  close(E04Unit)

  !npt = 2*nx+1
  !rhobeg = 0.1d0  ! initial lower boubnd on the value of trust region
  !rhoend = 1.0e-4  ! final lower bound on value of trust region
  !maxcal = 5000    ! maximum number of calls to objective function
 
  ierr = -1

  ! initialize last_x for use by monfun_E04JCF
  allocate(last_x(nx))
  last_x = x

  print *,'Begin maximization using E04JCF'

  !call  E04JCF(OBJFUN_E04JCF,NX,NPT,X,BL,BU,RHOBEG,RHOEND,E04JCP,MAXCAL,LVALUE,NF,IUSER,RUSER,ierr)
  call  E04JCF(OBJFUN_E04JCF,NX,NPT,X,BL,BU,RHOBEG,RHOEND,monfun_E04JCF,MAXCAL,LVALUE,NF,IUSER,RUSER,ierr)

  print *,'E04JCF complete.'
  write(6,98) 'NF','RHOEND','LVALUE','ierr'
  write(6,99) NF,RHOEND,LVALUE,ierr
98 format(a4,2a12,a4)
99 format(i4,2d12.5,i4)


! Compute gradient
allocate(x1(nx),x2(nx))
h = 1.e-4
do i1=1,nx
  x1 = x
  x1(i1) = x(i1) + h
  x2 = x
  x2(i1) = 2.0d0*x(i1) - x1(i1)
  call OBJFUN_E04JCF(nx,x1,L1 ,iuser,ruser,inform)
  call OBJFUN_E04JCF(nx,x2,L2 ,iuser,ruser,inform)
  Grad(i1) = (L1-L2)/(x1(i1)-x2(i1))
end do

deallocate(x1,x2)

!  call ComputeHess(x,LValue,Grad,Hess,iuser,ruser)


  ResultFile = trim(OutDir) // '/result.txt'
  open(unit = 130,file = ResultFile,action = 'write')
  write(130,'(a20,3a25)') 'Var. Name','x0','x','Gradient'
  do i1=1,nx
    write(130,'(a20,3d25.12)') trim(iFree%xlabels(i1)), x0(i1),x(i1), Grad(i1)
  end do
    write(130,'(a25,d25.12)') 'LValue',LValue
    write(130,'(a25,i4)')     'ierr',ierr

#ifdef PLOT_LIKE
    ntest=10
    allocate(x1(nx))
    do i1=1,nx
    do i2=1,ntest
      x1=x
      x1(i1) = x(i1) + 10.0d0 * real(i2,dp)/real(ntest,dp) - 5.0d0
      call OBJFUN_E04JCF(nx,x1,L1 ,iuser,ruser,inform)
      write(130,'(a20,i4,2d25.12)') trim(iFree%xlabels(i1)),i2,x1(i1),L1
    end do
    end do
    deallocate(x1)
#endif

close(130)

deallocate(last_x)

#if USE_MPI==1
! broadcast signal indicating that all mpi tasks are complete
task=0
call mpi_bcast(task,1,MPI_INTEGER,MasterID,MPI_COMM_WORLD,ierr)
#endif


  deallocate(BL,BU)
end subroutine Max_E04JCF

subroutine OBJFUN_E04JCF(n,x,L ,iuser,ruser,inform)
  use nrtype
  implicit NONE
  integer(i4b), intent(in)    :: n
  real(dp),     intent(in)    :: x(n)
  real(dp),     intent(out)   :: L
  integer(i4b), intent(inout) :: iuser(*)
  real(dp),     intent(inout) :: ruser(*)
  integer(i4b), intent(out)   :: inform

  integer(i4b) :: mode,nstate
  real(dp), allocatable :: gradL(:)

  allocate(GradL(n))
  GradL = 0.0d0
  mode  = 0   ! level of likelihood only
  nstate =  0
  call LikeFunc(mode,n,x,L,GradL,nstate,iuser,ruser)
  inform = 0
  deallocate(gradL)
end subroutine OBJFUN_E04JCF

subroutine monfun_E04JCF(n,nf,x,f,rho,iuser,ruser,inform)
  use nrtype
  use OutputModule, only : MakeFullFilename
  implicit none

  !       .. Scalar Arguments ..
  Real(dp),     Intent (In)      :: f, rho
  Integer(i4b), Intent (Out)     :: inform
  Integer(i4b), Intent (In)      :: n, nf

  !       .. Array Arguments ..
  Real (dp),    Intent (Inout)   :: ruser(*)
  Real (dp),    Intent (In)      :: x(n)
  Integer(i4b), Intent (Inout)   :: iuser(*)

  real(dp)                       :: delta_x

  ! unit number for output file
  integer(i4b)                     :: nout
  character(len=20)                :: filename

  !       .. Executable Statements ..
  inform = 0
  nout = 6   ! output = standard io
  filename = MakeFullFilename('e04jcf.log')
  if (nout .ne. 6) then
    open(unit=nout, &
         file=filename, &
         action='write')
  end if

  if (ruser(1)==1.0d0) then
    write(nout,100) 'NF','rho','F','Delta X'
    100 format(a5,3a12)
    ruser(1) = 0.0d0
  end if

  delta_x = sqrt(sum((x-last_x)*(x-last_x)))
  last_x  = x

  write(nout,101) NF,rho,F,delta_X
  101 format(i5,3d12.5)

  if (nout .ne. 6) Then
    close(nout)
  end if
end Subroutine monfun_E04JCF

end module LikelihoodModule
