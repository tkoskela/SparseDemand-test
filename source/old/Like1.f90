subroutine Like1(mode,nx,x,L,GradL,nstate,iuser,ruser)
  use nrtype
  use GlobalModule, only : parms,HHData,QuadRule,iFree
  implicit none

  integer(i4b), intent(inout)        :: mode
  integer(i4b), intent(in)           :: nx
  real(dp),     intent(in)           :: x(:)
  real(dp),     intent(out)          :: L
  real(dp),     intent(inout)        :: GradL(:)
  integer(i4b), intent(in)           :: nstate
  integer(i4b), intent(in)           :: iuser(*)
  real(dp),     intent(in)           :: ruser(*)

!
! Compute likelihood and its gradient w.r.t. x where utility is
!
!     u  = y - p'*q - 0.5 * (B*q - e)' * (B*q - e)
!
! 
! HHData%N        = (1 x 1)  number of households
! HHData%nNonZero = (N x 1)  number of distinct products purchased
! HHData%iNonZero = (N x 1)  cell array containing indexes of products purchased
!                            HHData%iNonZero{i1} list of products purchased
!                            by household i1
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
! QuadRule.v      = (dMax x 1)  cell array of quadrature rules
!                               QuadRule.v{d} = (n(d) x d) d-dimensional quad rule
! QuadRule.w      = (dMax x 1)  cell array of quadrature rule weights
!                               QuadRule.w{d} = (n(d) x 1) weights for
!                               d-dimensional rule
! Revision history
! 07DEC2012 LN  adapted from Like1.m
!
integer(i4b)          :: TestFlag
real(dp)              :: small
integer(i4b)          :: i1
integer(i4b)          :: d1,d2
real(dp), allocatable :: B1(:,:)
real(dp), allocatable :: e(:)

real(dp) :: F
real(dp) :: GradF_e(parms%K)
real(dp) :: GradF_MuE(parms%K)
real(dp) :: GradF_InvC(parms%K,parms%K)


! workspace for NAG matrix inversion and determinants
integer(i4b), allocatable :: iPivot(:)
real(dp)                  :: rCond,errBound
integer(i4b)              :: iFail
real(dp)                  :: D
integer(i4b)              :: ID

external F04BAF  ! matrix inversion using LU decomposition
external F07ADF  ! compute LU factorization of matrix
external F03BAF  ! compute determinant of matrix after F07ADF

call UpdateParms(x,iFree,parms)

! Initial values for likelihood and gradient
L     = 0.0d0
GradL = 0.0d0

small=1e-6   ! avoid log of zero
TestFlag=0   ! 0 no test
             ! 1 shut down inv(B1)
             ! 2 shut down B1 but turn on inv(B1)

! Loop through households
do i1=1,HHdata%N
  d1 = HHData%nNonZero(i1)
  d2 = parms%K-d1
  !------------------------------------------------------------------------  
  ! Case 1:   K products were purchased
  !           (HHData%nNonZero == K)
  !           Mapping from q(iNonZero) to e is one-to-one
  !------------------------------------------------------------------------
  if (d1==parms%K) then
    allocate(B1(d1,d1))
    allocate(e(d1))
    allocate(iPivot(d1))
    B1 = parms%B(:,HHData%iNonZero(1:d1,i1))
    if (TestFlag==0) then
      ! Compute  e = inv(B1')*p1 + B1*q1
      e = HHData%p(i1,HHData%iNonZero(1:d1,i1))
      call F04BAF(d1,1,transpose(B1),d1,iPivot,e,d1,rCond,errBound,ifail)
      e = e + matmul(B1,transpose(HHData%q(i1,HHData%iNonZero(1:d1,i1)))
    else if (TestFlag==1) then
      ! e = B1*q1
      e = matmul(B1,transpose(HHData%q(i1,HHData%iNonZero(1:d1,i1))))
    else if (TestFlag==2) then
      ! e = inv(B1')*p1
      e = HHData%p(i1,HHData%iNonZero(1:d1,i1))
      call F04BAF(d1,1,transpose(B1),d1,iPivot,e,d1,rCond,errBound,ifail)
    end if

    call LogDensityFunc(e,parms,F,GradF_e,GradF_MuE,GradF_InvC)

    ! F07ADF:  compute LU factorization of B1
    call F07ADF(d1,d1,B1,d1,iPivot,ifail)
    ! F03BAF:  compute determinant of B1 after facorizing
    call F03BAF(d1,B1,d1,iPivot,D,ID,ifail)
    ! log(det(B1)) = log(D) + ID*log(2.0)
    L = L + F + log(D) + real(ID)*log(2.0d0)
     
    ! Compute gradient w.r.t. elements of D
    integer(i4b) :: i2,i3
    integer(i4b) :: iXFree
    integer(i4b) :: j1,j2
    real(dp)     :: GradB1(parms%K,parms%K)
    real(dp)     :: GradInvB1T(parms%K,parms%K)
    real(dp)     :: InvB1T(parms%K,parms%K)

    do i2=1,length(iFree%D)
      iXFree = iFree%D(i2)
      j1     = iFree%DX(i2)
      if any(j1==HHData%iNonZero(1:d1,i1)) then
        ! j2 selects column in B1 corresponding to j1
        j2 = pack(HHData%iNonZero(1:d1,i1),(j1==HHData%iNonZero(1:d1,i1))
        GradB1 = 0.0d0
        GradB1(:,j2) = parms%GradB_D(:,j1)
    
        ! gradient of inverse of transpose of B1
        ! compute InvB1T = inv(B1')
        InvB1T = 0.d0
        do i3=1,d1
          InvB1T(i3,i3) = 1.0d0
        end do
        call F04BAF(d1,d1,transpose(B1),d1,iPivot,InvB1T,d1,rCond,errBound,ifail)
        GradInvB1T = -matmul(InvB1T,matmul(GradB1,InvB1T))
        GradLogDet = inv(B1)
        if (TestFlag==0) then
          GradL(iXFree) = GradL(iXFree)    &
            + matmul(transpose(GradF_e),   &
                     matmul(GradInvB1T,    &
                            transpose(HHData%p(i1,HHData%iNonZero(1:d1,i1)))) &
                    +matmul(GradB1,        &
                            transpose(HHData%q(i1,HHData%iNonZero(1:d1,i1))))) &
            + matmul(GradLogDet(j2,:),tranpose(GradB1(j2,:)))
        else if (TestFlag==1) then
          GradL(iXFree) = GradL(iXFree) ...
            + matmul(transpose(GradF_e), &
                     matmul(GradB1,      &
                            tranpose(HHData%q(i1,HHData%iNonZero(1:d1,i1)))) &
            + matmul(GradLogDet(j2,:),tranpose(GradB1(j2,:)))
        elseif (TestFlag==2) then
          GradL(iXFree) = GradL(iXFree) ...
            + matmul(transpose(GradF_e), &
                     matmul(GradInvB1T,  &
                            tranpose(HHData%p(i1,HHData%iNonZero(1:d1,i1))))) &
            + matmul(GradLogDet(j2,:),transpose(GradB1(j2,:)))
        end if ! if TestFlag==0
      end if   ! if any(j1==HHData%iNonZero(1:d1,i1))
    end if     ! for i2= 1:length(iFree.D)
      
    ! Compute gradient w.r.t. phi
    do i2=1,length(iFree%phi)
      iXFree = iFree%phi(i2)
      ! j1 = column in B corresponding to current element of phi   
      j1 = FindColumn(iFree.phi{2}(i2),parms%K);
      if any(j1==HHData%iNonZero{i1})
        ! j2 selects column in B1 corresponding to j1
        j2 =(j1==HHData%iNonZero{i1});
        GradB1 = zeros(parms%K,parms%K);
        GradB1(:,j2) = parms%GradB_phi(:,iFree.phi{2}(i2));
        ! gradient of inverse of transpose of B1
        GradInvB1T = -inv(B1.')*GradB1.'/B1.';
        GradLogDet = inv(B1);
        if TestFlag==0
          GradL(iXFree) = GradL(iXFree) ...
            + (GradF.e.'*(GradInvB1T*HHData%price(i1,HHData%iNonZero{i1}).' ...
                        + GradB1*HHData%q(i1,HHData%iNonZero{i1}).'))     ...
            + GradLogDet(j2,:)*GradB1(j2,:).';
        elseif TestFlag==1
          GradL(iXFree) = GradL(iXFree) ...
            + (GradF.e.'*(GradB1*HHData%q(i1,HHData%iNonZero{i1}).')) ...
            + GradLogDet(j2,:)*GradB1(j2,:).';
        elseif TestFlag==2
          GradL(iXFree) = GradL(iXFree) ...
            + (GradF.e.'*(GradInvB1T*HHData%price(i1,HHData%iNonZero{i1}).')) ...
            + GradLogDet(j2,:)*GradB1(j2,:).';
        end ! if TestFlag==0
      end   ! if any(j1==HHData%iNonZero{i1})
    end     ! for i2=1:length(iFree.phi{1})
    
    
    ! Gradient w.r.t. parms%MuE
    GradL(iFree.MuE{1}) = GradL(iFree.MuE{1})+GradF.MuE(iFree.MuE{2});
    
    ! Gradient w.r.t. parms%InvCSig_D
    for i2=1:length(iFree.InvCSig_D{1})
      iXFree = iFree.InvCSig_D{1}(i2);
      j1 = iFree.InvCSig_D{2}(i2);
      GradL(iXFree) = GradL(iXFree) ...
          + GradF.InvCSig(:,j1).'*parms%GradInvCSig_D(:,j1);
    end
    
    ! Gradient w.r.t. parms%InvCSig_phi
    for i2=1:length(iFree.InvCSig_phi{1})
      iXFree = iFree.InvCSig_phi{1}(i2);
      j1 = iFree.InvCSig_phi{2}(i2);
      GradL(iXFree) = GradL(iXFree) ...
          + GradF.InvCSig(:,j1).'*parms%GradInvCSig_phi(:,j1);
    end
    
  !------------------------------------------------------------------------  
  ! Case 2:   d < K products were purchased
  !           (HHData%nNonZero == d and d < K)
  !           Mapping from q(iNonZero) to e is NOT one-to-one
  !           need to integrate across region of e-space
  !           satisfying M1*e2<=M2
  !------------------------------------------------------------------------
  elseif d1<parms%K && d1>0
      
    index1 = (1:d1)';
    index2 = (d1+1:parms%K)';
    ! size(B1) = (K x d1)
    B1 = parms%B(:,HHData%iNonZero{i1});
    [U,S,V]=svd(B1);
    S1 = S(index1,index1);
    
    ! Define eTilda = U'*e
    ! 
    ! FOC can be solved for eTilda1 as function of q and p:
    !     eTilda1 = G0*p1 + G1*q1
    ! G0' = V*inv(S1)
    ! G1' = B1'*B1*V*inv(S1)
    !
   
    
    G0 = V/S1;
    G1 = (B1.'*B1)*G0;
    G0 = (HHData%price(i1,HHData%iNonZero{i1})*G0).';

    ! size(B2) = (K x J-d1)
    B2Tilda = U.'*parms%B(:,HHData%iZero{i1});
    B21 = B2Tilda(index1,:);
    B22 = B2Tilda(index2,:);
    
    ! M1*e2 <= M2
    M1  = B22.';
    M2  = HHData%price(i1,HHData%iZero{i1}).' ...
          - B21.'*G0 ...
          + B21.'*(S1*V.'-G1) ...
             *HHData%q(i1,HHData%iNonZero{i1}).';
    ! CSig = inv(InvCSig)
    ! Sig = CSig*CSig.';
    CSig = parms%InvCSig\eye(parms%K);
    SigTilda = U.'*(CSig*CSig.')*U;
    
    sig11 = SigTilda(index1,index1);
    sig22 = SigTilda(index2,index2);
    sig12 = SigTilda(index1,index2);
    C22   = chol(sig22,'upper');
    
    ! M1_tilda*z2 <= M2_tilda
    ! e2 = C2*z2+MuE(i2)
    MuTilda = U.'*parms%MuE;
    M1_tilda = M1*C22.';
    M2_tilda = M2 - M1*MuTilda(index2);
    [Q,R]=qr(M1_tilda.');
    omega11 = sig11 - sig12*(sig22\sig12.');
    [cc,pp]=chol(omega11);
    if pp>0
      keyboard
    end
    ! Prob = integral of DensityFunc over region of x satisfying R*x<=M2_tilda
    ! 
    Integrand = @(x)DensityFunc(x,G0,G1,Q.',               ...
                                MuTilda,sig12,C22,omega11, ...
                                d1,index1,                 ...
                                HHData%q(i1,HHData%iNonZero{i1}).');
    
    RowGroup = ComputeRowGroup(R.');

    Prob     = ComputeProb(QuadRule.x{d2},QuadRule.w{d2},RowGroup, ...
                           R.',M2_tilda,Integrand);

    L = L + log(Prob+small);
    
  !------------------------------------------------------------------------  
  ! Case 3:   d == 0 products were purchased
  !           (HHData%nNonZero == d and d == 0)
  !           Mapping from 0 to e is NOT one-to-one
  !           need to integrate across region of e-space
  !           satisfying B'*e<=p
  !------------------------------------------------------------------------
  elseif d1==0

    C = parms%InvCSig\eye(parms%K);
    M1 = parms%B.'*C;  
    [Q,R]=qr(M1.');
    M2 = HHData%price(i1,:)' - parms%B.'*parms%MuE;
      
    ! Prob = integral of DensityFunc over region of x satisfying R*x<=M2_tilda
    ! 
    Integrand = @(x)DensityFunc0(x);
    
    RowGroup = ComputeRowGroup(R.');

    if any(~isreal(R(:)))
      keyboard
    end
    Prob     = ComputeProb(QuadRule.x{d2},QuadRule.w{d2},RowGroup, ...
                           R.',M2,Integrand);

    L = L + log(Prob+small);

  end       ! if HHData%nNonZero(i1)==parms%k
  !------------------------------------------------------------------------

  deallocate(B1,e,iPivot)
  
end           ! for i1=1:HHData%n

L = -L;
GradL = -GradL;

