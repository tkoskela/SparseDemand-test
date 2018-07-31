! Revision history
! 2015AUG14 LN  add MatrixToSphere and MapToSpherical.
!               Update MatrixInverse to Nag Mark 25 libraries
!               edit SphereToMatrix to work on rectangular matrix
module NewTools
! Line No. | Procedure name
!----------|-------------------------------
! 29       | subroutine MatrixInverse(M,InvM,MatrixType)
! 138      | subroutine SphereToMatrix(phi,r,K,J,B,GradB_phi,GradB_r) 
! 215      | subroutine MapToCartesian(r,phi,x,dx)
! 271      | subroutine FindColumn(i1,K,j1)
! 292      | subroutine ComputeInverse_LU(K,A,InvA,ifail)
! 318      ! subroutine ComputeLQ(m,n,A,L,Q,ifail)
! 351      ! function det(A)
! 399      ! subroutine ComputeSVD(B,S,U,VT)

contains

subroutine ComputeMatrixType(M,MatrixType)
  use nrtype
  implicit none
  real(dp),         intent(in)  :: M(:,:)
  character(len=*), intent(out) :: MatrixType

  integer(i4b)                  :: N,i1,ix
  logical                       :: upper,lower
  ! Matrix type
  !  1) Lower triangular
  !  2) Upper triangular
  !  3) Other not symmetric
  !  4) Symmetric
  n = size(M,1)

  if (all((M - transpose(M))==0.0d0)) then
    MatrixType = 'Symmetric'
  else
    upper = .true.
    lower = .true.
    do i1=1,N
      upper = merge(upper,.false.,all(M(i1,(/(ix,ix=i1+1,N)/))==0.0d0))
      lower = merge(lower,.false.,all(M(i1+1:N,i1)==0.0d0))
    end do
    if (upper) then
      MatrixType = 'Upper triangular'
    elseif (lower) then
      MatrixType = 'Lower triangular'
    else
      MatrixType = 'Other not symmetric'
    end if
  end if
end subroutine ComputeMatrixType
!------------------------------------------------------------------------------
!
! subroutine MatrixInverse(M,InvM,MatrixType)
!    Compute InvM = inv(M)
!    MatrixType is one of 4 types:
!      'Lower triangular'
!      'Upper triangular'
!      'Other not symmetric'
!      'Symmetric'
!------------------------------------------------------------------------------
subroutine MatrixInverse(M,InvM,RawMatrixType)
  use nrtype
  use nag_library, only: X02AJF,F01ABF,F01BLF,F04AEF,F07AAF,F07MDF,F07MJF,F07TJF
  ! X02AJF :  compute machine epsilon
  ! F01ABF :  inverse of symmetric positive definite matrix, iterative
  ! F01BLF :  pseudo inverse of matrix
  ! F04AEF :  solve A*X = B
  ! F07AAF :  inverse of matrix
  ! F07MDF :  compute Bunch-Kaufman  factorization of real symmetric indefinite matrix
  ! F07MJF :  invert matrix after factoring
  ! F07TJF :  invert triangular matrix
  implicit none
  real(dp),                   intent(in)  :: M(:,:)
  real(dp),                   intent(out) :: InvM(:,:)
  character(len=*), optional, intent(in)  :: RawMatrixType
  character(len=50)                       :: MatrixType
  !variables used by F01ABF to compute inverse 
  !   F01ABF:  invert pos.def. symmetric matrix
  real(dp), allocatable :: A(:,:),A1(:,:),B(:,:),Z(:)
  integer(i4b)          :: n,IA,IB,ifail

  ! variables used by F07MDf and F07MJF to compute inverse
  !  F07MDF  invert symmetric indefinite matrix
  character(LEN=1)          :: UPLO
  integer(i4b)              :: LWORK,info
  integer(i4b), allocatable :: IPIV(:)
  real(dp),     allocatable :: WORK(:)

  ! variables used by F01BLF
  !   pseudo-inverse of matrix
  integer(i4b)               :: IRANK
  integer(i4b), allocatable  :: INC(:)
  real(dp)                   :: tol
  real(dp),     allocatable  :: AIJMAX(:),D(:),U(:,:),DU(:)

  ! variables used by F07AAF :matrix inversion
  integer(i4b), allocatable  :: pivot(:)

  integer(i4b)  :: i1,i2

  if (present(RawMatrixType)) then
    MatrixType = RawMatrixType
  else
    call ComputeMatrixType(M,MatrixType)
  end if
  ! Matrix Types:
  !     Lower triangular
  !     Upper triangular
  !     Other not symmetric
  !     Symmetric positive definite

  n = size(M,1)

  select case (MatrixType)

  case('Lower triangular')
    InvM = M
    call F07TJF('L','N',N,InvM,N,INFO)
  case('Upper triangular')
    InvM = M
    call F07TJF('U','N',N,InvM,N,INFO)
  case('Other not symmetric')
    allocate(B(N,N),Z(N),U(N,N),A1(N,N))
    B = 0.0d0
    do i1=1,N
      B(i1,i1) = 1.0d0
    end do
    ! B = identity
    ! z = workspace
    ! U = LU decomp of M
    ! A1 = B - A*InvM
    ! Computes InvM = inv(M)*B
    !
    call F04AEF(M,N,B,N,N,N,InvM,N,Z,U,N,A1,N,ifail)
    deallocate(B,Z,U,A1)
  case('Symmetric')
    ! Try to use F01ABF.
    IA = n+1
    IB = n
    allocate(A(N,N),A1(IA,N),B(N,N),Z(N))
    A = 0.0d0
    A1 = 0.0d0
    B = 0.0d0
    Z = 0.0d0
    A1(1:n,:) = M
    IFAIL=1
    call F01ABF(A1,IA,n,B,IB,Z,IFAIL)

    if (IFAIL==0) then
      InvM = B
    elseif (IFAIL==2) then
      A=M
      B=M
      tol=0.0d0
      do i1=1,n
      do i2=1,n
        tol=tol+B(i1,i2)*B(i1,i2)
      end do
      end do
      tol = dsqrt(tol)*X02AJF()

      ! Compute pseudo-inverse
      allocate(INC(N),AIJMAX(N),D(N),U(N,N),DU(N))
      call F01BLF(N,N,tol,B,N,AIJMAX,IRANK,INC,D,U,N,DU,IFAIL)
      deallocate(INC,AIJMAX,D,U,DU)
      B = transpose(B)
    
      A = matmul(B,A)
      ! Compute solution x to (B*A)*x = B
      allocate(pivot(n))
      call F07AAF(n,n,A,n,pivot,B,n,ifail)
      deallocate(pivot)
      if (ifail .ne. 0) then
        print *,'Both F01ABF and F07AAF failed to converge. Matrix is ill conditioned.'
        InvM = 0.0d0
        do i1=1,n
          InvM(i1,i1) =  1.0d0
        end do
      end if
    elseif (IFAIL==1) then
      ! M is not Positive definite, use F07 routines to invert
      B=M
      UPLO = 'L'
      if (n>10000) then
        LWORK=n
      elseif (n<=10000) then
        LWORK = min(10000,n*64)
      end if
      allocate(IPIV(n))
      allocate(WORK(LWORK))
      info=0

      ! compute Bunch-Kaufman factorization of M
      call F07MDF(UPLO,n,B,n,IPIV,WORK,LWORK,info)

      if (info<0) then
        print *, 'Invalid parameter values in F07MDF. info =',info
        ifail = info
        InvM = 0.0d0
        do i1=1,n
          InvM(i1,i1) = 1.0d0
        end do
      elseif (info>0) then
        print *, 'error in F07MDF. block diagonal matrix D is singular.'
        ifail = info
        InvM = 0.0d0
        do i1=1,n
          InvM(i1,i1) = 1.0d0
        end do
      elseif (info==0) then
        ! compute inverse of B
        call F07MJF(UPLO,n,B,N,IPIV,WORK,INFO)
        InvM = B
      end if
      deallocate(IPIV)
      deallocate(WORK)
    end if  ! if ifail==0
    deallocate(A,A1,B,Z)
  end select
end subroutine MatrixInverse

! subroutine SphereToMatrix(phi,r,K,J,B,GradB_phi,GradB_r) 
! 
! B is a (K x J) matrix.
!
! For i<=K, B(:,i) is upper triangular and 
!   B(1,1)   = r(1)
!   B(1:2,2) = MapToCartesian(r(2),phi(1))
!   B(1:3,3) = MapToCartesian(r(3),phi(2:3)
!                  ...
! For i>K
!   B(:,i)   = MapToCartesian(r(i),phi(index(i))
!
!  phi       = (n x 1)     angles of spherical coordinate representation of B
!                      stored as a vector with n = K*(K-1)/2 + (K-1)*(J-K).
!  r         = (J x 1)     scaling factors for each column of B
!  B         = (K x J)     utility matrix
!  GradB_phi = (K x n) GradB_phi(:,i1) = gradient of col(i1) of B w.r.t. phi(i1)
!                      all other columns have zero gradients
!                      col(i1) is the column of B corresponding to phi(i1)
!  GradB_r   = (K x J) GradB_r(:,i1) = gradient of column i1 w.r.t. r(i1)
!                      all other columns have gradient = 0
! modification history
! --------------------
! 09DEC2012 LN  translated from matlab file SphereToMatrix.m
subroutine SphereToMatrix(phi,r,K,J,B,GradB_phi,GradB_r) 
  use nrtype
  implicit none
  real(dp),     intent(in)  :: phi(:)
  real(dp),     intent(in)  :: r(:)
  integer(i4b), intent(in)  :: K,J
  real(dp),     intent(out) :: B(:,:)
  real(dp),     intent(out),optional :: GradB_phi(:,:)
  real(dp),     intent(out),optional :: GradB_r(:,:)

  integer(i4b)              :: n,nphi
  integer(i4b)              :: i1
  integer(i4b)              :: index(K-1)
  real(dp)                  :: GradB0(K,K)

  ! Determine size of x and check it is compatible with dimensions of M
  n = size(phi)
  
  if (n .ne. (K*(K-1)/2+ (K-1)*(J-K))) then
    nPhi = K*(K-1)/2+K*(J-K)
    print *, 'J   K   n   nPhi'
    print *, J,K,n,nPhi
    print *, 'Error in SphereToMatrix. Mismatch between size of phi and dimensions of B.'
    
    stop
  end if
  
  B         = 0.0d0
  if (present(GradB_phi)) then
    GradB_phi = 0.0d0
  end if
  if (present(GradB_r)) then
    GradB_r   = 0.0d0
  end if

  B(1,1)       = r(1)
  if (present(GradB_r)) then
    GradB_r(1,1) = 1
  end if

  do i1 = 2,K
    ! B is upper triangular
    ! B(1:i1,i1) = is set by elements of phi(index) and r(i1)
   
    index(1:i1-1) = (i1-1)*(i1-2)/2 + (/1:i1-1/)
    GradB0 = 0.0d0
    call MapToCartesian(r(i1),phi(index(1:i1-1)),B(1:i1,i1),GradB0(1:i1,1:i1))
    if (present(GradB_phi)) then
      GradB_phi(1:i1,index(1:i1-1)) = transpose(GradB0(1:i1-1,1:i1))
    end if
    if (present(GradB_r)) then
      GradB_r(1:i1,i1)              = GradB0(i1,1:i1)
    end if
  end do
  do i1=K+1,J
    index = K*(K-1)/2 + (K-1)*(i1-K-1)+(/1:K-1/)
    GradB0 = 0.0d0
    call MapToCartesian(r(i1),phi(index),B(:,i1),GradB0)
    if (present(GradB_phi)) then
      GradB_phi(:,index) = transpose(GradB0(1:K-1,:))
    end if
    if (present(GradB_r)) then
      GradB_r(:,i1)      = GradB0(K,:)
    end if
  end do
end subroutine SphereToMatrix

! Convert (upper triangular) matrix to spherical coordinates
! size(C)   = K x J    upper triangular matrix with K<J
! size(D)   = K x 1    magnitude of each column
! size(PHI) = N x 1    spherical coordinate represenation of column
!                      N = K*(K-1)/2 + (J-K)*(K-1)
subroutine MatrixToSphere(C,D,PHI)
  use nrtype
  implicit none
  real(dp), intent(in)  :: C(:,:)
  real(dp), intent(out) :: D(:),PHI(:)
  integer(i4b)          :: K,J,i1,j1,j2

  integer(i4b), allocatable :: index(:)

  K = size(C,1)
  J = size(C,2)

  D   = 0.0d0
  PHI = 0.0d0
  D(1) = C(1,1)

  allocate(index(K-1))

  do i1=2,K
    ! index used to put result in correct location in PHI(:)
    index(1:i1-1) = (i1-1)*(i1-2)/2 + (/1:i1-1/)
    j1= (i1-1)*(i1-2)/2+1
    j2 = j1+i1-2
    call MapToSpherical(C(1:i1,i1),D(i1),PHI(j1:j2))
  end do

  if (J>K) then
    do i1=K+1,J
      j1 = K*(K-1)/2 + (K-1)*(i1-K-1) + 1
      j2 = j1+K-2
     ! index = K*(K-1)/2 + (K-1)*(i1-K-1) + (/1:K-1/)
      call MapToSpherical(C(:,i1),D(i1),PHI(j1:j2))
    end do
  end if 

  deallocate(index)
end subroutine MatrixToSphere

subroutine MapToCartesian(r,phi,x,dx)
  use nrtype
  implicit none
! Map (r,phi) from spherical coordinates to Cartesian
! gradient of map
! r   (1 x 1)
! phi (n-1 x 1)   phi(i1) is in [0,pi] for i1<n-1
!                 phi(n-1) is in [0,2*pi) 
!                 However, since we only want positive elements for x(n)
!                 we restrict phi(n-1) to be in [0,pi]
! x   (n x 1)
! dx(i,:) = gradient of x w.r.t. phi(i) if i<n
! dx(n,:)) = gradient of x w.r.t. r      if i==n
!
! Revision history
! 09DEC2012  LN  translated from matlab file MapToCartesian.m
  real(dp), intent(in) :: r
  real(dp), intent(in) :: phi(:)
  real(dp), intent(out) :: x(:)
  real(dp), optional,intent(out) :: dx(:,:)
  
  integer(i4b) :: i1,i2
  integer(i4b) :: n
  real(dp), allocatable :: s1(:),c1(:),s1A(:)

  n=size(phi)+1

  allocate(s1(n),c1(n),s1A(n))
  s1=(/1.0d0,dsin(phi)/)
  c1=(/dcos(phi),1.0d0/)

  x  = 0.0d0
  if (present(dx)) then
    dx = 0.0d0
  end if

  do i1=1,n
    x(i1)=r*c1(i1)*product(s1(1:i1))

    if (present(dx)) then
      if (i1<n) then
        dx(i1,i1) = -r*product(s1(1:i1+1))
        do i2=1,i1-1
          s1A(1:i1-1) = pack(s1,(/1:i1/) .ne. (i2+1))
          dx(i2,i1) = r*c1(i1)*product(s1A(1:i1-1))*c1(i2)
        end do
      elseif (i1==n) then
        do i2=1,i1-1
          s1A(1:i1-1) = pack(s1,(/1:i1/) .ne. (i2+1))
          dx(i2,i1) = r*c1(i1)*product(s1A(1:i1-1))*c1(i2)
        end do
        dx(i1,:) = x/r
      end if
    end if ! if (present(dx)) then
  end do
  deallocate(s1,c1,s1A)

end subroutine MapToCartesian

subroutine MapToSpherical(x,R,PHI)
  use nrtype
  implicit none
  real(dp), intent(in)  :: x(:)
  real(dp), intent(out) :: R,PHI(:)
  integer(i4b)          :: n,i1
  real(dp), allocatable :: x2(:)
  real(dp)              :: RTemp

  n = size(x,1)

  allocate(x2(n))

  x2 = x*x
  R = sqrt(sum(x2))
  PHI = 0.0d0
  do i1=1,n-1
    RTemp = sum(x2(i1:n))
    if (RTemp==0) then
      PHI(i1) = merge(0.0d0,pi_d,x(i1)>0)
    else
      ! DACOS = double precision arccos(theta)
      PHI(i1) = dacos(x(i1)/sqrt(RTemp))
    end if
  end do
  !phi(n-1) = dacos(x(n-1)/sqrt(sum(x2(n-1:n))))
  if (x(n)<0.0d0) then
    phi(n-1) = 2.0d0*pi_d - phi(n-1)
  end if
  deallocate(x2)
end subroutine MapToSpherical

subroutine MapToSpherical_old(x,R,PHI)
  use nrtype
  implicit none
  real(dp), intent(in)  :: x(:)
  real(dp), intent(out) :: R,PHI(:)
  integer(i4b)          :: n,i1
  real(dp), allocatable :: x2(:)

  n = size(x,1)

  allocate(x2(n))

  x2 = x*x
  R = sqrt(sum(x2))
  do i1=1,n-2
    PHI(i1) = atan2(sqrt(sum(x2(i1+1:n))),x(i1))
  end do
  phi(n-1) = atan2(x(n),x(n-1))

  deallocate(x2)
end subroutine MapToSpherical_old

subroutine FindColumn(i1,K,j1)
  use nrtype
  implicit none
! Find column in B corresponding to element i1 in phi
! That is, phi(i1) maps into column j1
! B is K x J matrix. For example, when K=4
!  B = [1 phi(1) phi(2) phi(4) phi(7) ... ;
!       0 x2     phi(3) phi(5) phi(8) ... ;
!       0 0      x3     phi(6) phi(9) ... ;
!       0 0      0      x4     x5     ...
   integer(i4b), intent(in)  :: i1,K
   integer(i4b), intent(out) :: j1
if (i1<=K*(K-1)/2) then
  j1 = floor(0.5d0*(3.0d0+dsqrt(8.0d0*real(i1)-7.0d0)))    
else
  j1 = ceiling( real(i1)/real(K-1) + 0.5d0*real(K))
end if

end subroutine FindColumn

! Compute inverse of A using LU decomposition
subroutine ComputeInverse_LU(K,A,InvA,ifail)
  use nrtype
  use nag_library, only : F04BAF
  implicit none
  integer(i4b), intent(in)  :: K
  real(dp),     intent(in)  :: A(:,:)
  real(dp),     intent(out) :: InvA(:,:)
  integer(i4b), intent(out) :: ifail

  integer(i4b)              :: i1
  integer(i4b), allocatable :: iPivot(:)
  real(dp)                  :: rCond,errBound
  real(dp), allocatable     :: TempA(:,:)

  allocate(TempA(K,K))

  InvA = 0.0d0
  do i1=1,K
    InvA(i1,i1) = 1.0d0
  end do

  ! F04BAF: inverse of InvC
  allocate(iPivot(K))
  ! value of A is changed by this routine
  TempA = A
  call F04BAF(K,K,TempA,K,iPivot,InvA,K,rCond,errBound,ifail)
  deallocate(iPivot)
  deallocate(TempA)
end subroutine ComputeInverse_LU

subroutine ComputeLQ(m,n,A,L,Q,ifail)
! Compute  LQ decomposition: A = L*Q
!    size(A) = (m x n)
!    size(L) = (m x n)   lower triangular
!    size(Q) = (n x n)   orthogonal
  use nrtype
  use nag_library, only : F08AHF,F08AJF
  implicit none
  integer(i4b), intent(in)  :: m
  integer(i4b), intent(in)  :: n
  real(dp),     intent(in)  :: A(:,:)
  real(dp),     intent(out) :: L(:,:)
  real(dp),     intent(out) :: Q(:,:)
  integer(i4b), intent(out) :: ifail

  real(dp), allocatable :: work(:)
  integer(i4b)          :: lwork
  real(dp), allocatable :: tau(:)
  integer(i4b)          :: i1
  real(dp), allocatable :: QTemp(:,:)

!  external F08AHF
!  external F08AJF

  LWork =  5*m
  allocate(tau(min(m,n)),work(LWork))
  allocate(QTemp(m,n))
  L = A
  call F08AHF(m,n,L,m,tau,work,LWork,ifail)
  QTemp = L
  call F08AJF(n,n,min(m,n),QTemp,m,tau,work,lwork,ifail)
  do i1=1,min(m,n)-1
    L(i1,i1+1:n) = 0.0d0
  end do
  Q = QTemp(1:n,1:n)
end subroutine ComputeLQ

function det(A)
! Compute determinant of A
  use nrtype
  use nag_library, only : F07ADF,F03BAF

  implicit none
  real(dp), intent(in) :: A(:,:)
  real(dp)             :: det

  real(dp), allocatable :: TempA(:,:)
  integer(i4b)              :: n
  integer(i4b), allocatable :: iPivot(:)
  integer(i4b)              :: ifail
  integer(i4b)              :: ID


  n = size(A,1)
  allocate(iPivot(n))
  allocate(TempA(n,n))

  ifail = 0
  TempA = A
  call F07ADF(n,n,TempA,n,ipivot,ifail)
  ifail = 0
  call F03BAF(n,TempA,n,iPivot,det,ID,ifail)
  det = det*2.0d0**ID

  deallocate(iPivot,TempA)
end function det

! Compute singular value decomposition of B
!    B = U*[diag(S) 0]*VT
subroutine ComputeSVD(B,U,S,VT)
  use nag_library, only : F08KBF  ! singular value decomposition
  use nrtype
  implicit none
  real(dp), intent(in)  :: B(:,:)
  real(dp), intent(out) :: U(:,:)
  real(dp), intent(out) :: S(:)
  real(dp), intent(out) :: VT(:,:)
  
  real(dp), allocatable :: B1(:,:)
  ! size(B)  = (d1 x d2)
  ! size(S)  = d1  (non zero elements of a (d1 x d2) matrix
  ! size(U)  = (d1 x d1)
  ! size(VT) = (d2,d2)
  !    B = U*[diag(S) 0]*VT
  integer(i4b) :: LWORK
  integer(i4b) :: d1,d2
  integer(i4b) :: ifail
  real(dp), allocatable ::work(:)
  d1=size(B,1)
  d2 = size(B,2)
  allocate(B1(d1,d2))
  B1 = B
  LWORK = max(1,3*d2+d1,5*d2)
  allocate(work(LWORK))

  ! F08KBF : singular value decomposition of B1
  ! F08KBF(JOBU,JOBVT,d1,d2,B,d1,S,U,d1,VT,d2,work,LWORK,ifail)
  ifail = -1
  call F08KBF('A','A',d1,d2,B1,d1,S,U,d1,VT,d2,work,lwork,ifail)

  deallocate(work,B1)
end subroutine ComputeSVD
end module NewTools
