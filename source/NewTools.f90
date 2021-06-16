! Revision history
! 2015AUG14 LN  add MatrixToSphere and MapToSpherical.
!               edit SphereToMatrix to work on rectangular matrix
module NewTools

contains

function sample(i1,nall,n,inputSeed) result(ix)
  use ConstantsModule
  use Ranking

  implicit none
  integer(i4b), intent(in) :: i1,nall,n
  integer(i4b), intent(in),optional :: inputSeed
  integer(i4b)             :: ix(n)

  integer(i4b)              :: nseed,nuni,n1,i2
  integer(i4b), allocatable :: seed(:)
  real(dp),     allocatable :: temp(:)
  integer(i4b), allocatable :: isort(:)

  if (i1==1) then
    call random_seed(size=nseed)
    allocate(seed(nseed))
    if (present(inputSeed)) then
      seed = inputSeed
    else
      seed = 45728474
    end if
    call random_seed(put=seed)
  end if

  allocate(temp(n))
  call random_number(temp)

  ix = ceiling(temp*nall)
  allocate(isort(n))
  call unirank(ix,isort,nuni)
  ix = ix(isort)
  do while (nuni<n)
    n1 = n-nuni
    call random_number(temp((/(i2,i2=1,n1)/)))
    ix((/(i2,i2=nuni+1,n)/)) = ceiling(nall*temp((/(i2,i2=1,n1)/)))
    call unirank(ix,isort,nuni)
    ix = ix(isort)
  end do

end function Sample


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
  use ConstantsModule
  implicit none
  real(dp),     intent(in)  :: phi(:)
  real(dp),     intent(in)  :: r(:)
  integer(i4b), intent(in)  :: K,J
  real(dp),     intent(out) :: B(:,:)
  real(dp),     intent(out),optional :: GradB_phi(:,:)
  real(dp),     intent(out),optional :: GradB_r(:,:)

  integer(i4b)              :: n,nphi
  integer(i4b)              :: i1,ix
  integer(i4b)              :: tempindex(K-1)
  real(dp)                  :: GradB0(K,K)
  real(dp), allocatable     :: temp1(:),temp2(:,:)
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

    tempindex((/(ix,ix=1,i1-1)/)) = (i1-1)*(i1-2)/2 + (/(ix,ix=1,i1-1)/)
    GradB0 = 0.0d0
    allocate(temp1(i1),temp2(i1,i1))
    call MapToCartesian(r(i1),phi(tempindex((/(ix,ix=1,i1-1)/))), &
                        temp1,temp2)
    B((/(ix,ix=1,i1)/),i1) = temp1
    GradB0((/(ix,ix=1,i1)/),(/(ix,ix=1,i1)/)) = temp2
    if (present(GradB_phi)) then
      GradB_phi((/(ix,ix=1,i1)/),tempindex((/(ix,ix=1,i1-1)/))) &
         = transpose(GradB0((/(ix,ix=1,i1-1)/),(/(ix,ix=1,i1)/)))
    end if
    if (present(GradB_r)) then
      GradB_r((/(ix,ix=1,i1)/),i1) = GradB0(i1,(/(ix,ix=1,i1)/))
    end if
    deallocate(temp1,temp2)
  end do
  do i1=K+1,J
    tempindex = K*(K-1)/2 + (K-1)*(i1-K-1)+(/(ix,ix=1,K-1)/)
    GradB0 = 0.0d0
    call MapToCartesian(r(i1),phi(tempindex),B(:,i1),GradB0)
    if (present(GradB_phi)) then
      GradB_phi(:,tempindex) = transpose(GradB0((/(ix,ix=1,K-1)/),:))
    end if
    if (present(GradB_r)) then
      GradB_r(:,i1)      = GradB0(K,:)
    end if
  end do
end subroutine SphereToMatrix

! Convert (upper triangular) matrix to spherical coordinates
! size(C)   = K x J    upper triangular matrix with K<J
! size(D)   = K x 1    magnitude of each column
! size(PHI) = N x 1    spherical coordinate representation of column
!                      N = K*(K-1)/2 + (J-K)*(K-1)
subroutine MatrixToSphere(C,D,PHI)
  use ConstantsModule
  implicit none
  real(dp), intent(in)  :: C(:,:)
  real(dp), intent(out) :: D(:),PHI(:)
  integer(i4b)          :: K,J,i1,j1,j2,ix

  real(dp),     allocatable :: temp(:)

  K = size(C,1)
  J = size(C,2)

  D   = 0.0d0
  PHI = 0.0d0
  D(1) = C(1,1)

  do i1=2,K
    ! (j1:j2) are locations in phi
    j1= (i1-1)*(i1-2)/2+1
    j2 = j1+i1-2
    allocate(temp(i1-1))
    call MapToSpherical(C((/(ix,ix=1,i1)/),i1),D(i1),temp)
    PHI((/(ix,ix=j1,j2)/)) = temp
    deallocate(temp)
  end do

  if (J>K) then
    allocate(temp(K-1))
    do i1=K+1,J
      j1 = K*(K-1)/2 + (K-1)*(i1-K-1) + 1
      j2 = j1+K-2
      call MapToSpherical(C(:,i1),D(i1),temp)
      PHI((/(ix,ix=j1,j2)/)) = temp
    end do
    deallocate(temp)
  end if

end subroutine MatrixToSphere

subroutine MapToCartesian(r,phi,x,dx)
  use ConstantsModule
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
  real(dp), intent(in)  :: r
  real(dp), intent(in)  :: phi(:)
  real(dp), intent(out) :: x(:)
  real(dp), optional,intent(out) :: dx(:,:)

  integer(i4b) :: i1,i2,ix
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
    x(i1)=r*c1(i1)*product(s1((/(ix,ix=1,i1)/)))

    if (present(dx)) then
      if (i1<n) then
        dx(i1,i1) = -r*product(s1((/(ix,ix=1,i1+1)/)))
        do i2=1,i1-1
          s1A((/(ix,ix=1,i1-1)/)) = pack(s1,(/(ix,ix=1,i1)/) .ne. (i2+1))
          dx(i2,i1) = r*c1(i1)*product(s1A((/(ix,ix=1,i1-1)/)))*c1(i2)
        end do
      elseif (i1==n) then
        do i2=1,i1-1
          s1A((/(ix,ix=1,i1-1)/)) = pack(s1,(/(ix,ix=1,i1)/) .ne. (i2+1))
          dx(i2,i1) = r*c1(i1)*product(s1A((/(ix,ix=1,i1-1)/)))*c1(i2)
        end do
        dx(i1,:) = x/r
      end if
    end if ! if (present(dx)) then
  end do
  deallocate(s1,c1,s1A)

end subroutine MapToCartesian

subroutine MapToSpherical(x,R,PHI)
  use ConstantsModule
  implicit none
  real(dp), intent(in)  :: x(:)
  real(dp), intent(out) :: R,PHI(:)
  integer(i4b)          :: n,i1,ix
  real(dp), allocatable :: x2(:)
  real(dp)              :: RTemp

  n = size(x,1)

  allocate(x2(n))

  x2 = x*x
  R = sqrt(sum(x2))
  PHI = 0.0d0
  do i1=1,n-1
    RTemp = sum(x2((/(ix,ix=i1,n)/)))
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
  use ConstantsModule
  implicit none
  real(dp), intent(in)  :: x(:)
  real(dp), intent(out) :: R,PHI(:)
  integer(i4b)          :: n,i1
  real(dp), allocatable :: x2(:)
  integer(i4b)          :: ix
  n = size(x,1)

  allocate(x2(n))

  x2 = x*x
  R = sqrt(sum(x2))
  do i1=1,n-2
    PHI(i1) = atan2(sqrt(sum(x2((/(ix,ix=i1+1,n)/)))),x(i1))
  end do
  phi(n-1) = atan2(x(n),x(n-1))

  deallocate(x2)
end subroutine MapToSpherical_old

subroutine FindColumn(i1,K,j1)
  use ConstantsModule
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

end module NewTools
