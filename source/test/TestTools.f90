program TestTools

use nrtype
use NewTools

implicit none

! Test SphereToMatrix and MatrixToSphere

integer(i4b) :: K,J,i1
real(dp), allocatable :: B(:,:),D(:),PHI(:)
real(dp), allocatable :: B1(:,:),D1(:),PHI1(:)

K = 2
J = 3

allocate(B(K,J))
allocate(B1(K,J))

allocate(D(J))
allocate(D1(J))

if (K==J) then
  allocate(PHI(K*(K-1)/2))
  allocate(PHI1(K*(K-1)/2))
elseif (K<J) then
  allocate(PHI(K*(K-1)/2 + (J-K)*(K-1)))
  allocate(PHI1(K*(K-1)/2 + (J-K)*(K-1)))
elseif (K>J) then
  print *,'Error. K>J is not allowed.'
end if


B(1,:) = (/1.0d0,5.768d0,564.34d0/)
B(2,:) = (/0.0d0,4.2d0,3.2d0/)

call MatrixToSPhere(B,D,PHI)

call SphereToMatrix(PHI,D,K,J,B1)
call MatrixToSphere(B1,D1,PHI1)

print *, 'B = '
do i1=1,K
  print *,B(i1,:)
end do

print *, 'B1 = '
do i1=1,K
  print *,B1(i1,:)
end do

print *,'D = '
print *,D
print *,'D1 = '
print *,D1

print *,'PHI = '
print *,PHI
print *,'PHI1 = '
print *,PHI1

deallocate(B)
deallocate(D)
deallocate(PHI)
deallocate(B1)
deallocate(D1)
deallocate(PHI1)



end program TestTools
