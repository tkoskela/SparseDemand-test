program TestTools

use ConstantsModule
use NewTools

implicit none

! Test SphereToMatrix and MatrixToSphere

integer(i4b) :: K,J,i1,i2,nargs,nx
real(dp), allocatable :: B(:,:),D(:),PHI(:)
real(dp), allocatable :: B1(:,:),D1(:),PHI1(:)
character(len=5)      :: ctemp

print *,"Test MatrixToSphereSphere:"
print *,"Routine to convert upper-triangular matrix to spherical coordinates."

nargs = IARGC()
if (nargs==0) then
  print *,"Zero input arguments. Assuming default matrix size (K=2,J=3)."
  print *,"To set matrix size (K,J), type ./TestTools.dbg K J"
  K = 2
  J = 3
else if (nargs==1) then
  call GetArg(1,cTemp)
  read(cTemp,'(i3)') K
  J = K+1
elseif (nargs==2) then
  call getarg(1,cTemp)
  read(ctemp,'(i3)') K
  call getarg(2,ctemp)
  read(ctemp,'(i3)') J
end if

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

if (K==2 .and. J==3) then
  B(1,:) = (/1.0d0,5.768d0,564.34d0/)
  B(2,:) = (/0.0d0,4.2d0,3.2d0/)
else
  call random_number(B)
  B = 10*B-5
  do i1=2,K
    B(i1,1:i1-1) = 0.0d0
  end do
end if

call MatrixToSphere(B,D,PHI)

call SphereToMatrix(PHI,D,K,J,B1)
call MatrixToSphere(B1,D1,PHI1)

print *, 'B         = Initial (K,J) matrix.'
print *, "(D,PHI)   = MatrixToSphere(B)."
print *, "B1        = SphereToMatrix(PHI,D)."
print *, "(D1,PHI1) = MatrixToSphere(B1)."
do i1=1,K
  print *,B(i1,:)
end do

print *, "Test B1-B = 0."
print *, 'B1-B = '
do i1=1,K
  print *,B1(i1,:)-B(i1,:)
end do

print *, "Test D = D1."
print *,'D = '
print *,D
print *,'D1 = '
print *,D1

print *, "Test PHI = PHI1."
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
