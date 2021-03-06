! Revision history
! 2018JUL26 LN  convert to gfortran syntax (/(ix,ix=1,k)/)
! 2015AUG31 LN  add MultiKron_rectangular
! 2015MAY24 LN  rename from multi_module to MultiModule.f90
! 08OCT2013 LN  change notation (y1 becomes y and y becomes x)
! 15sep2007 lpn Added multi1a and multi1b for use when subset of
!               calculations in multi1 need to be done repeatedly.
! 11sep2007 lpn edited and tested multi_inverse. Added comments.
! 11sep2007 Both multi1 and multi_inverse work properly.
!            multi was not tested on this occasion.          
!
! subroutine multi(y,x,d,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
!      compute y = x * kron(s10,...,s1)
! subroutine multi1(y,x,d,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
!      compute y = x * kron(s10,...,s1)
! subroutine multi1a(SS,n,d,s1,s2,s3,s4,s5,s6,s7,s8,s9)
!   Compute SS = kron(s(d-1),...,s1)
!   This is subset of mulit1 calculations
! subroutine multi1b(y,x,SS,sd)
!   This completes multi1a
!   call multi1a + call multi1b = call multi1
! subroutine MultiInverse(x,y,d,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
!    compute x to solve x = y * inv( kron(s10,...,s1) )

module MultiModule

use nrtype

! type defined to compute kronecker product of following form
!     y = x * kron(A(d)%b,...,A(1)%b)
type KronType
  real(dp), allocatable :: b(:,:)
end type

contains


! Compute y = x * kron(A%b(d),...,A%b(1))
! x       = (1 x N)          N = product(size(A(:)%b,1))
! A       = (d x 1)          structure, each element is matrix
! A(i1)%b = (n(i1) x r(i1))  matrix of size n(i1) by r(i1)
! d       = (1 x 1)          dimension of A
! y       = (1 x R)          R = product(size(A(:)%b,2))
subroutine MultiKron_rectangular(x,A,d,y)
  use nrtype
  implicit none
  real(dp), intent(in)       :: x(:)
  type(KronType), intent(in) :: A(:)
  integer(i4b), intent(in)   :: d
  real(dp), intent(out)      :: y(:)

  integer(i4b)              :: i1,i2,m,ny,ix
  integer(i4b), allocatable :: n(:),r(:)
  real(dp), allocatable     :: x0(:),y0(:)

  allocate(n(d),r(d))

  do i1=1,d
    n(i1) = size(A(i1)%b,1)
    r(i1) = size(A(i1)%b,2)
  end do

  m = product(n((/(ix,ix=2,d)/)))
  do i1=1,d
    if (i1==1) then
      ny = r(1)*m
      allocate(x0(n(1)*m))
      x0 = x
      allocate(y0(ny))
      y0 = 0.0d0
    else
      m = ny/n(i1)
      allocate(x0(ny))
      x0 = y0
      deallocate(y0)
      ny = r(i1) *m
      allocate(y0(ny))
      y0 = 0.0d0
    end if

    do i2=1,m
      y0(m*(/(ix,ix=0,r(i1)-1)/)+i2) = matmul(x0(n(i1)*(i2-1)+(/(ix,ix=1,n(i1))/)),A(i1)%b)
    end do
    deallocate(x0)
  end do

  y = y0

  deallocate(y0)
  deallocate(n,r)

end subroutine MultiKron_rectangular

subroutine multi(y,x,d,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
  ! y  = real(dp) (1 x nall)   : y = x * kron(sd,...,s1)
  ! x  = real(dp) (1 x nall)   input vector
  ! d  = integer(i4b) number of matrices in kronecker product. d<=10
  ! sj = real(dp) (nj x mj) matrix for j = 1,...,d
  use nrtype
  implicit none

  real(dp), intent(in) :: x(:),s1(:,:)
  real(dp), intent(out) :: y(:)
  integer(i4b), intent(in) :: d
  real(dp), optional, intent(in) :: s2(:,:),s3(:,:),s4(:,:),s5(:,:)
  real(dp), optional, intent(in) :: s6(:,:),s7(:,:),s8(:,:),s9(:,:),s10(:,:)

  real(dp), allocatable :: y0(:),stemp(:,:)
  integer(i4b) i1,nall,m,n,ix

  nall=size(x,1)
  allocate(y0(nall))

  ! Multiply by s1
  n=size(s1,1)
  m=nall/n
  allocate(stemp(n,n))
  y0=x
  y=0.d0
  !stemp=transpose(s1)
  stemp=s1
  do i1=1,m
    y((/m*(/(ix,ix=0,n-1)/)+i1/))=matmul(stemp,y0((/(ix,ix=n*(i1-1)+1,i1*n)/)))
  enddo
  deallocate(stemp)

  ! Multiply by s2
  if (d<2) then
    deallocate(y0)
    return
  else
    n=size(s2,1)
    m=nall/n
	allocate(stemp(n,n))
    y0=y
    y=0.d0
    !stemp=transpose(s2)
    stemp=s2
    do i1=1,m
      y((/m*(/(ix,ix=0,n-1)/)+i1/))=matmul(stemp,y0((/(ix,ix=n*(i1-1)+1,i1*n)/)))
    enddo
	deallocate(stemp)
  endif

  ! Multiply by s3
  if (d<3) then
    deallocate(y0)
    return
  else
    n=size(s3,1)
    m=nall/n
	allocate(stemp(n,n))
    y0=y
    y=0.d0
    !stemp=transpose(s3)
    stemp=s3
    do i1=1,m
      y((/m*(/(ix,ix=0,n-1)/)+i1/))=matmul(stemp,y0((/(ix,ix=n*(i1-1)+1,i1*n)/)))
    enddo
	deallocate(stemp)
  endif

  ! Multiply by s4
  if (d<4) then
    deallocate(y0)
    return
  else
    n=size(s4,1)
    m=nall/n
	allocate(stemp(n,n))
    y0=y
    y=0.d0
    !stemp=transpose(s4)
    stemp=s4
    do i1=1,m
      y((/m*(/(ix,ix=0,n-1)/)+i1/))=matmul(stemp,y0((/(ix,ix=n*(i1-1)+1,i1*n)/)))
    enddo
	deallocate(stemp)
  endif

  ! Multiply by s5
  if (d<5) then
    deallocate(y0)
    return
  else
    n=size(s5,1)
    m=nall/n
	allocate(stemp(n,n))
    y0=y
    y=0.d0
    !stemp=transpose(s5)
    stemp=s5
    do i1=1,m
      y((/m*(/(ix,ix=0,n-1)/)+i1/))=matmul(stemp,y0((/(ix,ix=n*(i1-1)+1,i1*n)/)))
    enddo
	deallocate(stemp)
  endif

  ! Multiply by s6
  if (d<6) then
    deallocate(y0)
    return
  else
    n=size(s6,1)
    m=nall/n
	allocate(stemp(n,n))
    y0=y
    y=0.d0
    !stemp=transpose(s6)
    stemp=s6
    do i1=1,m
      y((/m*(/(ix,ix=0,n-1)/)+i1/))=matmul(stemp,y0((/(ix,ix=n*(i1-1)+1,i1*n)/)))
    enddo
	deallocate(stemp)
  endif

  ! Multiply by s7
  if (d<7) then
    deallocate(y0)
    return
  else
    n=size(s7,1)
    m=nall/n
	allocate(stemp(n,n))
    y0=y
    y=0.d0
    !stemp=transpose(s7)
    stemp=s7
    do i1=1,m
      y((/m*(/(ix,ix=0,n-1)/)+i1/))=matmul(stemp,y0((/(ix,ix=n*(i1-1)+1,i1*n)/)))
    enddo
	deallocate(stemp)
  endif

  ! Multiply by s8
  if (d<8) then
    deallocate(y0)
    return
  else
    n=size(s8,1)
    m=nall/n
	allocate(stemp(n,n))
    y0=y
    y=0.d0
    !stemp=transpose(s8)
    stemp=s8
    do i1=1,m
      y((/m*(/(ix,ix=0,n-1)/)+i1/))=matmul(stemp,y0((/(ix,ix=n*(i1-1)+1,i1*n)/)))
    enddo
	deallocate(stemp)
  endif

  ! Multiply by s9
  if (d<9) then
    deallocate(y0)
    return
  else
    n=size(s9,1)
    m=nall/n
	allocate(stemp(n,n))
    y0=y
    y=0.d0
    !stemp=transpose(s9)
    stemp=s9
    do i1=1,m
      y((/m*(/(ix,ix=0,n-1)/)+i1/))=matmul(stemp,y0((/(ix,ix=n*(i1-1)+1,i1*n)/)))
    enddo
	deallocate(stemp)
  endif

  ! Multiply by s10
  if (d<10) then
    deallocate(y0)
    return
  else
    n=size(s10,1)
    m=nall/n
	allocate(stemp(n,n))
    y0=y
    y=0.d0
    !stemp=transpose(s10)
    stemp=s10
    do i1=1,m
      y((/m*(/(ix,ix=0,n-1)/)+i1/))=matmul(stemp,y0((/(ix,ix=n*(i1-1)+1,i1*n)/)))
    enddo
	deallocate(stemp)
  endif
  deallocate(y0)
end subroutine multi

subroutine multi1(y,x,d,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
  ! Compute y such that y = x*kron(sd,...,s1)
  ! y  = real(dp) (1 x 1)   : y = x * kron(sd,...,s1)
  ! x  = real(dp) (1 x nall)   input vector
  ! d  = integer(i4b) number of matrices in kronecker product. d<=10
  ! sj = real(dp) (nj x 1) vector for j = 1,...,d
  use nrtype
  use ToolsModule, only : kron1
  implicit none

  real(dp), intent(in) :: x(:),s1(:)
  real(dp), intent(out) :: y
  integer(i4b), intent(in) :: d
  real(dp), optional, intent(in) :: s2(:),s3(:),s4(:),s5(:)
  real(dp), optional, intent(in) :: s6(:),s7(:),s8(:),s9(:),s10(:)

  real(dp), allocatable :: ss(:),sstemp(:),sstemp2(:)
  integer(i4b) nall,n,ix

  if (d>10) then
    print *, 'Error in "multi1". d is greater than 10.'
    stop
  endif
  
  nall=size(x)
  allocate(ss(nall))
  n=size(s1)
  ss((/(ix,ix=1,n)/))=s1

  if (d<2) then
    y=dot_product(x,ss)	! Multiply x*s1
    deallocate(ss)
    return
  else
    n=n*size(s2)
    allocate(sstemp(n))
    call kron1(s2,s1,sstemp) ! Compute kron(s2,s1)
    ss((/(ix,ix=1,n)/)) = sstemp
    deallocate(sstemp)
  endif
  
  if (d<3) then
    y=dot_product(x,ss)					! Multiply y*kron(s2,s1)
	deallocate(ss)
	return
  else	  
    allocate(sstemp(n))
	sstemp=ss((/(ix,ix=1,n)/))
    n=n*size(s3)
    allocate(sstemp2(n))
	call kron1(s3,sstemp,sstemp2)			! Compute kron(s3,sstemp)
    ss((/(ix,ix=1,n)/)) = sstemp2
	deallocate(sstemp,sstemp2)
  endif

  if (d<4) then
    y=dot_product(x,ss)					! Multiply y*ss
	deallocate(ss)
	return
  else
    allocate(sstemp(n))	
    sstemp=ss((/(ix,ix=1,n)/))	
    n=n*size(s4)
    allocate(sstemp2(n))
    call kron1(s4,sstemp,sstemp2) ! Compute kron(s4,sstemp)
    ss((/(ix,ix=1,n)/)) = sstemp2
    deallocate(sstemp,sstemp2)
   endif

  if (d<5) then
    y=dot_product(x,ss)					! Multiply y*ss
	deallocate(ss)
	return
  else
    allocate(sstemp(n))	
	sstemp=ss((/(ix,ix=1,n)/))	
	n=n*size(s5)
    allocate(sstemp2(n))
	call kron1(s5,sstemp,sstemp2)			! Compute kron(s5,sstemp)
    ss((/(ix,ix=1,n)/)) = sstemp2
	deallocate(sstemp,sstemp2)
   endif

  if (d<6) then
    y=dot_product(x,ss)					! Multiply y*ss
	deallocate(ss)
	return
  else
    allocate(sstemp(n))	
	sstemp=ss((/(ix,ix=1,n)/))	
	n=n*size(s6)
    allocate(sstemp2(n))
	call kron1(s6,sstemp,sstemp2)			! Compute kron(s6,sstemp)
    ss((/(ix,ix=1,n)/)) = sstemp2
	deallocate(sstemp,sstemp2)
   endif

  if (d<7) then
    y=dot_product(x,ss)					! Multiply y*ss
	deallocate(ss)
	return
  else
    allocate(sstemp(n))	
	sstemp=ss((/(ix,ix=1,n)/))	
	n=n*size(s7)
    allocate(sstemp2(n))
	call kron1(s7,sstemp,sstemp2)			! Compute kron(s7,sstemp)
    ss((/(ix,ix=1,n)/)) = sstemp2
	deallocate(sstemp,sstemp2)
   endif

  if (d<8) then
    y=dot_product(x,ss)					! Multiply y*ss
	deallocate(ss)
	return
  else
    allocate(sstemp(n))	
	sstemp=ss((/(ix,ix=1,n)/))	
	n=n*size(s8)
    allocate(sstemp2(n))
	call kron1(s8,sstemp,sstemp2)			! Compute kron(s8,sstemp)
    ss((/(ix,ix=1,n)/)) = sstemp2
	deallocate(sstemp,sstemp2)
   endif

  if (d<9) then
    y=dot_product(x,ss)					! Multiply y*ss
	deallocate(ss)
	return
  else
    allocate(sstemp(n))	
	sstemp=ss((/(ix,ix=1,n)/))	
	n=n*size(s9)
    allocate(sstemp2(n))
	call kron1(s9,sstemp,sstemp2)		! Compute kron(s9,sstemp)
    ss((/(ix,ix=1,n)/)) = sstemp2
	deallocate(sstemp,sstemp2)
   endif

  if (d<10) then
    y=dot_product(x,ss)					! Multiply y*ss
	deallocate(ss)
	return
  else
    allocate(sstemp(n))	
	sstemp=ss((/(ix,ix=1,n)/))	
	n=n*size(s10)
    allocate(sstemp2(n))
	call kron1(s10,sstemp,sstemp2)		! Compute kron(s10,sstemp)
    ss((/(ix,ix=1,n)/)) = sstemp2
	deallocate(sstemp,sstemp2)
   endif

end subroutine multi1

subroutine multi1a(ss,d,s1,s2,s3,s4,s5,s6,s7,s8,s9)
  ! Compute ss = kron(s(d),...,s(1))
  ! First step in computing y = x * kron(s(d+1),...,s1)
  ! Second step is completed in multi1b
  !
  ! call multi1a(ss,2,s1,s2)
  ! call multi1b(y,x,ss,s3)
  !
  ! d<=10
  ! only a maximum of 9 vectors (s1,...,s9) are input
  ! y  = real(dp) (1 x 1)   = x * kron(sd,...,s1)
  ! x  = real(dp) (1 x nall)   input vector
  ! d  = integer(i4b) number of matrices in kronecker product in first step. d<=10
  ! sj = real(dp) (nj x 1) vector for j = 1,...,d
  use nrtype
  use ToolsModule, only : kron1
  implicit none

  real(dp), intent(in)           :: s1(:)
  real(dp), intent(out)          :: ss(:)
  integer(i4b), intent(in)       :: d
  real(dp), optional, intent(in) :: s2(:),s3(:),s4(:),s5(:)
  real(dp), optional, intent(in) :: s6(:),s7(:),s8(:),s9(:)

  real(dp), allocatable :: sstemp(:),sstemp2(:)
  integer(i4b) ntemp,ix

  if (d>10) then
	print *, 'Error in "multi1". d is greater than 10.'
	stop
  endif
  
  ntemp=size(s1)
  ss((/(ix,ix=1,ntemp)/))=s1

  if (d<2) then
    return
  else
    ntemp=ntemp*size(s2)
    allocate(sstemp2(ntemp))
    call kron1(s2,s1,sstemp2)				! Compute kron(s2,s1)
    ss((/(ix,ix=1,ntemp)/)) = sstemp2
    deallocate(sstemp2)
  endif
  	
  if (d<3) then
    return
  else	  
    allocate(sstemp(ntemp))
	sstemp=ss((/(ix,ix=1,ntemp)/))
    ntemp=ntemp*size(s3)
    allocate(sstemp2(ntemp))
	call kron1(s3,sstemp,sstemp2)			! Compute kron(s3,sstemp)
    ss((/(ix,ix=1,ntemp)/)) = sstemp2
	deallocate(sstemp,sstemp2)
  endif

  if (d<4) then
	return
  else
    allocate(sstemp(ntemp))	
	sstemp=ss((/(ix,ix=1,ntemp)/))	
	ntemp=ntemp*size(s4)
    allocate(sstemp2(ntemp))
	call kron1(s4,sstemp,sstemp2)			! Compute kron(s4,sstemp)
    ss((/(ix,ix=1,ntemp)/)) = sstemp2
	deallocate(sstemp,sstemp2)
   endif

  if (d<5) then
	return
  else
    allocate(sstemp(ntemp))	
	sstemp=ss((/(ix,ix=1,ntemp)/))	
	ntemp=ntemp*size(s5)
    allocate(sstemp2(ntemp))
	call kron1(s5,sstemp,sstemp2)			! Compute kron(s5,sstemp)
    ss((/(ix,ix=1,ntemp)/)) = sstemp2
	deallocate(sstemp,sstemp2)
   endif

  if (d<6) then
    return
  else
    allocate(sstemp(ntemp))	
	sstemp=ss((/(ix,ix=1,ntemp)/))	
	ntemp=ntemp*size(s6)
    allocate(sstemp2(ntemp))
	call kron1(s6,sstemp,sstemp2)			! Compute kron(s6,sstemp)
    ss((/(ix,ix=1,ntemp)/)) = sstemp2
	deallocate(sstemp,sstemp2)
   endif

  if (d<7) then
	return
  else
    allocate(sstemp(ntemp))	
	sstemp=ss((/(ix,ix=1,ntemp)/))	
	ntemp=ntemp*size(s7)
    allocate(sstemp2(ntemp))
	call kron1(s7,sstemp,sstemp2)			! Compute kron(s7,sstemp)
    ss((/(ix,ix=1,ntemp)/)) = sstemp2
	deallocate(sstemp,sstemp2)
   endif

  if (d<8) then
	return
  else
    allocate(sstemp(ntemp))	
	sstemp=ss((/(ix,ix=1,ntemp)/))	
	ntemp=ntemp*size(s8)
    allocate(sstemp2(ntemp))
	call kron1(s8,sstemp,sstemp2)			! Compute kron(s8,sstemp)
    ss((/(ix,ix=1,ntemp)/)) = sstemp2
	deallocate(sstemp,sstemp2)
   endif

  if (d<9) then
	return
  else
    allocate(sstemp(ntemp))	
	sstemp=ss((/(ix,ix=1,ntemp)/))	
	ntemp=ntemp*size(s9)
    allocate(sstemp2(ntemp))
	call kron1(s9,sstemp,sstemp2)			! Compute kron(s9,sstemp)
    ss((/(ix,ix=1,ntemp)/)) = sstemp2
	deallocate(sstemp,sstemp2)
   endif

end subroutine multi1a

subroutine multi1b(y,x,ss,sd)
  ! Complete computations from multi1a
  ! Compute y = x*kron(sd,ss)
  ! y = real(dp) (1 x 1)   = x * kron(sd,ss)
  ! x  = real(dp) (nall x 1)   input vector
  ! ss = real(dp) (nall/nd x 1) input vector
  ! sd = real(dp) (nd x 1) vector 
  use nrtype
  use ToolsModule, only : kron1
  implicit none

  real(dp), intent(in)  :: x(:),ss(:),sd(:)
  real(dp), intent(out) :: y
  real(dp)              :: stemp(size(x))

  call kron1(sd,ss,stemp)
  y=dot_product(x,stemp)
  
end subroutine multi1b

subroutine MultiInverse(x,y,d,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
! subroutine MultiInverse(x,y,d,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10)
! 
! Find x that solves x = y * inv( kron(sd,...,s1) )
!
! x     = real(dp) (1 x n1*...*nd) solution to linear equation
! y     = real(dp) (1 x n1*...*nd) vector of data
! d		= integer(i4b) number of matrices in kronecker product
! sj    = real(dp) (nj x nj) for j=1,...,d. matrices to be inverted
  use nrtype
  use nr, only : ludcmp,lubksb

  implicit none
  real(dp),intent(in)  :: y(:)
  real(dp), intent(out) :: x(:)
  integer(i4b), intent(in) :: d
  real(dp), intent(in) :: s1(:,:)
  real(dp), optional, intent(in) :: s2(:,:)
  real(dp), optional, intent(in) :: s3(:,:),s4(:,:),s5(:,:),s6(:,:)
  real(dp), optional, intent(in) :: s7(:,:),s8(:,:),s9(:,:),s10(:,:)
  real(dp) ipivot
  real(dp), allocatable :: y0(:),stemp(:,:),ytemp(:)
  integer(i4b), allocatable :: pivot(:)
  integer(i4b) i1,n,m,nall,ix
  
  if (d<1 .or. d>10) then
    print *, 'Error in MultiInverse. d must be greater than 1 or less than 11.'
	print *, 'Current value of d is ',d
	stop
  endif

  nall=size(y)
  allocate(y0(nall))

  ! Invert s1
  ! To find x that solves b=a*x, use the following:
  ! call ludcmp(a,pivot,ipivot)
  ! x=b
  ! call lubksb(a,pivot,x)
  ! Below we want y1 that solves: y=y1*s
  ! So, to do this we first compute stemp=transpose(s)
  
  n=size(s1,1)
  allocate(pivot(n),stemp(n,n),ytemp(n))          ! For use in LU factorization
  m=nall/n
  y0=y
  x=0.d0
  stemp=transpose(s1)
  call ludcmp(stemp,pivot,ipivot)

  do i1=1,m
    ytemp=y0((/(ix,ix=n*(i1-1)+1,i1*n)/))
    call lubksb(stemp,pivot,ytemp)			! Solve x(section(1,i1)) =  s1*y0(section(1,i1))
	x(m*(/(ix,ix=0,n-1)/)+i1)=ytemp			! Solve for x.
  enddo
  deallocate(pivot,ytemp,stemp)

  ! Invert s2
  if (d<2) then
    return
  else
    n=size(s2,1)
    allocate(pivot(n),stemp(n,n),ytemp(n))          ! For use in LU factorization
    m=nall/n
    y0=x
    stemp=transpose(s2)
    call ludcmp(stemp,pivot,ipivot)
	
	do i1=1,m
      ytemp=y0((/(ix,ix=n*(i1-1)+1,i1*n)/))			! Solve x(section(2,i1)) =  s2*y0(section(2,i1))
	  call lubksb(stemp,pivot,ytemp)
	  x((/m*(/(ix,ix=0,n-1)/)+i1/))=ytemp			! Solve for x.
    enddo
    deallocate(pivot,stemp,ytemp)
  endif
  
  ! Invert s3
  if (d<3) then
    return
  else
    n=size(s3,1)
    allocate(pivot(n),stemp(n,n),ytemp(n))			! For use in LU factorisation
    m=nall/n
    y0=x
	stemp=transpose(s3)
    call ludcmp(stemp,pivot,ipivot)
	do i1=1,m
      ytemp=y0((/(ix,ix=n*(i1-1)+1,i1*n)/))			! Solve x(section(3,i1)) =  s3*y0(section(3,i1))
	  call lubksb(stemp,pivot,ytemp)
	  x((/m*(/(ix,ix=0,n-1)/)+i1/))=ytemp			! Solve for x.
    enddo
    deallocate(pivot,stemp,ytemp)
  endif

  ! Invert s4
  if (d<4) then
    return
  else
    n=size(s4,1)
    allocate(pivot(n),stemp(n,n),ytemp(n))			! For use in LU factorisation
    m=nall/n
    y0=x
	stemp=transpose(s4)
    call ludcmp(stemp,pivot,ipivot)
	do i1=1,m
      ytemp=y0((/(ix,ix=n*(i1-1)+1,i1*n)/))			! Solve x(section(4,i1)) =  s4*y0(section(4,i1))
	  call lubksb(stemp,pivot,ytemp)
	  x((/m*(/(ix,ix=0,n-1)/)+i1/))=ytemp			! Solve for x.
    enddo
    deallocate(pivot,stemp,ytemp)
  endif

  ! Invert s5
  if (d<5) then
    return
  else
    n=size(s5,1)
    allocate(pivot(n),stemp(n,n),ytemp(n))			! For use in LU factorisation
    m=nall/n
    y0=x
	stemp=transpose(s5)
	call ludcmp(stemp,pivot,ipivot)
	do i1=1,m
      ytemp=y0((/(ix,ix=n*(i1-1)+1,i1*n)/))			! Solve x(section(5,i1)) =  s5*y0(section(5,i1))
	  call lubksb(stemp,pivot,ytemp)
	  x((/m*(/(ix,ix=0,n-1)/)+i1/))=ytemp			! Solve for x.
    enddo
    deallocate(pivot,stemp,ytemp)
  endif

  ! Invert s6
  if (d<6) then
    return
  else
    n=size(s6,1)
    allocate(pivot(n),stemp(n,n),ytemp(n))			! For use in LU factorisation
    m=nall/n
    y0=x
	stemp=transpose(s6)
	call ludcmp(stemp,pivot,ipivot)
	do i1=1,m
      ytemp=y0((/(ix,ix=n*(i1-1)+1,i1*n)/))			! Solve x(section(6,i1)) =  s6*y0(section(6,i1))
	  call lubksb(stemp,pivot,ytemp)
	  x((/m*(/(ix,ix=0,n-1)/)+i1/))=ytemp			! Solve for x.
    enddo
    deallocate(pivot,stemp,ytemp)
  endif
  
  ! Invert s7
  if (d<7) then
    return
  else
    n=size(s7,1)
    allocate(pivot(n),stemp(n,n),ytemp(n))			! For use in LU factorisation
    m=nall/n
    y0=x
	stemp=transpose(s7)
    call ludcmp(stemp,pivot,ipivot)
	do i1=1,m
      ytemp=y0((/(ix,ix=n*(i1-1)+1,i1*n)/))			! Solve x(section(7,i1)) =  s7*y0(section(7,i1))
	  call lubksb(stemp,pivot,ytemp)
	  x((/m*(/(ix,ix=0,n-1)/)+i1/))=ytemp			! Solve for x.
    enddo
    deallocate(pivot,stemp,ytemp)
  endif

  ! Invert s8
  if (d<8) then
    return
  else
    n=size(s8,1)
    allocate(pivot(n),stemp(n,n),ytemp(n))			! For use in LU factorisation
    m=nall/n
    y0=x
	stemp=transpose(s8)
    call ludcmp(stemp,pivot,ipivot)
	do i1=1,m
      ytemp=y0((/(ix,ix=n*(i1-1)+1,i1*n)/))			! Solve x(section(8,i1)) =  s8*y0(section(8,i1))
	  call lubksb(stemp,pivot,ytemp)
	  x((/m*(/(ix,ix=0,n-1)/)+i1/))=ytemp			! Solve for x.
    enddo
    deallocate(pivot,stemp,ytemp)
  endif

  ! Invert s9
  if (d<9) then
    return
  else
    n=size(s9,1)
    allocate(pivot(n),stemp(n,n),ytemp(n))			! For use in LU factorisation
    m=nall/n
    y0=x
	stemp=transpose(s9)
	call ludcmp(stemp,pivot,ipivot)
	do i1=1,m
      ytemp=y0((/(ix,ix=n*(i1-1)+1,i1*n)/))			! Solve x(section(9,i1)) =  s9*y0(section(9,i1))
	  call lubksb(stemp,pivot,ytemp)
	  x((/m*(/(ix,ix=0,n-1)/)+i1/))=ytemp			! Solve for x.
    enddo
    deallocate(pivot,stemp,ytemp)
  endif

  ! Invert s10
  if (d<10) then
    return
  else
    n=size(s10,1)
    allocate(pivot(n),stemp(n,n),ytemp(n))			! For use in LU factorisation
    m=nall/n
    y0=x
	stemp=transpose(s10)
    call ludcmp(stemp,pivot,ipivot)
	do i1=1,m
      ytemp=y0((/(ix,ix=n*(i1-1)+1,i1*n)/))			! Solve x(section(10,i1)) =  s10*y0(section(10,i1))
	  call lubksb(stemp,pivot,ytemp)
	  x((/m*(/(ix,ix=0,n-1)/)+i1/))=ytemp			! Solve for x.
    enddo
    deallocate(pivot,stemp,ytemp)
  endif
end subroutine MultiInverse

end module MultiModule
