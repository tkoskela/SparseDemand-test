! ToolsModule
!
! Functions
!-----------------------------------------------------------------------------
! function gammln(xx)			          Log Gamma function
! function dmatmul(A,B,n1,n2,n3)          Compute A*B
! function trans1(a,n1,n2)		          Transpose of a. a has dimensions (n1 x n2).
! function kron2(x,n1,n2,nq)              Kronecker of x and identity
! function normcdf_mkl(x,nx) result(y)    Compute normal CDF using Intel MKL
! function normcdf(x,nx) result(y)        Compute normal CDF using erfc
! function InverseNormal_mkl(p,n,ifault) result(x)   Compute inverse CDF normal
!                                                    using Intel MKL
! function InverseNormal(P,IFAULT) result(x)  Compute normal inverse CDF
!                                             (less accurate than MKL)
! function normpdf(x,n) result(y)         Normal pdf.
! function InvertLower(x,d) result(InvX)  Invert lower triangular matrix.
!
!-----------------------------------------------------------------------------
! Subroutines
!-----------------------------------------------------------------------------
! subroutine bsolve1L(x,A,n,b)		Solve Ax=b using backsubstitution, A lower diagonal
! subroutine bsolve1U(x,A,n,b)          Solve Ax=b using backsubstitution, A upper diagonal
! subroutine bsolve2L(x,A,k,b) 	        Solve A*x=b, backsubstitution, A contains lower diagonal matrix
! subroutine bsolve2U(x,A,k,b)          Solve A*x=b using backsubstitution, A contains upper diagonal matrix
! subroutine chol1(L,A,n)		Cholesky decomposition of matrix A
! subroutine cholinv1(Linv,L,n)		Given L, lower diagonal matrix, compute L inverse
! subroutine cholinv2(Linv,L,n)		Given L, an n*(n+1)/2 lower diagonal matrix, compute L inverse
! subroutine invert1(Kinv,det,K,n)	Invert K using Cholesky deconp, K must be positive definite symmetric
! subroutine invert2(Kinv,det,K,n)	Invert K. K is the lower diagonal of positive definite symmetric
! subroutine chol2(L,A,n)		Cholesky decomposition of matrix A
! subroutine gaucheb(x,w)		Nodes and weights for Gauss-Chebyshev integration
! subroutine gauleg(x,w,n)		Nodes and weights for Gauss-Legendre integration
! subroutine gaulag(x,w,alpha0)		Nodes and weights for Gauss-Laguerre integration
! subroutine numdif00(func,nx,x0,df,iuser,user)
! SUBROUTINE CONFUN(MODE,NCNLN,N,LDCJ,NEEDC,X,C,CJAC,NSTATE,IUSER,USER)
! subroutine kron1(a,b,c)
! subroutine ikron1(a,b,c)
! subroutine linspace(a,b,n,x)
! subroutine kron(a,b,c)
! subroutine ikron(a,b,c)
!-----------------------------------------------------------------------------
module ToolsModule

implicit none

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Returns the value of ln[GAM(x)) for xx>0
! xx         =   real(dp) (1 x 1) input value
! gammln    =   real(dp) (1 x 1) output value
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gammln(xx)
  use nrtype
  implicit none
  real(dp), intent(in) :: xx
  real(dp) gammln
  real(dp), parameter :: cof(6)=(/76.18009172947146d0,-86.50532032941677d0, &
                                 24.01409824083091d0,-1.231739572450155d0, &
                                 0.1208650973866179d-2,-0.5395239384953d-5/)
  real(dp), parameter :: stp=2.5066282746310005d0
  real(dp) x,y,tmp,ser                               
  integer(i4b) j          
  
  x=xx
  y=xx
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  do j=1,6
    y=y+1.d0
    ser=ser+cof(j)/y
  end do
  gammln=tmp+log(stp*ser/x)
end function gammln    

function dmatmul(A,B,n1,n2,n3)
  use nrtype
  implicit none
! function dmatmul(A,B,n1,n2,n3)
! Matrix multiplication C=dmatmul(A,B)
! A(n1,n2)	real(8): matrix
! B(n2,n3)	real(8): matrix
! n1,n2,n3	integers: dimensions of matrices
! C(n1,n3)	real(8):  C=A*B
  integer, intent(in) :: n1,n2,n3
  real(dp), intent(in) :: A(n1,n2), B(n2,n3)
  real(dp)             :: dmatmul(n1,n3)
  integer i,j,k
  dmatmul=0.0d0
  do i=1,n3
    do j=1,n2
      do k=1,n1
        dmatmul(k,i)=dmatmul(k,i)+A(k,j)*B(j,i)
      end do
    end do
  end do
end function dmatmul

function trans1(a,n1,n2)
! Compute transpose of matrix a. a has dimensions (n1 x n2).
  use nrtype
  implicit none
  integer(i4b), intent(in) :: n1,n2
  real(dp), intent(in) :: a(n1,n2)
  real(dp) trans1(n2,n1)
  real(dp) temp
  integer(i4b) i1,i2
 
  do i1=1,n1
    if (i1<=n2) then
      trans1(i1,i1)=a(i1,i1)
    end if
    do i2=i1+1,n2
      if (i1<=n2 .and. i2<=n1) then        
        trans1(i1,i2)=a(i2,i1)
      end if
      trans1(i2,i1)=a(i1,i2)
    end do
  end do
end function trans1  

function kron2(x,n1,n2,nq)
  ! kron2 creates nq/n1 copies of the matrix x and creates the (nq x n2) matrix with values equal to
  !   [ x(1,1) x(1,2)...x(1,n2) ]
  !     .
  !     .
  !   [ x(1,1) x(1,2)...x(1,n2) ]
  !   [ x(2,1) x(2,2)...x(2,n2) ]
  !     .
  !     .
  !   [ x(n1,1) x(n1,2)...x(n1,n2) ]
  use nrtype
  implicit none  
  integer(i4b), intent(in) :: n1,n2,nq
  real(dp), intent(in) :: x(n1,n2)
  real(dp) kron2(nq,n2)
  real(dp), allocatable :: temp1a(:),temp1b(:),temp2a(:,:),temp2b(:,:),temp2c(:,:)
  
  allocate(temp1a(n1*n2),temp1b(n2*nq),temp2a(n2,n1),temp2b(n1*n2,nq/n1),temp2c(n2,nq))
  temp2a=transpose(x)
  temp1a=reshape(temp2a,(/n1*n2/))
  temp2b=spread(temp1a,2,nq/n1)
  temp1b=reshape(temp2b,(/nq*n2/))
  temp2c=reshape(temp1b,(/n2,nq/))
  kron2=transpose(temp2c)  
  deallocate(temp1a,temp1b,temp2a,temp2b,temp2c)
end function kron2      

function normcdf_mkl(x,nx) result(y)
  use nrtype
  implicit none
  integer(i4b), intent(in) :: nx
  real(dp),     intent(in) :: x(nx)
  real(dp)                 :: y(nx)

  include "mkl_vml.fi"
  
  ! normal cdf from intel_MKL
  call vdcdfnorm(nx,x,y)

end function normcdf_mkl

function normcdf(x,nx) result(y)
  use nrtype
  use nr, only : erfc,erf
  implicit none
  integer(i4b), intent(in) :: nx
  real(dp),     intent(in) :: x(nx)
  real(dp)                 :: y(nx)

  y = merge(1.0d0 - 0.5d0*erfc(x/sqrt(2.0d0)), &
            0.5d0+0.5d0*erf(x/sqrt(2.0d0)),    &
            x<=0.0d0)

end function normcdf

function InverseNormal_mkl(p,n,ifault) result(x)
  use nrtype
  implicit none
  integer(i4b), intent(out) :: ifault
  real(dp),     intent(in)   :: p(n)
  real(dp)                   :: x(n)
  integer(i4b)               :: n

  include "mkl_vml.fi"
  ifault = 0 
  ! inverse normal function from intel_mkl: double precision
  call vdcdfnorminv(n,p,x)
  !x = p*p
  
end function InverseNormal_mkl


FUNCTION InverseNormal(P,IFAULT) result(x)

use nrtype
implicit none
integer(i4b), intent(out) :: ifault
real(dp),     intent(in)   :: p
real(dp)                   :: x

REAL(dp) :: ZERO,ONE,HALF,SPLIT1,SPLIT2, CONST1, CONST2, &
A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6, B7, &
C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5, D6, D7, &
E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7
!real(dp), allocatable :: Q(:),R(:)
real(dp)          :: Q,R

!
PARAMETER (ZERO = 0.0E0, ONE = 1.0E0, HALF = ONE/2.0E0, &
SPLIT1 = 0.425E0, SPLIT2 = 5.0E0, &
CONST1 = 0.180625E0, CONST2 = 1.6E0)
! COEFFICIENTS FOR P CLOSE TO 1/2
PARAMETER (A0 = 3.3871328727963666080E0, &
           A1 = 1.3314166789178437745E2, &
           A2 = 1.9715909503065514427E3, &
           A3 = 1.3731693765509461125E4, &
           A4 = 4.5921953931549871457E4, &
           A5 = 6.7265770927008700853E4, &
           A6 = 3.3430575583588128105E4, &
           A7 = 2.5090809287301226727E3, &
           B1 = 4.2313330701600911252E1, &
           B2 = 6.8718700749205790830E2, &
           B3 = 5.3941960214247511077E3, &
           B4 = 2.1213794301586595867E4, &
           B5 = 3.9307895800092710610E4, &
           B6 = 2.8729085735721942674E4, &
           B7 = 5.2264952788528545610E3)
! HASH SUM AB 55.88319 28806 14901 4439
!
! COEFFICIENTS FOR P NEITHER CLOSE TO 1/2 NOR 0 OR 1
PARAMETER (C0 = 1.42343711074968357734E0, &
C1 = 4.63033784615654529590E0, &
C2 = 5.76949722146069140550E0, &
C3 = 3.64784832476320460504E0, &
C4 = 1.27045825245236838258E0, &
C5 = 2.41780725177450611770E-1, &
C6 = 2.27238449892691845833E-2, &
C7 = 7.74545014278341407640E-4, &
D1 = 2.05319162663775882187E0, &
D2 = 1.67638483018380384940E0, &
D3 = 6.89767334985100004550E-1, &
D4 = 1.48103976427480074590E-1, &
D5 = 1.51986665636164571966E-2, &
D6 = 5.47593808499534494600E-4, &
D7 = 1.05075007164441684324E-9)
! HASH SUM CD 49.33206 50330 16102 89036
!
! COEFFICIENTS FOR P NEAR 0 OR 1
PARAMETER (E0 = 6.65790464350110377720E0, &
E1 = 5.46378491116411436990E0, &
E2 = 1.78482653991729133580E0, &
E3 = 2.96560571828504891230E-1, &
E4 = 2.65321895265761230930E-2, &
E5 = 1.24266094738807843860E-3, &
E6 = 2.71155556874348757815E-5, &
E7 = 2.01033439929228813265E-7, &
F1 = 5.99832206555887937690E-1, &
F2 = 1.36929880922735805310E-1, &
F3 = 1.48753612908506148525E-2, &
F4 = 7.86869131145613259100E-4, &
F5 = 1.84631831751005468180E-5, &
F6 = 1.42151175831644588870E-7, &
F7 = 2.04426310338993978564E-15)
! HASH SUM EF 47.52583 31754 92896 71629
!
!allocate(Q(np),R(np))

IFAULT = 0
Q = P - HALF
IF (ABS(Q) .LE. SPLIT1) THEN
  R = CONST1 - Q * Q
  X = Q * (((((((A7 * R + A6) * R + A5) * R + A4) * R + A3) &
        * R + A2) * R + A1) * R + A0) / (((((((B7 * R + B6) * R + B5) &
        * R + B4) * R + B3) * R + B2) * R + B1) * R + ONE)
  RETURN
ELSE
  IF (Q .LT. 0.0d0) THEN
    R= P
  ELSE
    R = ONE - P
  ENDIF
  IF (R .LE. ZERO) THEN
    IFAULT = 1
    X = ZERO
    RETURN
  ENDIF
  R = SQRT(-LOG(R))
  IF (R .LE. SPLIT2) THEN
    R = R - CONST2
    X = (((((((C7 * R + C6) * R + C5) * R + C4) * R &
        + C3) * R + C2) * R + C1) * R + C0) / (((((((D7 * R &
        + D6) * R + D5) * R + D4) * R + D3) * R + D2) * R   &
        + D1) * R + ONE)
  ELSE
    R = R - SPLIT2
    X = (((((((E7 * R + E6) * R + E5) * R + E4) * R &
      + E3) * R + E2) * R + E1) * R + E0) / (((((((F7 * R &
      + F6) * R + F5) * R + F4) * R + F3) * R + F2) * R &
      + F1) * R + ONE)
  ENDIF
  IF (Q .LT. 0.0d0) X = -X
  RETURN
ENDIF

!deallocate(Q,R)
END function InverseNormal

function normpdf(x,n) result(y)

use nrtype
implicit none
integer(i4b), intent(in) :: n 
real(dp),     intent(in) :: x(n)
real(dp)                 :: y(n)

y = exp( -0.5*x*x)/sqrt(2.0d0*pi_d)

end function normpdf

! X (d x d)      lower triangular matrix
! InvX (d x d)   inverse of X
function InvertLower(x,d) result(InvX)
  use nrtype
  implicit none
  integer(i4b), intent(in) :: d
  real(dp),     intent(in) :: x(d,d)
  real(dp)                 :: InvX(d,d)

  integer(i4b)             :: i1,i2

  InvX = 0.0d0

  !  x(1,1) * InvX(1,1) + x(1,2)*InvX(2,1) + x(1,3)*InvX(3,1) = 1.0
 
  !  x(2,1) * InvX(1,2) + x(2,2)*InvX(2,2) + x(2,3)*InvX(3,2) = 1.0
  !  x(2,1) * InvX(1,1) + x(2,2)*InvX(2,1) + x(2,3)*InvX(3,1) = 0.0


  !  x(3,1) * InvX(1,3) + x(3,2)*InvX(2,3) + x(3,3)*InvX(3,3) = 1.0
  !  x(3,1) * InvX(1,1) + x(3,2)*InvX(2,1) + x(3,3)*InvX(3,1) = 0.0
  !  x(3,1) * InvX(1,2) + x(3,2)*InvX(2,2) + x(3,3)*InvX(3,2) = 0.0


  !  x(1,1) * InvX(1,2) + x(1,2)*InvX(2,2) + x(1,3)*InvX(3,2) = 0.0
  !  x(1,1) * InvX(1,3) + x(1,2)*InvX(2,3) + x(1,3)*InvX(3,3) = 0.0
  !  x(2,1) * InvX(1,3) + x(2,2)*InvX(2,3) + x(2,3)*InvX(3,3) = 0.0

  do i1=1,d
    InvX(i1,i1) = 1.0d0 / x(i1,i1)
    do i2=1,i1-1
      InvX(i1,i2) = -dot_product(x(i1,1:i1-1),InvX(1:i1-1,i2)) / x(i1,i1)
    end do
  end do
end function InvertLower

!-------------------------------------------------------------------------
! Subroutines
!-------------------------------------------------------------------------
subroutine bsolve1L(x,A,n,b)
! solve Ax=b using backsubstitution
! A is (n x n), lower diagonal matrix 
  use nrtype
  implicit none
  integer(i4b), intent(in) :: n
  real(dp), intent(in) :: A(n,n),b(n)
  real(dp), intent(out) :: x(n)
  integer(i4b) i1,i2
  real(dp) sum
	
  do i1=1,n
    sum=b(i1)
	do i2=i1-1,1,-1
	  sum=sum-A(i1,i2)*x(i2)
	end do
	x(i1)=sum/A(i1,i1)
  end do
end subroutine
	   
subroutine bsolve1U(x,A,n,b)
! solve Ax=b using backsubstitution
! A is (n x n), upper diagonal matrix
  use nrtype
  implicit none
  integer(i4b), intent(in) :: n
  real(dp), intent(in) :: A(n,n),b(n) 
  real(dp), intent(out) :: x(n)
  integer(i4b) i1,i2
  real(dp) sum
 
  do i1=n,1,-1
    sum=b(i1)
	do i2=n,i1+1,-1
	  sum=sum-A(i2,i1)*x(i2)
	end do
	x(i1)=sum/A(i1,i1)
  end do
end subroutine 

subroutine bsolve2L(x,A,k,b)
! solve A*x=b using backsubstitution
! A is a vector containing terms of a lower diagonal matrix of
! length k (in row order)
  use nrtype
  implicit none
  integer(i4b), intent(in) :: k
  real(dp), intent(in) :: A(k*(k+1)/2),b(k)
  real(dp), intent(out) :: x(k)
  integer(i4b) i1,i2,i3
  real(dp) sum

  do i1=1,k
    sum=b(i1)
	do i2=1,i1-1
	  i3=i1*(i1-1)/2+i2
	  sum=sum-A(i3)*x(i2)
	end do
	i3=i1*(i1-1)/2+i1
	x(i1)=sum/A(i3)
  end do
end subroutine

subroutine bsolve2U(x,A,k,b)
! solve A*x=b using backsubstitution
! A is a vector containing terms of an upper diagonal matrix of
! length k (in column order)
  use nrtype
  implicit none
  integer(i4b), intent(in) :: k
  real(dp), intent(in) :: A(k*(k+1)/2),b(k)
  real(dp), intent(out) :: x(k)
  integer(i4b) i1,i2,i3
  real(dp) sum
    
  do i1=k,1,-1
    sum=b(i1)
	do i2=k,i1+1,-1
	  i3=i2*(i2-1)/2+i1
	  sum=sum-a(i3)*x(i2)
	end do
    i3=i1*(i1-1)/2+i1
	x(i1)=sum/A(i3)
  end do
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! subroutine chol1(A,n,L)
! Compute the Cholesky decomposition of matrix A
! A(n,n)	real(dp): positive definite symmetric matrix
! n		integer: dimension of A
! L(n,n)  real(dp): satisfies A=L*L'
subroutine chol1(L,A,n)
  use nrtype
  implicit none
  integer(i4b), intent(in) :: n   
  real(dp), intent(in) :: A(n,n)
  real(dp), intent(out) :: L(n,n)
  real(dp) sum
  integer(i4b) i,j,k
  do i=1,n
    do j=i,n
      sum=A(i,j)
      do k=i-1,1,-1
	    sum=sum-L(i,k)*L(j,k)
	  end do
	  if (i==j) then
	    if (sum<0.0) then
	      print *, 'Sum < 0 in chol1'
	   !   pause
   	    end if
		L(i,i)=sqrt(sum)
	  else
	    L(j,i)=sum/L(i,i)
	    L(i,j)=0.0
	  end if
	end do ! j=i,n
  end do   ! i=1,n
end subroutine chol1

subroutine cholinv1(Linv,L,n)
! Given L, an n x n lower diagonal matrix,
! compute L inverse
  use nrtype
  implicit none
  integer(i4b), intent(in) :: n
  real(dp), intent(in) :: L(n,n)
  real(dp), intent(out) :: Linv(n,n)
  integer(i4b) i1,i2
  real(dp) eye(n)

  eye=0.0
  do i1=1,n
    eye(i1)=1.0
	call bsolve1L(Linv(:,i1),L,n,eye)
    eye(i1)=0.0
  end do
end subroutine cholinv1

subroutine cholinv2(Linv,L,n)
! Given L, an n*(n+1)/2 matrix containing the lower diagonal matrix L
! compute L inverse
! L(i,j) = L(i*(i-1)/2+j) for i>=j
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: L(n*(n+1)/2)
  real(8), intent(out) :: Linv(n*(n+1)/2)
  real(8) sum
  integer i1,i2,i3,i4,i5

  do i1=1,n
    do i2=1,i1
	  i3=i1*(i1-1)/2+i2
	  if (i1==i2) then
	    Linv(i3)=1.0/L(i3)
	  else
	    sum=0.0
	    do i4=0,i1-i2-1
	      i5=(i2+i4)*(i2+i4-1)/2+i2
	   	  sum=sum-Linv(i5)*L(i3+i4)
		end do
		Linv(i3)=sum/L(i3+i4)
	  end if
    end do
  end do
end subroutine cholinv2
	
subroutine invert1(Kinv,det,K,n)
! Given K, an (n x n) positive definite symmetric matrix,
! invert it using cholesky decomp. Also, returns the
! cholesky decomposition of K and its inverse as well as the
! the determinant of Kinv.
!	K	 (n x n) real symmetric pos.def. matrix
!   n    (1 x 1) integer(i4b) dimension of K
!   Kinv (n x n) real inverse of K 
!   L	 (n x n) real cholesky decomposition of K. K = L*L'
!   Linv (n x n) real inverse of L.  L*Linv=I
!   det  (1 x 1) real determinant of Kinv
  use nrtype
  implicit none
  integer(i4b), intent(in) :: n
  real(dp), intent(in) :: K(n,n)
  real(dp), intent(out) :: Kinv(n,n),det
  real(dp) L(n,n),Linv(n,n)
  integer(i4b) i1,i2,i3

  call chol1(L,K,n)
  call cholinv1(Linv,L,n)

  Kinv=0.0
  det=1.0
  do i1=1,n
    det=det*Linv(i1,i1)
    do i2=1,i1
	  do i3=1,n
	    Kinv(i1,i2)=Kinv(i1,i2)+Linv(i3,i1)*Linv(i3,i2)
	  end do
	  if (i2<i1) Kinv(i2,i1)=Kinv(i1,i2)
	end do
  end do
  det=det*det
end subroutine  

subroutine invert2(Kinv,det,K,n)
! Given K, the lower diagonal of an (n x n) positive definite
! symmetric matrix, invert it using cholesky decomp. Also, 
! calculate the determinant of Kinv.
  implicit none
  integer, intent(in) :: n
  real(8), intent(in) :: K(n*(n+1)/2)
  real(8), intent(out) :: Kinv(n*(n+1)/2),det
  real(8) L(n*(n+1)/2), Linv(n*(n+1)/2)
  integer i1,i2,i3,i11,i21,i31,i32

  call chol2(L,K,n)
  call cholinv2(Linv,L,n)

  Kinv=0.0
  det=1.0
  do i1=1,n
    i11=i1*(i1-1)/2+i1
    det=det*Linv(i11)
    do i2=1,i1
	  do i3=n,i1,-1
	    i21=i1*(i1-1)/2+i2
		i31=i3*(i3-1)/2+i1
		i32=i3*(i3-1)/2+i2
	    Kinv(i21)=Kinv(i21)+Linv(i31)*Linv(i32)
	  end do
	end do
  end do
  det=det*det
end subroutine  

subroutine chol2(L,A,n)
! subroutine chol2(L,A,n)
! Compute the Cholesky decomposition of matrix A, where the lower 
! diagonal of A is stored as a long vector n*(n+1)/2 and 
! A(i,j) = A(i*(i-1)/2+j) for i>=j
!
! A(n*(n+1)/2) real(8): positive definite symmetric matrix 
! n		     integer: dimension of A
! L(n*(n+1)/2) real(8): satisfies A=L*L'
! L(i,j) =  L(i*(i-1)/2+j
  implicit none
  integer, intent(in) :: n   
  real(8), intent(in) :: A(n*(n+1)/2)
  real(8), intent(out) :: L(n*(n+1)/2)
  real(8) sum,p
  integer i,j,k,ij,ik,jk,ji
  do i=1,n
    do j=i,n
	  ij=j*(j-1)/2+i
	  sum=a(ij)
	  do k=i-1,1,-1
	    ik=i*(i-1)/2+k
		jk=j*(j-1)/2+k
		sum=sum-l(ik)*l(jk)
	  end do
	  if (i==j) then
	    p=sqrt(sum)
		L(ij)=p
	  else
	    L(ij)=sum/p
	  end if
	end do ! j=i,n
  end do   ! i=1,n
end subroutine chol2

subroutine gaucheb(x,w)
  use nrtype
  implicit none
  ! Compute nodes and weights for the n-point Gauss-Chebyshev integration
  ! formula
  ! x	= real(dp) (n x 1)	vector of points
  ! w	= real(dp) (n x 1)  vector of weights
  real(dp), intent(out) :: x(:),w(:)
  integer(i4b) i1,n

  n=size(x)
  x=cos(pi_d*(dble((/(i1,i1=n,1,-1)/))-0.5d0)/dble(n))
  w=pi_d/dble(n)
end subroutine gaucheb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Calculate nodes and weights for Gauss-Legendre integration on [-1,1]
!
! n	= integer (1 x 1)    number of nodes to use 
! x	= real(8) (n x 1)    nodes
! w	= real(8) (n x 1)    weights
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine gauleg(x,w)
  use nrtype
  
  implicit none
  real(dp), intent(out) :: x(:),w(:)
  real(dp) tol,xm,x1,z,z1,p1,p2,p3,pp
  integer(i4b) m,i,j,n
 
  n=size(x) 
  tol=3.0d-11
  x=0.0d0
  w=0.0d0

  m=(n+1)/2
  xm=0.0d0
  x1=1.0d0
  do i=1,m
    z=cos(pi*(dble(i)-0.25d0)/(dble(n)+0.5d0))
    z1=z+1.0d0
    do 
      if (abs(z-z1)<=tol) then
        exit
      end if
      p1=1.0d0
      p2=0.0d0
      do j=1,n
        p3=p2
        p2=p1
        p1=((2.0d0*dble(j)-1.0d0)*z*p2-(dble(j)-1.0d0)*p3)/dble(j)
      end do
      pp=dble(n)*(z*p1-p2)/(z*z-1.0d0)
      z1=z
      z=z1-p1/pp
    end do
    x(i)=xm-x1*z
    x(n+1-i)=xm+x1*z
    w(i)=2.0d0*x1/((1.0d0-z*z)*pp*pp)
    w(n+1-i)=w(i)
  end do
end subroutine gauleg

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! gaulag(x,w,alpha)
!
! Create nodes and weights for Gauss-Laguerre integration
! That is, nodes and weights for integration against x**alpha * exp(-x) on
! the domain (0,inf)
!
! x     =   real(dp) (n x 1) nodes
! w	    =   real(dp) (n x 1) weights
! alpha0 =   real(dp) (1 x 1) optional parameter
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine gaulag(x,w,alpha0)
  use nrtype
  implicit none
  real(dp), intent(out) :: x(:),w(:)
  real(dp), optional :: alpha0
  real(16), parameter :: epp=3.0d-14
  integer(i4b), parameter :: maxit=10
  integer(i4b) i,j,j1,n
  real(16) z,z1,ai,p1,p2,p3,pp,alpha,w1

  if (present(alpha0)) then
    alpha=real(alpha0,16)
  else
    alpha=real(0.0,16)
  end if

  n=size(x)
  x=0.0d0
  w=0.0d0
  
  do i=1,n
    if (i==1) then
      z=(real(1.0,16)+alpha)*(real(3.0,16)+real(0.92,16)*alpha)/   &
        (real(1.0,16)+real(2.4,16)*real(n,16)+real(1.8,16)*alpha)
    elseif (i==2) then
      z=z+(real(15.0,16)+real(6.25,16)*alpha)/(real(1.0,16)+real(0.9,16)*alpha+real(2.5,16)*real(n,16))
    else
      ai=real(i,16)-real(2.0,16)
      z=z+((real(1.0,16)+real(2.55,16)*ai)/ &
           (real(1.9,16)*ai)+real(1.26,16)*ai*alpha/(real(1.0,16)+real(3.5,16)*ai))*(z-x(i-2))/(real(1.0,16)+real(0.3,16)*alpha)
    end if
    do j=1,maxit
      p1=real(1.0,16)
      p2=real(0.0,16)
      do j1=1,n
        p3=p2
        p2=p1
        p1=((real(2.0,16)*real(j1,16)-real(1.0,16)+alpha-z)*p2-(real(j1,16)-real(1.0,16)+alpha)*p3)/real(j1,16)
      end do
      pp=(real(n,16)*p1-(real(n,16)+alpha)*p2)/z
      z1=z
      z=z1-p1/pp
      if (abs(z-z1)<epp) then
        exit
      end if
    end do
    if (j>maxit) then
      print *,z-z1,j
      stop 'Too many iterations in gaulag'
    end if
    x(i)=real(z,dp)
    w(i)=-exp(gammln(real(alpha,dp)+real(n,dp))-gammln(real(n,dp)))/(real(pp,dp)*real(n,dp)*real(p2,dp))
    w(i)=real(w1,dp)
  end do
end subroutine gaulag


! Compute numerical gradient of func

subroutine numdif00(func,nx,x0,df,iuser,user)
  use nrtype
  implicit none

  integer(i4b), intent(in) :: nx
  real(dp), intent(in) :: x0(nx)
  real(dp), intent(out) :: df(nx)
  integer(i4b), intent(in) :: iuser(*)
  real(dp), intent(in) :: user(*)
  external func
  !interface
  !  subroutine func(mode,nx,x,f,df,nstate,iuser,user)
  !    use nrtype
  !!                                         ! 1 = value of gradient only
  !                                         ! 2 = values of objective and gradient
  !    integer(i4b), intent(in) :: nx       ! length of x
  !    real(dp), intent(in) :: x(nx)        ! input x
  !    real(dp), intent(out) :: f           ! value of objective function
  !    real(dp), intent(out) :: df(nx)      ! value of gradient
  !	  integer(i4b), intent(in) :: nstate   ! 1 if first evaluation, 0 otherwise
  !    integer(i4b), intent(in) :: iuser(*) ! integer(i4b) inputs
  !    real(dp), intent(in) :: user(*)      ! real(dp) inputs
  !  end subroutine
  !end interface

  real(dp) x1(nx),x2(nx),dx,f0,f1,f2,df0(nx)
  integer(i4b) i1,mode,nstate

  mode=0
  nstate=1
  
  call func(mode,nx,x0,f0,df0,nstate,iuser,user)

  df=0.d0

  do i1=1,nx
    x1=x0
    x2=x0
    dx=1.d-4*max(0.1d0,abs(x0(i1)))
    x1(i1)=x0(i1)+dx
    x2(i1)=2.d0*x0(i1)-x1(i1)
    call func(mode,nx,x1,f1,df0,nstate,iuser,user)
    call func(mode,nx,x2,f2,df0,nstate,iuser,user)
    df(i1)=(f1-f2)/(x1(i1)-x2(i1))
  enddo
end subroutine numdif00

SUBROUTINE CONFUN(MODE,NCNLN,N,LDCJ,NEEDC,X,C,CJAC,NSTATE,IUSER,USER)
USE NRTYPE

IMPLICIT NONE

INTEGER(I4B), INTENT(INOUT) :: MODE
INTEGER(I4B), INTENT(IN) :: NCNLN,N,LDCJ,NEEDC(NCNLN),NSTATE,IUSER(*)
REAL(DP), INTENT(IN) :: X(N)
REAL(DP), INTENT(OUT) :: C(NCNLN),CJAC(LDCJ,N)
REAL(DP), INTENT(IN) :: USER(*)
c=0.d0
cjac=0.d0
END SUBROUTINE confun

subroutine kron1(a,b,c)
  ! subroutine kron1(a,b,c)
  ! Kronecker product of a and b where a is size (na x 1) and b is size(nb x 1)
  ! a	= real(dp) (na x 1)
  ! b	= real(dp) (nb x 1)
  ! c	= real(dp) (na*nb x 1)
  use nrtype
  implicit none
  real(dp), intent(in) :: a(:),b(:)
  real(dp), intent(out) :: c(:)
  
  real(dp), allocatable :: atemp1(:,:),atemp2(:),btemp1(:,:),btemp2(:)
  integer(i4b) na,nb

  na=size(a)
  nb=size(b)

  allocate(atemp1(nb,na),atemp2(na*nb))
  allocate(btemp1(nb,na),btemp2(na*nb))

  atemp1=spread(a,1,nb)
  atemp2=reshape(atemp1,(/na*nb/))
  btemp1=spread(b,2,na)
  btemp2=reshape(btemp1,(/na*nb/))

  c=atemp2*btemp2
  deallocate(atemp1,atemp2,btemp1,btemp2)
end subroutine kron1

subroutine ikron1(a,b,c)
  ! subroutine kron1(a,b,c)
  ! Kronecker product of a and b where a is size (na x 1) and b is size(nb x 1)
  ! a	= integer(i4b) (na x 1)
  ! b	= integer(i4b) (nb x 1)
  ! c	= integer(i4b) (na*nb x 1)
  use nrtype
  implicit none
  integer(i4b), intent(in) :: a(:),b(:)
  integer(i4b), intent(out) :: c(:)
  
  integer(i4b), allocatable :: atemp1(:,:),atemp2(:),btemp1(:,:),btemp2(:)
  integer(i4b) na,nb

  na=size(a)
  nb=size(b)

  allocate(atemp1(nb,na),atemp2(na*nb))
  allocate(btemp1(nb,na),btemp2(na*nb))

  atemp1=spread(a,1,nb)
  atemp2=reshape(atemp1,(/na*nb/))
  btemp1=spread(b,2,na)
  btemp2=reshape(btemp1,(/na*nb/))

  c=atemp2*btemp2
  deallocate(atemp1,atemp2,btemp1,btemp2)
end subroutine ikron1


subroutine linspace(a,b,n,x)
  use nrtype
  implicit none
  real(dp), intent(in) :: a,b
  integer(i4b), intent(in) :: n
  real(dp), intent(out) :: x(n)
 real(dp) step
 integer(i4b) i1

  step=(b-a)/real(n-1,dp)
  x(1)=a
  x(n)=b

  do i1=2,n-1
   x(i1)=x(i1-1)+step
 enddo
!  x(2:n-1)=a+(b-a)*real((/1:n-2/),dp)/real(n-1,dp)
!  x(2:n-1)=a+(b-a)*real((/1,n-2/),dp)/real(n-1,dp)
end subroutine linspace

! Compute c = kron(a,b)
subroutine kron(a,b,c)

use nrtype
implicit none
real(dp), intent(in) :: a(:,:)
real(dp), intent(in) :: b(:,:)
real(dp), intent(out) :: c(:,:)

integer(i4b) :: na1,na2,nb1,nb2
real(dp), allocatable :: btemp(:,:)

na1=size(a,1)
na2=size(a,2)
nb1=size(b,1)
nb2=size(b,2)

allocate(btemp(na1*nb1,na2))

btemp = reshape(spread(a,1,nb1),(/na1*nb1,na2/))
c = reshape(spread(btemp,2,nb2),(/na1*nb1,na2*nb2/))

deallocate(btemp)

allocate(btemp(nb1,na2*nb2))
btemp = reshape(spread(b,3,na2),(/nb1,na2*nb2/))

c = c * &
    reshape(spread(btemp,2,na1),(/na1*nb1,na2*nb2/))
deallocate(btemp)

end subroutine kron

! Compute c = ikron(a,b)
! kronecker product of two integer-valued matrixes (a,b)
subroutine ikron(a,b,c)

use nrtype
implicit none
integer(i4b), intent(in)  :: a(:,:)
integer(i4b), intent(in)  :: b(:,:)
integer(i4b), intent(out) :: c(:,:)

integer(i4b)              :: na1,na2,nb1,nb2
integer(i4b), allocatable :: btemp(:,:)

na1=size(a,1)
na2=size(a,2)
nb1=size(b,1)
nb2=size(b,2)

allocate(btemp(na1*nb1,na2))

btemp = reshape(spread(a,1,nb1),(/na1*nb1,na2/))
c = reshape(spread(btemp,2,nb2),(/na1*nb1,na2*nb2/))

deallocate(btemp)

allocate(btemp(nb1,na2*nb2))
btemp = reshape(spread(b,3,na2),(/nb1,na2*nb2/))

c = c * &
    reshape(spread(btemp,2,na1),(/na1*nb1,na2*nb2/))
deallocate(btemp)

end subroutine ikron


end module ToolsModule


