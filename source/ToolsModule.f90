! ToolsModule
!
! Functions
!-----------------------------------------------------------------------------
! function normcdf_mkl(x,nx) result(y)    Compute normal CDF using Intel MKL
! function InverseNormal_mkl(p,n,ifault) result(x)   Compute inverse CDF normal
!                                                    using Intel MKL
! function InvertLower(x,d) result(InvX)  Invert lower triangular matrix.
!
!-----------------------------------------------------------------------------
! Subroutines
!-----------------------------------------------------------------------------
! subroutine gaucheb(x,w)		Nodes and weights for Gauss-Chebyshev integration
! subroutine gauleg(x,w,n)		Nodes and weights for Gauss-Legendre integration
! subroutine kron1(a,b,c)
! subroutine linspace(a,b,n,x)
!-----------------------------------------------------------------------------
module ToolsModule

implicit none

contains


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
  
end function InverseNormal_mkl


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
end subroutine linspace

end module ToolsModule

