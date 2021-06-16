! ToolsModule
!
! Functions
!-----------------------------------------------------------------------------
! function normcdf_mkl(x,nx) result(y)    Compute normal CDF using Intel MKL
! function InverseNormal_mkl(p,n,ifault) result(x)   Compute inverse CDF normal
!                                                    using Intel MKL
!
!-----------------------------------------------------------------------------
! Subroutines
!-----------------------------------------------------------------------------
! subroutine kron1(a,b,c)
! subroutine linspace(a,b,n,x)
!-----------------------------------------------------------------------------
module ToolsModule

implicit none

contains


function normcdf_mkl(x,nx) result(y)
  use ConstantsModule
  implicit none
  integer(i4b), intent(in) :: nx
  real(dp),     intent(in) :: x(nx)
  real(dp)                 :: y(nx)

  include "mkl_vml.fi"

  ! normal cdf from intel_MKL
  call vdcdfnorm(nx,x,y)

end function normcdf_mkl


function InverseNormal_mkl(p,n,ifault) result(x)
  use ConstantsModule
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

!-------------------------------------------------------------------------
! Subroutines
!-------------------------------------------------------------------------

subroutine kron1(a,b,c)
  ! subroutine kron1(a,b,c)
  ! Kronecker product of a and b where a is size (na x 1) and b is size(nb x 1)
  ! a	= real(dp) (na x 1)
  ! b	= real(dp) (nb x 1)
  ! c	= real(dp) (na*nb x 1)
  use ConstantsModule
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
  use ConstantsModule
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

