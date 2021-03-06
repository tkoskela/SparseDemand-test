!
! Modification history
! 2018JUL26 LN  convert (1:k) to (/(ix,ix=1,k)/).  i.e. convert to gfortran
!               syntax
! 2014DEC22 LN  add ComputeKnotAverages to compute a grid of points defined as the averages of the knots
! 30jun2011 LN  change specification statement to use allocatable variables instead of
!               assumed shape
! 10jan09 LPN Fix bug in bspline3 and clean up comments in bspline2 and bspline3.
! 13jun08 lpn Delete obsolete code and edit bspline3 to output nonzero indexes
! 22apr08 lpn Edit bspline3.
! 18apr08 Finish editing of sparse-capable versions of bspline1 and bspline2.
! 17apr08 Begin edit to take advantage of sparsity.
! 19sep07 Fixed bug in calling sequence when sorted=1.
! 13sep06 Accuracy test passed. option sorted=1 not tested.
! Contains:
!
!  Subroutines that compute (b,db,ddb) where
!  b is the B-Spline collocation matrix of degree k, 
!  db is its gradient
!  ddb is its hessian.
!
!   1) subroutine bspline1(x,knots,k,sorted,b,nonzero,loknots,hiknots)
!   2) subroutine bspline2(x,knots,k,sorted,b,db,nonzero,loknots,hiknots)
!   3) subroutine bspline3(x,knots,k,sorted,b,db,ddb,nonzero,loknots,hiknots)
!   4) subroutine bspline4(b,x,xg,k,i,m,xg0,xg1)
!      Compute i'th B-spline of degree k (order k+1) at 
!      the vector of points x. If k==3, the splines are 
!      cubic splines. The values of x lie in the
!      interval [xg(m),xg(m+1)].
!   5) subroutine ComputeKnotAverages(knots,k,x)


module SplineTools
    
contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! subroutine bspline1(x,knots,k,sorted,b,nonzero,loknots,hiknots)
! 
! Compute B-splines of degree k (order k+1) at the vector of points x.
! 
! Required inputs:
!
!   real(8) x(nx)                  Vector of points.
!   real(8) knots(ngrid)           Knots for spline.
!   integer k			   Degree of B-spline.
!			           If k==3, the splines are cubic splines.
!   integer sorted		   1 if x is sorted, 0 otherwise.
!
! Required outputs
!
!   real(8) b(ngrid-1+k,nx)        b(i,j) =  B^k_i(x(j)).
!   integer(i4b) nonzero(nx)       nonzero(i) = first index of nonzero elements 
!                                  in column i
!
! Optional inputs:
!
!   real(8) loknots(k)             Extension of grid below knots(1).
!   real(8) hiknots(k)             Extension of grid above knots(n).
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine bspline1(x,knots,k,sorted,b,nonzero,loknots,hiknots)
  use nrtype
  use nr, only: locate,indexx
  implicit none
  real(dp), intent(in)           :: x(:)
  real(dp), intent(in)           :: knots(:)
  integer(i4b), intent(in)       :: k
  integer(i4b), intent(in)       :: sorted
  real(dp), intent(out)          :: b(:,:)
  integer(i4b), intent(out)      :: nonzero(:)
  real(dp), optional, intent(in) :: loknots(:),hiknots(:)

  integer(i4b)              :: nx,nknots,n1,n2,jhi,j1,j2,ix
  integer(i4b), allocatable :: i00(:),i01(:),indx(:)
  logical, allocatable      :: temp(:,:)
  real(dp), allocatable     :: allknots(:)
  real(dp), allocatable     :: t1d(:,:),t2d(:,:),t1n(:,:),t2n(:,:),t1(:,:),t2(:,:)
  real(dp), allocatable     :: temp2(:,:),temp3(:,:)
  real(dp), allocatable     :: x1(:)

  nx=size(x)
  nknots=size(knots)
  allocate(allknots(nknots+2*k))
  allocate(x1(nx))
  x1=x

  if (sorted==0) then
    allocate(i00(nx),i01(nx))
 
    ! Sort x in ascending order.
    call indexx(x,i00)		
    x1=x(i00)

    ! Indices for unsorting x.
    call indexx(i00,i01)		
  endif

  if (present(loknots) .and. present(hiknots)) then
    allknots((/(ix,ix=1,k)/))=loknots
    allknots((/(ix,ix=k+1,nknots+k)/))=knots
    allknots((/(ix,ix=nknots+k+1,nknots+2*k)/))=hiknots
  elseif (present(loknots) .and. .not.present(hiknots)) then
    allknots((/(ix,ix=1,k)/))=loknots
    allknots((/(ix,ix=k+1,nknots+k)/))=knots
    allknots((/(ix,ix=nknots+k+1,nknots+2*k)/))=knots(nknots)  
  elseif (.not.present(loknots) .and. present(hiknots)) then
    allknots((/(ix,ix=1,k)/))=knots(1)
    allknots((/(ix,ix=k+1,nknots+k)/))=knots
    allknots((/(ix,ix=nknots+k+1,nknots+2*k)/))=hiknots
  else
    allknots((/(ix,ix=1,k)/))=knots(1)
    allknots((/(ix,ix=k+1,nknots+k)/))=knots
    allknots((/(ix,ix=nknots+k+1,nknots+2*k)/))=knots(nknots)
  endif

  n2=nknots+k-1
  b=0.d0
  nonzero=0

  jhi=1

  do j1=k+2,nknots+k
    if (jhi>nx) exit
    n1=locate(x1((/(ix,ix=jhi,nx)/)),allknots(j1))
    nonzero((/(ix,ix=jhi,nx)/))=nonzero((/(ix,ix=jhi,nx)/))+1
    
    if (n1>0) then
      allocate(indx(n1),t1d(n2,n1),t2d(n2,n1),t1n(n2,n1),t2n(n2,n1),t1(n2,n1),t2(n2,n1))
      allocate(temp(n2,n1),temp2(n2,n1),temp3(n2,n1))

      indx=(/(ix,ix=jhi,jhi+n1-1)/)
      b(j1-1,indx)=1.d0
      do j2=1,k
        t1d=spread(allknots((/(ix,ix=1+j2,j2+n2)/))-allknots((/(ix,ix=1,n2)/)),dim=2,ncopies=n1)
        t2d=spread(allknots((/(ix,ix=2+j2,j2+n2+1)/))-allknots((/(ix,ix=2,n2+1)/)),2,n1)
        t1n=transpose(spread(x1(indx),dim=2,ncopies=n2))-spread(allknots((/(ix,ix=1,n2)/)),dim=2,ncopies=n1)
        t2n=spread(allknots((/(ix,ix=2+j2,j2+n2+1)/)),dim=2,ncopies=n1)-transpose(spread(x1(indx),2,n2))

        t1=0.d0
        t2=t1
      
        ! If t1d~=0, t1=t1n/t1d
        ! If t1d==0, t1=0.d0
        temp=(t1d==0.d0)			
        t1d=merge(1.d0,t1d,temp)      
        t1=t1n/t1d
        t1=merge(0.d0,t1,temp)

        ! If t2d~=0, t2=t2n/t2d
        ! If t2d==0, t2=0.d0
        temp=(t2d==0.d0)
        t2d=merge(1.d0,t2d,temp)
        t2=t2n/t2d					
        t2=merge(0.d0,t2,temp)	

        temp2=0.d0
        temp3=0.d0
        temp2((/(ix,ix=1,n2-1)/),:)=t2((/(ix,ix=1,n2-1)/),:)
        temp3((/(ix,ix=1,n2-1)/),:)=b((/(ix,ix=2,n2)/),indx)
        b(:,indx)=t1*b(:,indx)+temp2*temp3
      enddo
      jhi=jhi+n1
      deallocate(indx,t1d,t2d,t1n,t2n,t1,t2,temp,temp2,temp3)
    endif
  enddo

  ! If some values of x, still not processed, then compute b for these last x's.
  if (jhi<=nx) then
    n1=nx+1-jhi
    nonzero((/(ix,ix=jhi,nx)/))=nonzero((/(ix,ix=jhi,nx)/))+1
    allocate(indx(n1),t1d(n2,n1),t2d(n2,n1),t1n(n2,n1),t2n(n2,n1),t1(n2,n1),t2(n2,n1))
    allocate(temp(n2,n1),temp2(n2,n1),temp3(n2,n1))
    
    indx=(/(ix,ix=jhi,jhi+n1-1)/)

    b(nknots+k-1,indx)=1.d0
    
    do j2=1,k
      t1d=spread(allknots((/(ix,ix=1+j2,j2+n2)/))-allknots((/(ix,ix=1,n2)/)),dim=2,ncopies=n1)
      t2d=spread(allknots((/(ix,ix=2+j2,j2+n2+1)/))-allknots((/(ix,ix=2,n2+1)/)),2,n1)
      t1n=transpose(spread(x1(indx),dim=2,ncopies=n2))-spread(allknots((/(ix,ix=1,n2)/)),dim=2,ncopies=n1)
      t2n=spread(allknots((/(ix,ix=2+j2,j2+n2+1)/)),dim=2,ncopies=n1)-transpose(spread(x1(indx),2,n2))

      t1=0.d0
      t2=t1
      
      ! If t1d~=0, t1=t1n/t1d
      ! If t1d==0, t1=0.d0
      temp=(t1d==0.d0)				
      t1d=merge(1.d0,t1d,temp)      
      t1=t1n/t1d
      t1=merge(0.d0,t1,temp)

      ! If t2d~=0, t2=t2n/t2d
      ! If t2d==0, t2=0.d0
      temp=(t2d==0.d0)
      t2d=merge(1.d0,t2d,temp)
      t2=t2n/t2d					
      t2=merge(0.d0,t2,temp)	
	
      temp2=0.d0
      temp3=0.d0
      temp2((/(ix,ix=1,n2-1)/),:)=t2((/(ix,ix=1,n2-1)/),:)
      temp3((/(ix,ix=1,n2-1)/),:)=b((/(ix,ix=2,n2)/),indx)
      b(:,indx)=t1*b(:,indx)+temp2*temp3
    enddo
    deallocate(indx,t1d,t2d,t1n,t2n,t1,t2,temp,temp2,temp3)
  endif

  deallocate(x1)
  deallocate(allknots)

  if (sorted==0) then
    b=b(:,i01)
    nonzero=nonzero(i01)
    deallocate(i00,i01)
  endif
end subroutine bspline1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! subroutine bspline2(x,knots,k,sorted,b,db,nonzero,loknots,hiknots)
! 
! Compute B-splines of degree k (order k+1) at the vector of points x.
! Also compute 1st derivative
! 
! Required inputs:
!
!   real(dp) x(nx)                 Vector of points.
!   real(dp) knots(ngrid)          Knots for spline.
!   integer k		 	   Degree of B-spline.
!				   If k==3, the splines are cubic splines.
!   integer sorted		   1 if x is sorted, 0 otherwise.
!
! Required outputs
!
!   real(dp) b(ngrid-1+k,nx)       b(i,j) =  B^k_i(x(j)).
!   real(dp) db(ngrid-1+k,nx)	   db(i,j) = dB^k_i(x(j))/dx.
!   integer(i4b) nonzero(nx)       nonzero(i) = index of first nonzero element
!                                  in column i
!
! Optional inputs:
!
!   real(dp) loknots(k)             Extension of grid below knots(1).
!   real(dp) hiknots(k)             Extension of grid above knots(n).
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine bspline2(x,knots,k,sorted,b,db,nonzero,loknots,hiknots)
  use nrtype
  use nr, only: locate,indexx
  implicit none
  real(dp),     intent(in)          :: x(:)
  real(dp),     intent(in)          :: knots(:)
  integer(i4b), intent(in)          :: k
  integer(i4b), intent(in)          :: sorted
  real(dp),     intent(out)         :: b(:,:)
  real(dp),     intent(out)         :: db(:,:)
  integer(i4b), intent(out)         :: nonzero(:)
  real(dp),     intent(in),optional :: loknots(:),hiknots(:)

  integer(i4b)              :: nx,nknots,n1,n2,jhi,j1,j2,ix
  integer(i4b), allocatable :: i00(:),i01(:),indx(:)
  logical, allocatable      :: temp(:,:)
  real(dp), allocatable     :: t1d(:,:),t2d(:,:),t1n(:,:),t2n(:,:),t1(:,:),t2(:,:)
  real(dp), allocatable     :: t4(:,:),t5(:,:),temp2(:,:),temp3(:,:)
  real(dp), allocatable     :: x1(:)
  real(dp), allocatable     :: allknots(:)

  nx=size(x)
  nknots=size(knots)
  allocate(x1(nx))
  allocate(allknots(nknots+2*k))

  x1=x

  if (sorted==0) then 
    allocate(i00(nx),i01(nx))
    i00=0

    ! Sort x in ascending order.
    call indexx(x,i00)		
    x1=x(i00)

    ! Indices for unsorting x.
    call indexx(i00,i01)		
  endif

  if (present(loknots) .and. present(hiknots)) then
    allknots((/(ix,ix=1,k)/))=loknots
    allknots((/(ix,ix=k+1,nknots+k)/))=knots
    allknots((/(ix,ix=nknots+k+1,nknots+2*k)/))=hiknots
  elseif (present(loknots) .and. .not.present(hiknots)) then
    allknots((/(ix,ix=1,k)/))=loknots
    allknots((/(ix,ix=k+1,nknots+k)/))=knots
    allknots((/(ix,ix=nknots+k+1,nknots+2*k)/))=knots(nknots)  
  elseif (.not.present(loknots) .and. present(hiknots)) then
    allknots((/(ix,ix=1,k)/))=knots(1)
    allknots((/(ix,ix=k+1,nknots+k)/))=knots
    allknots((/(ix,ix=nknots+k+1,nknots+2*k)/))=hiknots
  else
    allknots((/(ix,ix=1,k)/))=knots(1)
    allknots((/(ix,ix=k+1,nknots+k)/))=knots
    allknots((/(ix,ix=nknots+k+1,nknots+2*k)/))=knots(nknots)
  endif

  n2=nknots+k-1
  b=0.d0
  db=b
  nonzero=0

  jhi=1

  do j1=k+2,nknots+k
 
    if (jhi>nx) exit
 
    n1=locate(x1((/(ix,ix=jhi,nx)/)),allknots(j1))
    nonzero((/(ix,ix=jhi,nx)/))=nonzero((/(ix,ix=jhi,nx)/))+1
    if (n1>0) then
      allocate(indx(n1),t1d(n2,n1),t2d(n2,n1),t1n(n2,n1),t2n(n2,n1),t1(n2,n1),t2(n2,n1))
      allocate(t4(n2,n1),t5(n2,n1),temp(n2,n1),temp2(n2,n1),temp3(n2,n1))
      indx=(/(ix,ix=jhi,jhi+n1-1)/)
      b(j1-1,indx)=1.d0
      do j2=1,k
        t1d=spread(allknots((/(ix,ix=1+j2,j2+n2)/))-allknots((/(ix,ix=1,n2)/)),dim=2,ncopies=n1)
        t2d=spread(allknots((/(ix,ix=2+j2,j2+n2+1)/))-allknots((/(ix,ix=2,n2+1)/)),2,n1)
        t1n=transpose(spread(x1(indx),dim=2,ncopies=n2))-spread(allknots((/(ix,ix=1,n2)/)),dim=2,ncopies=n1)
        t2n=spread(allknots((/(ix,ix=2+j2,j2+n2+1)/)),dim=2,ncopies=n1)-transpose(spread(x1(indx),2,n2))

        t1=0.d0
        t2=t1
        if (j2==k) then
          t4=t1
          t5=t1
        endif
    
        ! If t1d~=0, t1=t1n/t1d
        ! If t1d==0, t1=0.d0
        temp=(t1d==0.d0)				
        t1d=merge(1.d0,t1d,temp)      
        t1=t1n/t1d
        t1=merge(0.d0,t1,temp)

        ! If t1d~=0, t4=3/t1d
        ! If t1d==0, t4=0.d0
        if (j2==k) then				
          t4=dble(k)/t1d					
          t4=merge(0.d0,t4,temp)
        endif
         
        ! If t2d~=0, t2=t2n/t2d
        ! If t2d==0, t2=0.d0
        temp=(t2d==0.d0)
        t2d=merge(1.d0,t2d,temp)
        t2=t2n/t2d					
        t2=merge(0.d0,t2,temp)		

        ! If t2d~=0, t5=-3/t2d
        ! If t2d==0, t5=0.d0
        if (j2==k) then
          t5=-dble(k)/t2d					
          t5=merge(0.d0,t5,temp)	
        endif

        temp2=0.d0
        temp3=0.d0
        if (j2==k) then
          temp2((/(ix,ix=1,n2-1)/),:)=t5((/(ix,ix=1,n2-1)/),:)
          temp3((/(ix,ix=1,n2-1)/),:)=b((/(ix,ix=2,n2)/),indx)
          db(:,indx)=t4*b(:,indx)+temp2*temp3
        endif

        temp2((/(ix,ix=1,n2-1)/),:)=t2((/(ix,ix=1,n2-1)/),:)
        temp3((/(ix,ix=1,n2-1)/),:)=b((/(ix,ix=2,n2)/),indx) 
        b(:,indx)=t1*b(:,indx)+temp2*temp3
      enddo
      jhi=jhi+n1
      deallocate(indx,t1d,t2d,t1n,t2n,t1,t2,t4,t5,temp,temp2,temp3)
    endif
  enddo

  ! If some values of x, still not processed, then compute b for these last x's.
  if (jhi<=nx) then
    n1=nx+1-jhi
    nonzero((/(ix,ix=jhi,nx)/))=nonzero((/(ix,ix=jhi,nx)/))+1
    allocate(indx(n1),t1d(n2,n1),t2d(n2,n1),t1n(n2,n1),t2n(n2,n1),t1(n2,n1),t2(n2,n1))
    allocate(t4(n2,n1),t5(n2,n1),temp(n2,n1),temp2(n2,n1),temp3(n2,n1))
    
    indx=(/(ix,ix=jhi,jhi+n1-1)/)

    b(nknots+k-1,indx)=1.d0
    
    do j2=1,k
      t1d=spread(allknots((/(ix,ix=1+j2,j2+n2)/))-allknots((/(ix,ix=1,n2)/)),dim=2,ncopies=n1)
      t2d=spread(allknots((/(ix,ix=2+j2,j2+n2+1)/))-allknots((/(ix,ix=2,n2+1)/)),2,n1)
      t1n=transpose(spread(x1(indx),dim=2,ncopies=n2))-spread(allknots((/(ix,ix=1,n2)/)),dim=2,ncopies=n1)
      t2n=spread(allknots((/(ix,ix=2+j2,j2+n2+1)/)),dim=2,ncopies=n1)-transpose(spread(x1(indx),2,n2))

      t1=0.d0
      t2=t1
      if (j2==k) then
        t4=t1
        t5=t1
      endif
      
      ! If t1d~=0, t1=t1n/t1d
      ! If t1d==0, t1=0.d0
      temp=(t1d==0.d0)				
      t1d=merge(1.d0,t1d,temp)      
      t1=t1n/t1d
      t1=merge(0.d0,t1,temp)

      ! If t1d~=0, t4=3/t1d
      ! If t1d==0, t4=0.d0
      if (j2==k) then				
        t4=dble(k)/t1d					
        t4=merge(0.d0,t4,temp)
      endif
    
      ! If t2d~=0, t2=t2n/t2d
      ! If t2d==0, t2=0.d0
      temp=(t2d==0.d0)
      t2d=merge(1.d0,t2d,temp)
      t2=t2n/t2d					
      t2=merge(0.d0,t2,temp)	
	
      ! If t2d~=0, t5=-3/t2d
      ! If t2d==0, t5=0.d0
      if (j2==k) then
        t5=-dble(k)/t2d					
    	t5=merge(0.d0,t5,temp)		
      endif

      temp2=0.d0
      temp3=0.d0

      if (j2==k) then
        temp2((/(ix,ix=1,n2-1)/),:)=t5((/(ix,ix=1,n2-1)/),:)
    	temp3((/(ix,ix=1,n2-1)/),:)=b((/(ix,ix=2,n2)/),indx)
        db(:,indx)=t4*b(:,indx)+temp2*temp3
      endif
      temp2((/(ix,ix=1,n2-1)/),:)=t2((/(ix,ix=1,n2-1)/),:)
      temp3((/(ix,ix=1,n2-1)/),:)=b((/(ix,ix=2,n2)/),indx)
      b(:,indx)=t1*b(:,indx)+temp2*temp3
    enddo
    deallocate(indx,t1d,t2d,t1n,t2n,t1,t2,t4,t5,temp,temp2,temp3)
  endif
  deallocate(x1)
  deallocate(allknots)

  if (sorted==0) then
    b=b(:,i01)
    db=db(:,i01)
    nonzero=nonzero(i01)
    deallocate(i00,i01)
  endif
end subroutine bspline2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! subroutine bspline3(x,knots,k,sorted,b,db,ddb,nonzero,loknots,hiknots)
! 
! Compute B-splines of degree k (order k+1) at the vector of points x.
! Also compute 1st and 2nd derivatives.
! 
! Required inputs:
!
!   real(dp) x(nx)                  Vector of points.
!   real(dp) knots(ngrid)           Knots for spline.
!   integer k		  	    Degree of B-spline.
!				    If k==3, the splines are cubic splines.
!   integer sorted		    1 if x is sorted, 0 otherwise.
!
! Required outputs:
!
!   real(dp) b(ngrid-1+k,nx)        b(i,j) =  B^k_i(x(j)).
!   real(dp) db(ngrid-1+k,nx)	    db(i,j) = dB^k_i(x(j))/dx.
!   real(dp) ddb(ngrid-1+k,nx)	    ddb(i,j) = d^2 B^k_i(x(j))/dx*dx.
!   integer(i4b) nonzero(nx,1)      nonzero(i) = index of first nonzero element in column i
!
! Optional inputs:
!
!   real(dp) loknots(k)             Extension of grid below knots(1).
!   real(dp) hiknots(k)             Extension of grid above knots(n).
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine bspline3(x,knots,k,sorted,b,db,ddb,nonzero,loknots,hiknots)
  use nrtype
  use nr, only: locate,indexx
  implicit none
  real(dp), intent(in)           :: knots(:)
  real(dp), intent(in)           :: x(:)
  integer(i4b), intent(in)       :: k
  integer(i4b), intent(in)       :: sorted
  real(dp), intent(out)          :: b(:,:),db(:,:),ddb(:,:)
  integer(i4b), intent(out)      :: nonzero(:)
  real(dp), optional, intent(in) :: loknots(:),hiknots(:)

  integer(i4b)              :: nx,nknots,n1,n2,jhi,j1,j2,ix
  integer(i4b), allocatable :: i00(:),i01(:),indx(:)
  logical, allocatable      :: temp(:,:)
  real(dp), allocatable     :: t1d(:,:),t2d(:,:),t1n(:,:),t2n(:,:),t1(:,:),t2(:,:)
  real(dp), allocatable     :: t4(:,:),t5(:,:),temp2(:,:),temp3(:,:)
  real(dp), allocatable     :: x1(:)
  real(dp), allocatable     :: allknots(:)

  nx=size(x)
  nknots=size(knots)
  allocate(x1(nx))
  allocate(allknots(nknots+2*k))
  x1=x

  if (sorted==0) then 
    allocate(i00(nx),i01(nx))
    i00=0

    ! Sort x in ascending order.
    call indexx(x,i00)		
    x1=x(i00)

    ! Indices for unsorting x.
    call indexx(i00,i01)		
  endif

  if (present(loknots) .and. present(hiknots)) then
    allknots((/(ix,ix=1,k)/))=loknots
    allknots((/(ix,ix=k+1,nknots+k)/))=knots
    allknots((/(ix,ix=nknots+k+1,nknots+2*k)/))=hiknots
  elseif (present(loknots) .and. .not.present(hiknots)) then
    allknots((/(ix,ix=1,k)/))=loknots
    allknots((/(ix,ix=k+1,nknots+k)/))=knots
    allknots((/(ix,ix=nknots+k+1,nknots+2*k)/))=knots(nknots)  
  elseif (.not.present(loknots) .and. present(hiknots)) then
    allknots((/(ix,ix=1,k)/))=knots(1)
    allknots((/(ix,ix=k+1,nknots+k)/))=knots
    allknots((/(ix,ix=nknots+k+1,nknots+2*k)/))=hiknots
  else
    allknots((/(ix,ix=1,k)/))=knots(1)
    allknots((/(ix,ix=k+1,nknots+k)/))=knots
    allknots((/(ix,ix=nknots+k+1,nknots+2*k)/))=knots(nknots)
  endif

  n2=nknots+k-1
  b=0.d0
  db=b
  ddb=b
  nonzero=0

  jhi=1

  do j1=k+2,nknots+k
    if (jhi>nx) exit
      n1=locate(x1((/(ix,ix=jhi,nx)/)),allknots(j1))
      nonzero((/(ix,ix=jhi,nx)/))=nonzero((/(ix,ix=jhi,nx)/))+1
      if (n1>0) then
        allocate(indx(n1),t1d(n2,n1),t2d(n2,n1),t1n(n2,n1),t2n(n2,n1),t1(n2,n1),t2(n2,n1))
        allocate(t4(n2,n1),t5(n2,n1),temp(n2,n1),temp2(n2,n1),temp3(n2,n1))
        indx=(/(ix,ix=jhi,jhi+n1-1)/)
        b(j1-1,indx)=1.d0
        do j2=1,k
          t1d=spread(allknots((/(ix,ix=1+j2,j2+n2)/))-allknots((/(ix,ix=1,n2)/)),dim=2,ncopies=n1)
          t2d=spread(allknots((/(ix,ix=2+j2,j2+n2+1)/))-allknots((/(ix,ix=2,n2+1)/)),2,n1)
          t1n=transpose(spread(x1(indx),dim=2,ncopies=n2))-spread(allknots((/(ix,ix=1,n2)/)),dim=2,ncopies=n1)
          t2n=spread(allknots((/(ix,ix=2+j2,j2+n2+1)/)),dim=2,ncopies=n1)-transpose(spread(x1(indx),2,n2))

          t1=0.d0
          t2=t1
          if (j2>k-2) then
            t4=t1
            t5=t1
          endif
      
          ! If t1d~=0, t1=t1n/t1d
          ! If t1d==0, t1=0.d0
          temp=(t1d==0.d0)				
          t1d=merge(1.d0,t1d,temp)      
          t1=t1n/t1d
          t1=merge(0.d0,t1,temp)


          if (j2==k-1) then

            ! If t1d~=0, t4=(k-1)/t1d
            ! If t1d==0, t4=0.d0
            t4=dble(k-1)/t1d
            t4=merge(0.d0,t4,temp)

          elseif (j2==k) then				

            ! If t1d~=0, t4=k/t1d
            ! If t1d==0, t4=0.d0
            t4=dble(k)/t1d					
            t4=merge(0.d0,t4,temp)

  	  endif
          
          ! If t2d~=0, t2=t2n/t2d
          ! If t2d==0, t2=0.d0
          temp=(t2d==0.d0)
          t2d=merge(1.d0,t2d,temp)
          t2=t2n/t2d					
          t2=merge(0.d0,t2,temp)		

          if (j2==k-1) then

            ! If t2d~=0, t5=-(k-1)/t2d
            ! If t2d==0, t5=0.d0
            t5=-dble(k-1)/t2d
            t5=merge(0.d0,t5,temp)

          elseif (j2==k) then

            ! If t2d~=0, t5=-k/t2d
            ! If t2d==0, t5=0.d0
            t5=-dble(k)/t2d					
            t5=merge(0.d0,t5,temp)	

          endif

          temp2=0.d0
          temp3=0.d0
          if (j2==k-1) then
            temp2((/(ix,ix=1,n2-1)/),:)=t5((/(ix,ix=1,n2-1)/),:)
            temp3((/(ix,ix=1,n2-1)/),:)=b((/(ix,ix=2,n2)/),indx)
            ddb(:,indx)=t4*b(:,indx)+temp2*temp3
          elseif (j2==k) then
            temp2((/(ix,ix=1,n2-1)/),:)=t5((/(ix,ix=1,n2-1)/),:)
	        temp3((/(ix,ix=1,n2-1)/),:)=b((/(ix,ix=2,n2)/),indx)
            db(:,indx)=t4*b(:,indx)+temp2*temp3
            temp3((/(ix,ix=1,n2-1)/),:)=ddb((/(ix,ix=2,n2)/),indx)
            ddb(:,indx)=t4*ddb(:,indx)+temp2*temp3
          endif

          temp2((/(ix,ix=1,n2-1)/),:)=t2((/(ix,ix=1,n2-1)/),:)
          temp3((/(ix,ix=1,n2-1)/),:)=b((/(ix,ix=2,n2)/),indx) 
          b(:,indx)=t1*b(:,indx)+temp2*temp3
        enddo
      jhi=jhi+n1
      deallocate(indx,t1d,t2d,t1n,t2n,t1,t2,t4,t5,temp,temp2,temp3)
    endif
  enddo

  if (jhi<=nx) then
    n1=nx+1-jhi
    nonzero((/(ix,ix=jhi,nx)/))=nonzero((/(ix,ix=jhi,nx)/))+1
    allocate(indx(n1),t1d(n2,n1),t2d(n2,n1),t1n(n2,n1),t2n(n2,n1),t1(n2,n1),t2(n2,n1))
    allocate(t4(n2,n1),t5(n2,n1),temp(n2,n1),temp2(n2,n1),temp3(n2,n1))
    
    indx=(/(ix,ix=jhi,jhi+n1-1)/)

    b(nknots+k-1,indx)=1.d0
    
    do j2=1,k
      t1d=spread(allknots((/(ix,ix=1+j2,j2+n2)/))-allknots((/(ix,ix=1,n2)/)),dim=2,ncopies=n1)
      t2d=spread(allknots((/(ix,ix=2+j2,j2+n2+1)/))-allknots((/(ix,ix=2,n2+1)/)),2,n1)
      t1n=transpose(spread(x1(indx),dim=2,ncopies=n2))-spread(allknots((/(ix,ix=1,n2)/)),dim=2,ncopies=n1)
      t2n=spread(allknots((/(ix,ix=2+j2,j2+n2+1)/)),dim=2,ncopies=n1)-transpose(spread(x1(indx),2,n2))

      t1=0.d0
      t2=t1
      if (j2>k-2) then
        t4=t1
        t5=t1
      endif
      
      ! If t1d~=0, t1=t1n/t1d
      ! If t1d==0, t1=0.d0
      temp=(t1d==0.d0)				
      t1d=merge(1.d0,t1d,temp)      
      t1=t1n/t1d
      t1=merge(0.d0,t1,temp)

      if (j2==k-1) then

        ! If t1d~=0, t4=(k-1)/t1d
        ! If t1d==0, t4=0.d0
        t4=dble(k-1)/t1d
        t4=merge(0.d0,t4,temp)

      elseif (j2==k) then				

        ! If t1d~=0, t4=k/t1d
        ! If t1d==0, t4=0.d0
        t4=dble(3)/t1d					
        t4=merge(0.d0,t4,temp)

      endif
    
      ! If t2d~=0, t2=t2n/t2d
      ! If t2d==0, t2=0.d0
      temp=(t2d==0.d0)
      t2d=merge(1.d0,t2d,temp)
      t2=t2n/t2d					
      t2=merge(0.d0,t2,temp)		

      if (j2==k-1) then

        ! If t2d~=0, t5=-(k-1)/t2d
        ! If t2d==0, t5=0.d0
        t5=-dble(k-1)/t2d					
    	t5=merge(0.d0,t5,temp)
		
      elseif (j2==k) then

        ! If t2d~=0, t5=-k/t2d
        ! If t2d==0, t5=0.d0
        t5=-dble(k)/t2d					
    	t5=merge(0.d0,t5,temp)		

      endif

      temp2=0.d0
      temp3=0.d0

      if (j2==k-1) then
        temp2((/(ix,ix=1,n2-1)/),:)=t5((/(ix,ix=1,n2-1)/),:)
    	temp3((/(ix,ix=1,n2-1)/),:)=b((/(ix,ix=2,n2)/),indx)
        db(:,indx)=t4*b(:,indx)+temp2*temp3
      elseif (j2==k) then
        temp2((/(ix,ix=1,n2-1)/),:)=t5((/(ix,ix=1,n2-1)/),:)
    	temp3((/(ix,ix=1,n2-1)/),:)=b((/(ix,ix=2,n2)/),indx)
        db(:,indx)=t4*b(:,indx)+temp2*temp3
    	temp3((/(ix,ix=1,n2-1)/),:)=ddb((/(ix,ix=2,n2)/),indx)
        ddb(:,indx)=t4*ddb(:,indx)+temp2*temp3
      endif
      temp2((/(ix,ix=1,n2-1)/),:)=t2((/(ix,ix=1,n2-1)/),:)
      temp3((/(ix,ix=1,n2-1)/),:)=b((/(ix,ix=2,n2)/),indx)
      b(:,indx)=t1*b(:,indx)+temp2*temp3
    enddo
    deallocate(indx,t1d,t2d,t1n,t2n,t1,t2,t4,t5,temp,temp2,temp3)
  endif

  deallocate(x1)
  deallocate(allknots)
  if (sorted==0) then
    b=b(:,i01)
    db=db(:,i01)
    ddb=ddb(:,i01)
    nonzero=nonzero(i01)
    deallocate(i00,i01)
  endif
end subroutine bspline3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!  subroutine bspline4(b,x,xg,k,i,m,xg0,xg1)
!
! Compute i'th B-spline of degree k (order k+1) at the vector of points x.
! If k==3, the splines are cubic splines.
! The values of x lie in the interval [xg(m),xg(m+1)].
! 
! b     =   real(dp) (nx x 1)    b =  B^k_i(x)
! x     =   real(dp) (nx x 1)    vector of points
! xg    =   real(dp) (ng x 1)    grid for spline
! k     =   integer (1 x 1)     degree of B-spline
! i     =   integer (1 x 1)     index of B-Spline to compute. i is in {1,...,ng+k-1}
! m     =   integer (1 x 1)	    index for x interval. All values of x lie in interval [xg(m),xg(m+1)].
! xg0   =   real(dp) (k x 1)     optional, extension of grid below xg(1).
! xg1   =   real(dp) (k x 1)     optional, extension of grid above xg(n).
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine bspline4(b,x,xg,k,i,m,xg00,xg10)
  use nrtype
  implicit none
  real(dp),     intent(in)          :: x(:),xg(:)
  integer(i4b), intent(in)          :: k,i,m
  real(dp),     intent(out)         :: b(:)
  real(dp),     intent(in),optional :: xg00(:),xg10(:)  

  integer(i4b)          :: nx,ng,i1,i2,n2,ix
  real(dp)              :: dx
  real(dp), allocatable :: xg0(:)
  real(dp), allocatable :: xg1(:)
  real(dp), allocatable :: xg2(:)
  real(dp), allocatable :: b1(:)
  real(dp), allocatable :: t1(:)
  real(dp), allocatable :: t2(:)
  
  nx=size(x) 
  ng=size(xg)
  allocate(xg0(k))
  allocate(xg1(k))
  allocate(xg2(ng+k+k))
  allocate(b1(k+1))
  allocate(t1(k))
  allocate(t2(k))
  
  if (present(xg00)) then
    xg0=xg00
  else
    dx=xg(2)-xg(1)
    do i1=k,1,-1
      xg0(1+k-i1)=xg(1)-dx*dble(i1)
    end do      
  end if
    
  if (present(xg10)) then
    xg1=xg10
  else
    dx=xg(ng)-xg(ng-1)
    do i1=1,k
      xg1(i1)=xg(ng)+dx*dble(i1)
    end do      
  end if
  
  n2=k+1
  xg2((/(ix,ix=1,k)/))=xg0
  xg2((/(ix,ix=k+1,k+ng)/))=xg
  xg2((/(ix,ix=k+1+ng,ng+k+k)/))=xg1

  if ((1+m+k-i<1) .or. (1+m+k-i>k+1)) then
    b=0.d0
    return
  end if
  
  do i1=1,nx
    b1=0.d0
    b1(1+m+k-i)=1.d0
    do i2=1,k
      t1=1.d0
      t2=1.d0
      t1((/(ix,ix=1,n2-i2)/)) &
           = (x(i1)-xg2((/(ix,ix=i,i+k-i2)/))) &
           / (xg2((/(ix,ix=i+i2,i+k)/))-xg2((/(ix,ix=i,i+k-i2)/)))    
      t2((/(ix,ix=1,n2-i2)/)) &
           = (xg2((/(ix,ix=i+i2+1,i+k+1)/))-x(i1)) &
            / (xg2((/(ix,ix=i2+i+1,i+k+1)/))-xg2((/(ix,ix=i+1,i+k+1-i2)/)))
      b1((/(ix,ix=1,n2-i2)/)) = t1((/(ix,ix=1,n2-i2)/))*b1((/(ix,ix=1,n2-i2)/)) &
                              + t2((/(ix,ix=1,n2-i2)/))*b1((/(ix,ix=2,n2-i2+1)/))
    end do
    b(i1)=b1(1)
  end do

  deallocate(xg0)
  deallocate(xg1)
  deallocate(xg2)
  deallocate(b1)
  deallocate(t1)
  deallocate(t2)
end subroutine bspline4

subroutine ComputeKnotAverages(knots,k,x)

use nrtype

implicit none
real(dp),     intent(in)  :: knots(:)
integer(i4b), intent(in)  :: k
real(dp),     intent(out) :: x(:)

integer(i4b)          :: nknots,i1,ix
real(dp), allocatable :: temp(:)

nknots = size(knots,1)

allocate(temp(nknots+2*k))

temp = 0.0d0
temp((/(ix,ix=1,k)/)) = knots(1)
temp(k+(/(ix,ix=1,nknots)/)) = knots
temp(k+nknots+(/(ix,ix=1,k)/)) = knots(nknots)

x = 0.0d0
do i1=1,nknots+k-1
x(i1) = sum(temp((/(ix,ix=i1,i1+k)/)))/dble(k+1)
end do

deallocate(temp)

end subroutine ComputeKnotAverages

end module SplineTools
