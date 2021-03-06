!Modification history
! 19sep2007 lpn deleted indexx_sp which was causing problems.
!
! SUBROUTINE indexx_dp(arr,index)		Sort data (double precision)
! SUBROUTINE indexx_i4b(arr,index)		Sort data (integer)
! FUNCTION locate(xx,x)				
! SUBROUTINE ludcmp(a,indx,d)	 		LU decomposition
! SUBROUTINE lubksb(a,indx,b)  			Inverse after LU decomposition
! FUNCTION ran_d(idum)
! FUNCTION erfc_d(x)
! FUNCTION erfc_v(x)
! FUNCTION erf_d(x)
! FUNCTION erf_v(x)
! FUNCTION gammp_d(a,x)
! FUNCTION gammp_v(a,x)
! FUNCTION gammq_d(a,x)
! FUNCTION gammq_v(a,x)
! FUNCTION gcf_d(a,x,gln)
! FUNCTION gcf_v(a,x,gln)
! FUNCTION gser_d(a,x,gln)
! FUNCTION gser_v(a,x,gln)
! FUNCTION gammln_d(xx)
! FUNCTION gammln_v(xx)
! SUBROUTINE gasdev_d(harvest)
! SUBROUTINE gasdev_v(harvest)
! SUBROUTINE ran1_d(harvest)
! SUBROUTINE ran1_v(harvest)
! SUBROUTINE gauher(x,w)			Gauss-Hermite Quadrature nodes
! SUBROUTINE gauleg(x1,x2,x,w)			Gauss-Legendre Quadrature nodes
! SUBROUTINE choldc(a,p)			Cholesky decomposition	
! subroutine invchol(a,inva)			Inverse of Cholesky decomposition
! subroutine mnewt(ntrial,x,n,tolx,tolf)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 	SUBROUTINE indexx_dp(arr,indx)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE indexx_dp(arr,indx)
use nrtype
USE nrutil, ONLY : arth,assert_eq,nrerror,swap
IMPLICIT NONE
REAL(DP), DIMENSION(:), INTENT(IN) :: arr
INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
REAL(DP) :: a
INTEGER(I4B) :: n,k,i,j,indxt,jstack,l,r
INTEGER(I4B), DIMENSION(NSTACK) :: istack
n=assert_eq(size(indx),size(arr),'indexx_dp')
indx=arth(1,1,n)
jstack=0
l=1
r=n
do
  if (r-l < NN) then
    do j=l+1,r
      indxt=indx(j)
   	  a=arr(indxt)
	  do i=j-1,l,-1
		  if (arr(indx(i)) <= a) exit
		  indx(i+1)=indx(i)
		end do
		indx(i+1)=indxt
	  end do
	  if (jstack == 0) RETURN
	  r=istack(jstack)
	  l=istack(jstack-1)
	  jstack=jstack-2
	else
	  k=(l+r)/2
	  call swap(indx(k),indx(l+1))
	  call icomp_xchg(indx(l),indx(r))
	  call icomp_xchg(indx(l+1),indx(r))
	  call icomp_xchg(indx(l),indx(l+1))
	  i=l+1
	  j=r
	  indxt=indx(l+1)
	  a=arr(indxt)
	  do
	    do
		  i=i+1
		  if (arr(indx(i)) >= a) exit
		end do
		do
		  j=j-1
	      if (arr(indx(j)) <= a) exit
		end do
		if (j < i) exit
		call swap(indx(i),indx(j))
	  end do
	  indx(l+1)=indx(j)
	  indx(j)=indxt
	  jstack=jstack+2
	  if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
	  if (r-i+1 >= j-l) then
	    istack(jstack)=r
		istack(jstack-1)=i
		r=j-1
	  else
	    istack(jstack)=j-1
		istack(jstack-1)=l
		l=i
	  end if
	end if
  end do
  CONTAINS
!BL
  SUBROUTINE icomp_xchg(i,j)
    INTEGER(I4B), INTENT(INOUT) :: i,j
	INTEGER(I4B) :: swp
	if (arr(j) < arr(i)) then
	  swp=i
	  i=j
	  j=swp
	end if
  END SUBROUTINE icomp_xchg
END SUBROUTINE indexx_dp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 	SUBROUTINE indexx_i4b(arr,indx)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE indexx_i4b(iarr,indx)
use nrtype
USE nrutil, ONLY : arth,assert_eq,nrerror,swap
IMPLICIT NONE
INTEGER(I4B), DIMENSION(:), INTENT(IN) :: iarr
INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50
INTEGER(I4B) :: a
INTEGER(I4B) :: n,k,i,j,indxt,jstack,l,r
INTEGER(I4B), DIMENSION(NSTACK) :: istack

n=assert_eq(size(indx),size(iarr),'indexx_i4b')
	indx=arth(1,1,n)
	jstack=0
	l=1
	r=n
	do
		if (r-l < NN) then
			do j=l+1,r
				indxt=indx(j)
				a=iarr(indxt)
				do i=j-1,l,-1
					if (iarr(indx(i)) <= a) exit
					indx(i+1)=indx(i)
				end do
				indx(i+1)=indxt
			end do
			if (jstack == 0) RETURN
			r=istack(jstack)
			l=istack(jstack-1)
			jstack=jstack-2
		else
			k=(l+r)/2
			call swap(indx(k),indx(l+1))
			call icomp_xchg(indx(l),indx(r))
			call icomp_xchg(indx(l+1),indx(r))
			call icomp_xchg(indx(l),indx(l+1))
			i=l+1
			j=r
			indxt=indx(l+1)
			a=iarr(indxt)
			do
				do
					i=i+1
					if (iarr(indx(i)) >= a) exit
				end do
				do
					j=j-1
					if (iarr(indx(j)) <= a) exit
				end do
				if (j < i) exit
				call swap(indx(i),indx(j))
			end do
			indx(l+1)=indx(j)
			indx(j)=indxt
			jstack=jstack+2
			if (jstack > NSTACK) call nrerror('indexx: NSTACK too small')
			if (r-i+1 >= j-l) then
				istack(jstack)=r
				istack(jstack-1)=i
				r=j-1
			else
				istack(jstack)=j-1
				istack(jstack-1)=l
				l=i
			end if
		end if
	end do
	CONTAINS
!BL
	SUBROUTINE icomp_xchg(i,j)
	INTEGER(I4B), INTENT(INOUT) :: i,j
	INTEGER(I4B) :: swp
	if (iarr(j) < iarr(i)) then
		swp=i
		i=j
		j=swp
	end if
	END SUBROUTINE icomp_xchg
	END SUBROUTINE indexx_i4b

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 	FUNCTION locate(xx,x)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FUNCTION locate(xx,x)
use nrtype
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: xx
  REAL(DP), INTENT(IN) :: x
  INTEGER(I4B) :: locate
  INTEGER(I4B) :: n,jl,jm,ju
  LOGICAL :: ascnd
  n=size(xx)
  ascnd = (xx(n) >= xx(1))
  jl=0
  ju=n+1
  do
    if (ju-jl <= 1) exit
	jm=(ju+jl)/2
	if (ascnd .eqv. (x >= xx(jm))) then
	  jl=jm
	else
	  ju=jm
	end if
  end do
  if (x == xx(1)) then
    locate=1
  else if (x == xx(n)) then
	locate=n-1
  else
	locate=jl
  end if
END FUNCTION locate

SUBROUTINE ludcmp(a,indx,d)
  USE nrtype
  USE nrutil, ONLY : assert_eq,imaxloc,nrerror,outerprod,swap
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
  INTEGER(I4B), DIMENSION(:), INTENT(OUT) :: indx
  REAL(DP), INTENT(OUT) :: d
  REAL(DP), DIMENSION(size(a,1)) :: vv
  REAL(DP), PARAMETER :: TINY=1.0e-20_dp
  INTEGER(I4B) :: j,n,imax
  n=assert_eq(size(a,1),size(a,2),size(indx),'ludcmp')
  d=1.0d0
  vv=maxval(abs(a),dim=2)
  if (any(vv == 0.0)) call nrerror('singular matrix in ludcmp')
  vv=1.0_dp/vv
  do j=1,n
  	imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j)))
	if (j /= imax) then
		call swap(a(imax,:),a(j,:))
		d=-d
		vv(imax)=vv(j)
	end if
	indx(j)=imax
	if (a(j,j) == 0.0d0) a(j,j)=TINY
    a(j+1:n,j)=a(j+1:n,j)/a(j,j)
	a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
  end do
END SUBROUTINE ludcmp

SUBROUTINE lubksb(a,indx,b)
  USE nrtype
  USE nrutil, ONLY : assert_eq
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
  INTEGER(I4B), DIMENSION(:), INTENT(IN) :: indx
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: b
  INTEGER(I4B) :: i,n,ii,ll
  REAL(DP) :: summ
  n=assert_eq(size(a,1),size(a,2),size(indx),'lubksb')
  ii=0
  do i=1,n
 	ll=indx(i)
	summ=b(ll)
	b(ll)=b(i)
	if (ii /= 0) then
		summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
	else if (summ /= 0.0) then
		ii=i
	end if
	b(i)=summ
  end do
  do i=n,1,-1
 	b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
  end do
END SUBROUTINE lubksb

FUNCTION ran_d(idum)
  use nrtype
  IMPLICIT NONE
  INTEGER, PARAMETER :: K4B=selected_int_kind(9)
  INTEGER(K4B), INTENT(INOUT) :: idum
  REAL(DP) :: ran_d
  INTEGER(K4B), PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
  REAL(DP), SAVE :: am
  INTEGER(K4B), SAVE :: ix=-1,iy=-1,k
  if (idum <= 0 .or. iy < 0) then
    am=nearest(1.0,-1.0)/IM
    iy=ior(ieor(888889999,abs(idum)),1)
    ix=ieor(777755555,abs(idum))
    idum=abs(idum)+1
  end if
  ix=ieor(ix,ishft(ix,13))
  ix=ieor(ix,ishft(ix,-17))
  ix=ieor(ix,ishft(ix,5))
  k=iy/IQ
  iy=IA*(iy-k*IQ)-IR*k
  if (iy < 0) iy=iy+IM
  ran_d=am*ior(iand(IM,ieor(ix,iy)),1)
END FUNCTION ran_d

FUNCTION erfc_d(x)
  USE nrtype
  use nr, only : gammp,gammq
  IMPLICIT NONE
  REAL(dP), INTENT(IN) :: x
  REAL(dP) :: erfc_d
  erfc_d=merge(1.0_dp+gammp(0.5_dp,x**2),gammq(0.5_dp,x**2), x < 0.0d0)
END FUNCTION erfc_d

FUNCTION erfc_v(x)
  USE nrtype
  use nr, only : gammp,gammq
  IMPLICIT NONE
  REAL(DP), DIMENSION(:), INTENT(IN) :: x
  REAL(DP), DIMENSION(size(x)) :: erfc_v
  LOGICAL(LGT), DIMENSION(size(x)) :: mask
  mask = (x < 0.0d0)
  erfc_v=merge(1.0_dp+gammp(spread(0.5_dp,1,size(x)), &
 		merge(x,0.0_dp,mask)**2),gammq(spread(0.5_dp,1,size(x)), &
 		merge(x,0.0_dp,.not. mask)**2),mask)
END FUNCTION erfc_v

FUNCTION erf_d(x)
	USE nrtype
	USE nr, ONLY : gammp
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x
	REAL(DP) :: erf_d
	erf_d=gammp(0.5_dp,x**2)
	if (x < 0.0d0) erf_d=-erf_d
END FUNCTION erf_d

FUNCTION erf_v(x)
	USE nrtype
	USE nr, ONLY : gammp
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: x
	REAL(DP), DIMENSION(size(x)) :: erf_v
	erf_v=gammp(spread(0.5_dp,1,size(x)),x**2)
	where (x < 0.0d0) erf_v=-erf_v
END FUNCTION erf_v

FUNCTION gammp_d(a,x)
	USE nrtype; USE nrutil, ONLY : assert
    USE nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,x
	REAL(DP) :: gammp_d
	call assert( x >= 0.0d0,  a > 0.0d0, 'gammp_d args')
	if (x<a+1.0_dp) then
		gammp_d=gser(a,x)
	else
		gammp_d=1.0_dp-gcf(a,x)
	end if
END FUNCTION gammp_d

FUNCTION gammp_v(a,x)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(dP), DIMENSION(:), INTENT(IN) :: a,x
	REAL(dP), DIMENSION(size(x)) :: gammp_v
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(a),size(x),'gammp_v')
	call assert( all(x >= 0.0d0),  all(a > 0.0d0), 'gammp_v args')
	mask = (x<a+1.0_dp)
	gammp_v=merge(gser(a,merge(x,0.0_dp,mask)), &
		1.0_dp-gcf(a,merge(x,0.0_dp,.not. mask)),mask)
END FUNCTION gammp_v

FUNCTION gammq_d(a,x)
	USE nrtype; USE nrutil, ONLY : assert
	USE nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,x
	REAL(DP) :: gammq_d
	call assert( x >= 0.0d0,  a > 0.0d0, 'gammq_d args')
	if (x<a+1.0_dp) then
		gammq_d=1.0_dp-gser(a,x)
	else
		gammq_d=gcf(a,x)
	end if
END FUNCTION gammq_d

FUNCTION gammq_v(a,x)
	USE nrtype; USE nrutil, ONLY : assert,assert_eq
	USE nr, ONLY : gcf,gser
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
	REAL(DP), DIMENSION(size(a)) :: gammq_v
	LOGICAL(LGT), DIMENSION(size(x)) :: mask
	INTEGER(I4B) :: ndum
	ndum=assert_eq(size(a),size(x),'gammq_v')
	call assert( all(x >= 0.0d0),  all(a > 0.0d0), 'gammq_v args')
	mask = (x<a+1.0_dp)
	gammq_v=merge(1.0_dp-gser(a,merge(x,0.0_dp,mask)), &
		gcf(a,merge(x,0.0_dp,.not. mask)),mask)
END FUNCTION gammq_v

FUNCTION gcf_d(a,x,gln)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,x
	REAL(DP), OPTIONAL, INTENT(OUT) :: gln
	REAL(DP) :: gcf_d
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(DP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
	INTEGER(I4B) :: i
	REAL(DP) :: an,b,c,d,del,h
	if (x == 0.0d0) then
		gcf_d=1.0d0
		RETURN
	end if
	b=x+1.0_dp-a
	c=1.0_dp/FPMIN
	d=1.0_dp/b
	h=d
	do i=1,ITMAX
		an=-i*(i-a)
		b=b+2.0_dp
		d=an*d+b
		if (abs(d) < FPMIN) d=FPMIN
		c=b+an/c
		if (abs(c) < FPMIN) c=FPMIN
		d=1.0_dp/d
		del=d*c
		h=h*del
		if (abs(del-1.0_dp) <= EPS) exit
	end do
	if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_d')
	if (present(gln)) then
		gln=gammln(a)
		gcf_d=exp(-x+a*log(x)-gln)*h
	else
		gcf_d=exp(-x+a*log(x)-gammln(a))*h
	end if
END FUNCTION gcf_d

FUNCTION gcf_v(a,x,gln)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
	REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
	REAL(DP), DIMENSION(size(a)) :: gcf_v
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(DP), PARAMETER :: EPS=epsilon(x),FPMIN=tiny(x)/EPS
	INTEGER(I4B) :: i
	REAL(DP), DIMENSION(size(a)) :: an,b,c,d,del,h
	LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
	i=assert_eq(size(a),size(x),'gcf_v')
	zero=(x == 0.0d0)
	where (zero)
		gcf_v=1.0_dp
	elsewhere
		b=x+1.0_dp-a
		c=1.0_dp/FPMIN
		d=1.0_dp/b
		h=d
	end where
	converged=zero
	do i=1,ITMAX
		where (.not. converged)
			an=-i*(i-a)
			b=b+2.0_dp
			d=an*d+b
			d=merge(FPMIN,d, abs(d)<FPMIN )
			c=b+an/c
			c=merge(FPMIN,c, abs(c)<FPMIN )
			d=1.0_dp/d
			del=d*c
			h=h*del
			converged = (abs(del-1.0_dp)<=EPS)
		end where
		if (all(converged)) exit
	end do
	if (i > ITMAX) call nrerror('a too large, ITMAX too small in gcf_v')
	if (present(gln)) then
		if (size(gln) < size(a)) call &
			nrerror('gser: Not enough space for gln')
		gln=gammln(a)
		where (.not. zero) gcf_v=exp(-x+a*log(x)-gln)*h
	else
		where (.not. zero) gcf_v=exp(-x+a*log(x)-gammln(a))*h
	end if
END FUNCTION gcf_v

FUNCTION gser_d(a,x,gln)
	USE nrtype; USE nrutil, ONLY : nrerror
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: a,x
	REAL(DP), OPTIONAL, INTENT(OUT) :: gln
	REAL(DP) :: gser_d
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(dP), PARAMETER :: EPS=epsilon(x)
	INTEGER(I4B) :: n
	REAL(dP) :: ap,del,summ
	if (x == 0.0d0) then
		gser_d=0.0d0
		RETURN
	end if
	ap=a
	summ=1.0_dp/a
	del=summ
	do n=1,ITMAX
		ap=ap+1.0_dp
		del=del*x/ap
		summ=summ+del
		if (abs(del) < abs(summ)*EPS) exit
	end do
	if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_d')
	if (present(gln)) then
		gln=gammln(a)
		gser_d=summ*exp(-x+a*log(x)-gln)
	else
		gser_d=summ*exp(-x+a*log(x)-gammln(a))
	end if
END FUNCTION gser_d

FUNCTION gser_v(a,x,gln)
	USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
	USE nr, ONLY : gammln
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(IN) :: a,x
	REAL(DP), DIMENSION(:), OPTIONAL, INTENT(OUT) :: gln
	REAL(DP), DIMENSION(size(a)) :: gser_v
	INTEGER(I4B), PARAMETER :: ITMAX=100
	REAL(DP), PARAMETER :: EPS=epsilon(x)
	INTEGER(I4B) :: n
	REAL(DP), DIMENSION(size(a)) :: ap,del,summ
	LOGICAL(LGT), DIMENSION(size(a)) :: converged,zero
	n=assert_eq(size(a),size(x),'gser_v')
	zero=(x == 0.0d0)
	where (zero) gser_v=0.0d0
	ap=a
	summ=1.0_dp/a
	del=summ
	converged=zero
	do n=1,ITMAX
		where (.not. converged)
			ap=ap+1.0_dp
			del=del*x/ap
			summ=summ+del
			converged = (abs(del) < abs(summ)*EPS)
		end where
		if (all(converged)) exit
	end do
	if (n > ITMAX) call nrerror('a too large, ITMAX too small in gser_v')
	if (present(gln)) then
		if (size(gln) < size(a)) call &
			nrerror('gser: Not enough space for gln')
		gln=gammln(a)
		where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gln)
	else
		where (.not. zero) gser_v=summ*exp(-x+a*log(x)-gammln(a))
	end if
END FUNCTION gser_v

FUNCTION gammln_d(xx)
	USE nrtype; USE nrutil, ONLY : arth,assert
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: xx
	REAL(DP) :: gammln_d
	REAL(DP) :: tmp,x
	REAL(DP) :: stp = 2.5066282746310005_dp
	REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
		-86.50532032941677_dp,24.01409824083091_dp,&
		-1.231739572450155_dp,0.1208650973866179e-2_dp,&
		-0.5395239384953e-5_dp/)
	call assert(xx > 0.0d0, 'gammln_d arg')
	x=xx
	tmp=x+5.5_dp
	tmp=(x+0.5_dp)*log(tmp)-tmp
	gammln_d=tmp+log(stp*(1.000000000190015_dp+&
		sum(coef(:)/arth(x+1.0_dp,1.0_dp,size(coef))))/x)
END FUNCTION gammln_d

FUNCTION gammln_v(xx)
	USE nrtype; USE nrutil, ONLY: assert
	IMPLICIT NONE
	INTEGER(I4B) :: i
	REAL(DP), DIMENSION(:), INTENT(IN) :: xx
	REAL(DP), DIMENSION(size(xx)) :: gammln_v
	REAL(DP), DIMENSION(size(xx)) :: ser,tmp,x,y
	REAL(DP) :: stp = 2.5066282746310005_dp
	REAL(DP), DIMENSION(6) :: coef = (/76.18009172947146_dp,&
		-86.50532032941677_dp,24.01409824083091_dp,&
		-1.231739572450155_dp,0.1208650973866179e-2_dp,&
		-0.5395239384953e-5_dp/)
	if (size(xx) == 0) RETURN
	call assert(all(xx > 0.0d0), 'gammln_v arg')
	x=xx
	tmp=x+5.5_dp
	tmp=(x+0.5_dp)*log(tmp)-tmp
	ser=1.000000000190015_dp
	y=x
	do i=1,size(coef)
		y=y+1.0_dp
		ser=ser+coef(i)/y
	end do
	gammln_v=tmp+log(stp*ser/x)
END FUNCTION gammln_v

SUBROUTINE gasdev_d(harvest)
	USE nrtype
	USE nr, ONLY : ran1
	IMPLICIT NONE
	REAL(DP), INTENT(OUT) :: harvest
	REAL(DP) :: rsq,v1,v2
	REAL(DP), SAVE :: g
	LOGICAL, SAVE :: gaus_stored=.false.
	if (gaus_stored) then
		harvest=g
		gaus_stored=.false.
	else
		do
			call ran1(v1)
			call ran1(v2)
			v1=2.0_dp*v1-1.0_dp
			v2=2.0_dp*v2-1.0_dp
			rsq=v1**2+v2**2
			if (rsq > 0.0d0 .and. rsq < 1.0d0) exit
		end do
		rsq=sqrt(-2.0_dp*log(rsq)/rsq)
		harvest=v1*rsq
		g=v2*rsq
		gaus_stored=.true.
	end if
END SUBROUTINE gasdev_d

SUBROUTINE gasdev_v(harvest)
	USE nrtype; USE nrutil, ONLY : array_copy
	USE nr, ONLY : ran1
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(OUT) :: harvest
	REAL(DP), DIMENSION(size(harvest)) :: rsq,v1,v2
	REAL(DP), ALLOCATABLE, DIMENSION(:), SAVE :: g
	INTEGER(I4B) :: n,ng,nn,m
	INTEGER(I4B), SAVE :: last_allocated=0
	LOGICAL, SAVE :: gaus_stored=.false.
	LOGICAL, DIMENSION(size(harvest)) :: mask
	n=size(harvest)
	if (n /= last_allocated) then
		if (last_allocated /= 0) deallocate(g)
		allocate(g(n))
		last_allocated=n
		gaus_stored=.false.
	end if
	if (gaus_stored) then
		harvest=g
		gaus_stored=.false.
	else
		ng=1
		do
			if (ng > n) exit
			call ran1(v1(ng:n))
			call ran1(v2(ng:n))
			v1(ng:n)=2.0_dp*v1(ng:n)-1.0_dp
			v2(ng:n)=2.0_dp*v2(ng:n)-1.0_dp
			rsq(ng:n)=v1(ng:n)**2+v2(ng:n)**2
			mask(ng:n)=(rsq(ng:n)>0.0d0 .and. rsq(ng:n)<1.0)
			call array_copy(pack(v1(ng:n),mask(ng:n)),v1(ng:),nn,m)
			v2(ng:ng+nn-1)=pack(v2(ng:n),mask(ng:n))
			rsq(ng:ng+nn-1)=pack(rsq(ng:n),mask(ng:n))
			ng=ng+nn
		end do
		rsq=sqrt(-2.0_dp*log(rsq)/rsq)
		harvest=v1*rsq
		g=v2*rsq
		gaus_stored=.true.
	end if
END SUBROUTINE gasdev_v

SUBROUTINE ran1_d(harvest)
	USE nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
		iran0,jran0,kran0,nran0,mran0,rans
	IMPLICIT NONE
	REAL(DP), INTENT(OUT) :: harvest
	if (lenran < 1) call ran_init(1)
	rans=iran0-kran0
	if (rans < 0) rans=rans+2147483579_k4b
	iran0=jran0
	jran0=kran0
	kran0=rans
	nran0=ieor(nran0,ishft(nran0,13))
	nran0=ieor(nran0,ishft(nran0,-17))
	nran0=ieor(nran0,ishft(nran0,5))
	if (nran0 == 1) nran0=270369_k4b
	mran0=ieor(mran0,ishft(mran0,5))
	mran0=ieor(mran0,ishft(mran0,-13))
	mran0=ieor(mran0,ishft(mran0,6))
	rans=ieor(nran0,rans)+mran0
	harvest=amm*merge(rans,not(rans), rans<0 )
END SUBROUTINE ran1_d

SUBROUTINE ran1_v(harvest)
	USE nrtype
	USE ran_state, ONLY: K4B,amm,lenran,ran_init, &
		iran,jran,kran,nran,mran,ranv
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(OUT) :: harvest
	INTEGER(K4B) :: n
	n=size(harvest)
	if (lenran < n+1) call ran_init(n+1)
	ranv(1:n)=iran(1:n)-kran(1:n)
	where (ranv(1:n) < 0) ranv(1:n)=ranv(1:n)+2147483579_k4b
	iran(1:n)=jran(1:n)
	jran(1:n)=kran(1:n)
	kran(1:n)=ranv(1:n)
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),13))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),-17))
	nran(1:n)=ieor(nran(1:n),ishft(nran(1:n),5))
	where (nran(1:n) == 1) nran(1:n)=270369_k4b
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),5))
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),-13))
	mran(1:n)=ieor(mran(1:n),ishft(mran(1:n),6))
	ranv(1:n)=ieor(nran(1:n),ranv(1:n))+mran(1:n)
	harvest=amm*merge(ranv(1:n),not(ranv(1:n)), ranv(1:n)<0 )
END SUBROUTINE ran1_v

SUBROUTINE gauher(x,w)
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror
	IMPLICIT NONE
	REAL(DP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(DP), PARAMETER :: EPS=3.0e-13_dp,PIM4=0.7511255444649425_dp
	INTEGER(I4B) :: its,j,m,n
	INTEGER(I4B), PARAMETER :: MAXIT=10
	REAL(DP) :: anu
	REAL(DP), PARAMETER :: C1=9.084064e-01_dp,C2=5.214976e-02_dp,&
		C3=2.579930e-03_dp,C4=3.986126e-03_dp
	REAL(DP), DIMENSION((size(x)+1)/2) :: rhs,r2,r3,theta
	REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
	LOGICAL(LGT), DIMENSION((size(x)+1)/2) :: unfinished
	n=assert_eq(size(x),size(w),'gauher')
	m=(n+1)/2
	anu=2.0_dp*n+1.0_dp
	rhs=arth(3,4,m)*PI/anu
	r3=rhs**(1.0_dp/3.0_dp)
	r2=r3**2
	theta=r3*(C1+r2*(C2+r2*(C3+r2*C4)))
	z=sqrt(anu)*cos(theta)
	unfinished=.true.
	do its=1,MAXIT
		where (unfinished)
			p1=PIM4
			p2=0.0
		end where
		do j=1,n
			where (unfinished)
				p3=p2
				p2=p1
				p1=z*sqrt(2.0_dp/j)*p2-sqrt(real(j-1,dp)/real(j,dp))*p3
			end where
		end do
		where (unfinished)
			pp=sqrt(2.0_dp*n)*p2
			z1=z
			z=z1-p1/pp
			unfinished=(abs(z-z1) > EPS)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == MAXIT+1) call nrerror('too many iterations in gauher')
	x(1:m)=-z
	x(n:n-m+1:-1)=z
	w(1:m)=2.0_dp/pp**2
	w(n:n-m+1:-1)=w(1:m)
END SUBROUTINE gauher

SUBROUTINE gauleg(x1,x2,x,w)
	! x1	= lower limit of integration
	! x2	= upper limit of integration
	USE nrtype; USE nrutil, ONLY : arth,assert_eq,nrerror
	IMPLICIT NONE
	REAL(DP), INTENT(IN) :: x1,x2
	REAL(DP), DIMENSION(:), INTENT(OUT) :: x,w
	REAL(DP), PARAMETER :: EPS=3.0e-14_dp
	INTEGER(I4B) :: its,j,m,n
	INTEGER(I4B), PARAMETER :: MAXIT=10
	REAL(DP) :: xl,xm
	REAL(DP), DIMENSION((size(x)+1)/2) :: p1,p2,p3,pp,z,z1
	LOGICAL(LGT), DIMENSION((size(x)+1)/2) :: unfinished
	n=assert_eq(size(x),size(w),'gauleg')
	m=(n+1)/2
	xm=0.5_dp*(x2+x1)
	xl=0.5_dp*(x2-x1)
	z=cos(PI_D*(arth(1,1,m)-0.25_dp)/(n+0.5_dp))
	unfinished=.true.
	do its=1,MAXIT
		where (unfinished)
			p1=1.0
			p2=0.0
		end where
		do j=1,n
			where (unfinished)
				p3=p2
				p2=p1
				p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
			end where
		end do
		where (unfinished)
			pp=n*(z*p1-p2)/(z*z-1.0_dp)
			z1=z
			z=z1-p1/pp
			unfinished=(abs(z-z1) > EPS)
		end where
		if (.not. any(unfinished)) exit
	end do
	if (its == MAXIT+1) call nrerror('too many iterations in gauleg')
	x(1:m)=xm-xl*z
	x(n:n-m+1:-1)=xm+xl*z
	w(1:m)=2.0_dp*xl/((1.0_dp-z**2)*pp**2)
	w(n:n-m+1:-1)=w(1:m)
END SUBROUTINE gauleg

! Compute cholesky decomposition of a
SUBROUTINE choldc(a)
USE nrtype; USE nrutil, ONLY : nrerror
IMPLICIT NONE
REAL(DP), DIMENSION(:,:), INTENT(INOUT) :: a
INTEGER(I4B) :: i,n
real(dp), ALLOCATABLE ::  p(:)
REAL(DP) :: summ
n=size(a,1)
allocate(p(n))
do i=1,n
  summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
  if (summ <= 0.0) call nrerror('choldc failed')
  p(i)=sqrt(summ)
  a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
  a(i,i)=p(i)
  if (i<n) then
    a(i,i+1:n)=0.d0
  end if
end do
deallocate(p)
END SUBROUTINE choldc

! Given lower triangular matrix a and vector b, compute x to solve
! a*x=b
SUBROUTINE cholsl(a,b,x)
  USE nrtype; USE nrutil, ONLY : assert_eq
  IMPLICIT NONE
  REAL(DP), DIMENSION(:,:), INTENT(IN) :: a
  REAL(DP), DIMENSION(:), INTENT(IN) :: b
  REAL(DP), DIMENSION(:), INTENT(INOUT) :: x
  INTEGER(I4B) :: i,n
  n=assert_eq((/size(a,1),size(a,2),size(b),size(x)/),'cholsl')
  do i=1,n
	x(i)=(b(i)-dot_product(a(i,1:i-1),x(1:i-1)))/a(i,i)
  end do
 ! do i=n,1,-1
! 	x(i)=(x(i)-dot_product(a(i+1:n,i),x(i+1:n)))/a(i,i)
!  end do
END SUBROUTINE cholsl

! Given a, the lower triangular cholesky decomposition of a
! This subroutine computes the inva=inv(a)
subroutine invchol(a,inva)
  use nrtype
  use nr, only : cholsl
  implicit none
  real(dp), intent(in) :: a(:,:)
  real(dp), intent(out) :: inva(:,:)
  real(dp), allocatable :: eye(:)
  integer(i4b) i,n

  n=size(a,1)
  allocate(eye(n))
  inva=0.d0

  do i=1,n
    eye=0.d0
    eye(i)=1.d0
    call cholsl(a,eye,inva(:,i))
  end do
  deallocate(eye)
end subroutine invchol 

subroutine mnewt(usrfun,ntrial,x,n,tolx,tolf)
  use nrtype
  use nr, only : ludcmp,lubksb
  implicit none
  interface
    subroutine usrfun(x,n,np,fvec,fjac)
      use nrtype
      implicit none
      real(dp),     intent(in)  :: x(:)
      integer(i4b), intent(in)  :: n,np
      real(dp),     intent(out) :: fvec(:),fjac(:,:)
    end subroutine usrfun
  end interface
  INTEGER(i4b) n,ntrial,NP
  REAL(dp) tolf,tolx,x(n)
  PARAMETER (NP=15)
!    USES lubksb,ludcmp,usrfun
  INTEGER(i4b) i,k,indx(NP)
  REAL(dp) d,errf,errx,fjac(NP,NP),fvec(NP),p(NP)
  
  do k=1,ntrial
    call usrfun(x,n,NP,fvec,fjac)
    errf = sum(abs(fvec))
    if (errf .le. tolf) return
    p = fvec
    call ludcmp(fjac,indx,d)
    call lubksb(fjac,indx,p)
    errx = sum(abs(p))
    x    = x+p
    if(errx.le.tolx) return
  end do
end subroutine mnewt

subroutine lnsrch(n,xold,fold,g,p,x,f,fvec,func,iuser,ruser,stpmax,MaxIter,ifail)
  use nrtype
  implicit none
  integer(i4b), intent(in)    :: n
  real(dp),     intent(in)    :: xold(:)
  real(dp),     intent(in)    :: fold
  real(dp),     intent(in)    :: g(:)
  real(dp),     intent(inout) :: p(:)
  real(dp),     intent(out)   :: x(:)
  real(dp),     intent(out)   :: f
  real(dp),     intent(out)   :: fvec(:)
  integer(i4b), intent(inout) :: iuser(:)
  real(dp),     intent(inout) :: ruser(:)
  real(dp),     intent(in)    :: stpmax
  integer(i4b), intent(in)    :: MaxIter
  integer(i4b), intent(out)   ::  ifail
  interface
    subroutine func(nx,x,f,df,iuser,ruser,iflag)
      use nrtype
      integer(i4b), intent(in)    :: nx
      real(dp),     intent(in)    :: x(:)
      real(dp),     intent(inout) :: f(:)
      real(dp),     intent(inout) :: df(:,:)
      integer(i4b), intent(inout) :: iuser(:)
      real(dp),     intent(inout) :: ruser(:)
      integer(i4b), intent(inout) :: iflag
    end subroutine func
  end interface

  real(dp), allocatable  :: TempDF(:,:)
  integer(i4b)           :: iflag
 
  real(dp)       :: alf,tolx
  parameter (alf=1.e-4,tolx=1.e-7)

  INTEGER(i4b)   :: i1
  REAL(dp)       :: a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope
  real(dp)       :: total,test,tmplam
 
  allocate(TempDF(n,n))
  iflag=1  ! func only computes f, df is an empty placeholder

  ifail = 0

  ! total = length of p
  total = sqrt(sum(p*p))

  if (total .gt. stpmax) then
    p = p*stpmax/total 
  endif
 
  ! slope = dot-product of g and p
  slope = dot_product(g,p) 
 
  if (slope .ge. 0.0d0) then
    print *, 'roundoff problem in lnsrch'
  end if
 
  test = maxval( abs(p)/max(abs(xold),1.0d0) )
 
  alamin=TOLX/test
  alam=1.0d0
  do i1=1,MaxIter
    x = xold + alam*p
    call func(n,x,fvec,TempDF,iuser,ruser,iflag)
    f = 0.5d0*dot_product(fvec,fvec)

    if (alam .lt. alamin) then
      x = xold
      ifail = 1
      exit
    else if (f .le. fold+ALF*alam*slope) then
      exit
    else
      if (alam .eq. 1.0d0) then
        tmplam=-slope/(2.0d0*(f-fold-slope))
      else
        rhs1=f-fold-alam*slope
        rhs2=f2-fold-alam2*slope
        a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
        b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
        if (a .eq. 0.0d0) then
          tmplam=-slope/(2.0d0*b)
        else
          disc=b*b-3.0d0*a*slope
          if (disc .lt. 0.0d0) then
            tmplam=.5d0*alam
          else if (b .le. 0.d0) then
            tmplam=(-b+sqrt(disc))/(3.0d0*a)
          else
            tmplam=-slope/(b+sqrt(disc))
          endif
        endif
        if (tmplam .gt. .5d0*alam) tmplam=.5d0*alam
      endif
    endif
    alam2=alam
    f2=f
    alam=max(tmplam,.1d0*alam)
  end do

  if (i1==MaxIter) then
    ifail=1
  end if
  deallocate(TempDF)
end subroutine lnsrch

subroutine newt(func,x,n,iuser,ruser,fvec,ifail)
  use nrtype
  use nr, only : ludcmp,lubksb,lnsrch
  implicit none
  real(dp),     intent(inout) :: x(:)
  integer(i4b), intent(in)    :: n
  integer(i4b), intent(inout) :: iuser(:)
  real(dp),     intent(inout) :: ruser(:)
  real(dp),     intent(out)   :: fvec(:)
  integer(i4b), intent(out)   :: ifail
  interface
    subroutine func(n,x,f,df,iuser,ruser,iflag)
      use nrtype
      integer(i4b), intent(in)    :: n
      real(dp),     intent(in)    :: x(:)
      real(dp),     intent(inout) :: f(:)
      real(dp),     intent(inout) :: df(:,:)
      integer(i4b), intent(inout) :: iuser(:)
      real(dp),     intent(inout) :: ruser(:)
      integer(i4b), intent(inout) :: iflag
    end subroutine func
  end interface
  
  INTEGER(i4b)    :: MAXITS
  integer(i4b)    :: MaxIter_line_search
  REAL(dp)        :: TOLF,TOLMIN,TOLX,STPMX
  PARAMETER (MAXITS=200,TOLF=1.e-4,TOLMIN=1.e-6,TOLX=1.e-7,stpmx=100.0d0)
  parameter (MaxIter_line_search=100)
  
  INTEGER(i4b)              ::  i,its,j,iflag
  integer(i4b), allocatable :: indx(:)
  REAL(dp)                  ::  d,den,f,fold,stpmax,total,temp,testx,testf
  real(dp),     allocatable :: fjac(:,:),g(:),p(:),xold(:)
 
  allocate(indx(n))
  allocate(fjac(n,n),g(n),p(n),xold(n))

  iflag=1  ! iflag=1 compute fvec only 
  call func(n,x,fvec,fjac,iuser,ruser,iflag)
  f = 0.5d0 * dot_product(fvec,fvec)

  testf = maxval(abs(fvec))
  
  if (testf .lt. .01d0*TOLF) then
    ifail= 0
    deallocate(indx,fjac,g,p,xold)
    return
  endif

  ! length of x 
  total = sqrt(sum(x*x))

  stpmax=STPMX*max(total,real(n,dp))
  do its=1,MAXITS
    iflag=2
    ! compute jacobian of fvec
    call func(n,x,fvec,fjac,iuser,ruser,iflag) 
   
    g = matmul(transpose(fjac),fvec)
    xold = x 
    fold = f
    p    = -fvec

    call ludcmp(fjac,indx,d)
    call lubksb(fjac,indx,p)
    call lnsrch(n,xold,fold,g,p,x,f,fvec,func,iuser,ruser,stpmax,MaxIter_line_search,ifail)

    testf=maxval(abs(fvec))
    if (testf .lt. TOLF) then
      ifail = 0
      exit
    endif
 
    if (ifail==1) then
      den=max(f,.5d0*real(n,dp))
      testf = maxval(   abs(g)*max(abs(x),1.0d0)/den)

      if (testf .lt. TOLMIN) then
        ifail = 1
      else
        ifail = 0
      endif
      exit
    endif
    
    testx = maxval( abs(x-xold)/max(abs(x),1.0d0))
    if (testx .lt. TOLX) exit
    print *,its,testx,testf
  end do

  deallocate(indx)
  deallocate(fjac,g,p,xold)

  if (its==MaxIts) then
    ifail = 1
  end if

end subroutine newt


