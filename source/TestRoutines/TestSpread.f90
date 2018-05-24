program TestSpread
  implicit none
  real(8), allocatable :: e(:,:),mu(:,:),temp_mu1(:),temp_mu2(:,:)
  integer(4), allocatable :: month(:)
  integer(4)    :: k,n,i1,i2

  n = 20
  k = 5
  
  allocate(e(n,k))
  allocate(temp_mu1(k),temp_mu2(n,k))
  allocate(mu(k,12))
  allocate(month(n))

  e = 0.0d0
  mu = 0.0d0
  month = 0

  do i1=1,k
  do i2=1,12
    mu(i1,i2) = real(i1,8)+real(i2,8)
  end do
  end do

  do i1=1,n
    month(i1) = i1 / (n/12+1) +1
  end do

  do i1=1,12
    temp_mu1 = mu(:,i1)
    temp_mu2 = spread(mu(:,i1),1,n)
    e = e + merge(spread(mu(:,i1),1,n),0.0d0,spread(month,2,12)==i1)
  end do

  print *,e
end program
