program TestMaxArray
  implicit none
  real(8), allocatable :: x(:,:),grad(:),hess(:,:)
  integer(4)           :: n,k,i1,i2
  n = 20000
  k = 285
  allocate(x(n,k))
  allocate(grad(k))
  allocate(hess(k,k))

  call random_number(x)
  do i1=1,k
    grad(i1) = sum(x(:,i1))/real(n,8)
    do i2=1,i1
      hess(i1,i2) = sum((x(:,i1)-grad(i1))*(x(:,i2)-grad(i2)))/real(n,8)
      hess(i2,i1) = hess(i1,i2)
    end do
  end do   
end program
