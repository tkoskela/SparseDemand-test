program test_write
  implicit none
  integer  :: i1,i2,i3,n1,n2
  real(8), allocatable :: x(:,:)
  integer  :: ix1,ix2
  integer  :: unit1
  character(len=13) :: filename

  unit1 = 10
  filename = "parms.txt"
 
  i1 = 4
  i2 = 10
  i3 = 20
  n1 = 3
  n2 = 3

allocate(x(n1,n2))

  do ix1=1,n1
  do ix2=1,n2
    x(ix1,ix2) = dble(ix1) + dble(ix2)
  end do
  end do

  ! write parameters
  open(unit=unit1,  &
       file=filename, &
       action='write')

  write(unit1,99) i1,i2,i3
99 format(2(i4,','),i4)
  do ix1=1,n1
    write(unit1,100) x(ix1,:)
  end do
100 format(<n2-1>(d25.16,','),d25.16)

  close(unit1)

deallocate(x)

end program
