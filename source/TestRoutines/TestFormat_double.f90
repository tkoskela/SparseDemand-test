program TestFormat_double
implicit none
real(8), allocatable  :: x(:),y(:)
integer(4)            :: i1,n,file_unit
character(len=20)     :: file_name

n=6
file_unit=20
file_name = 'double.txt'
allocate(x(n))
allocate(y(n))

call random_number(x)
x = x-1.0d0
do i1=2,n
  x(i1) = x(1)
end do

write(6,10) x(1)
write(6,11) x(1)
write(6,12) x(1)
write(6,13) x(1)
write(6,14) x(1)
write(6,15) x(1)

open(unit=file_unit, &
     file=file_name, &
     action='write')
write(file_unit,10) x(1)
write(file_unit,11) x(1)
write(file_unit,12) x(1)
write(file_unit,13) x(1)
write(file_unit,14) x(1)
write(file_unit,15) x(1)
close(file_unit)
open(unit=file_unit, &
     file=file_name, &
     action='read')
read(file_unit,10) y(1)
read(file_unit,11) y(2)
read(file_unit,12) y(3)
read(file_unit,13) y(4)
read(file_unit,14) y(5)
read(file_unit,15) y(6)
close(file_unit)


10 format(g30.16)
11 format(g30.17)
12 format(g30.18)
13 format(g30.19)
14 format(g30.20)
15 format(g30.21)

do i1=1,n
  write(6,'(2g30.21)') x(i1),y(i1)
end do

deallocate(x,y)
end program TestFormat_double
