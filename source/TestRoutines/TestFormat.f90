program TestFormat
  implicit none
  integer i1,i2,n
  real    x1,x2
  i1 = 5
  i2 = 6
  x1 = 3.2d0
  x2 = 4.3d0
  n = 1

  write(6,110) i1
  write(6,111) i1
  i1 = 12
  print 110,i1
  print 111,i1
110 format(i2.2)
111 format(i2)

  write(6,112) i1,i2,i1,x1,x2
112 format(2(i2,","),i2,",", &
           <n+1>(g12.4,:,","))
end program
