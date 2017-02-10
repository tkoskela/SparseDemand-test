program TestFormat
  implicit none
  integer i1
  i1 = 5

  write(6,110) i1
  write(6,111) i1
  i1 = 12
  print 110,i1
  print 111,i1
110 format(i2.2)
111 format(i2)
 
end program
