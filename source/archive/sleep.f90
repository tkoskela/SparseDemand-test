program sleep0
  use IFPORT
  implicit none

  integer :: i1
 
  i1=0 
  do while (i1==0)
    call sleep(5)
  end do


end program
