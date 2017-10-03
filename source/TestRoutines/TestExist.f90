program TestExist
  implicit none
  logical exflag1,exflag2

  open(unit=10,file='test1.txt',action='write')
  write(10,*) 'Hello'
  close(10)

  inquire(file='test1.txt',exist=exflag1)
  inquire(file='test2.txt',exist=exflag2)
  print *,'test1.txt',exflag1
  print *,'test2.txt',exflag2
end program
