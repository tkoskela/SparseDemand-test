program TestLoadData
  implicit none
  integer(4)         :: J,N,DataUnit,i1
  character(len=100) :: DataFile

  integer(4), allocatable :: HHID(:),date1(:),shopid(:),nNonZero(:),day(:)
  real(8),    allocatable :: qp(:,:),err(:,:)
  

  J = 24
  N = 5
  allocate(HHID(N))
  allocate(date1(N))
  allocate(shopid(N))
  allocate(qp(N,2*J))
  allocate(nNonZero(N))
  allocate(day(N))
  allocate(err(N,J))

  DataUnit = 20
  DataFile = 'TestData.csv'
  open(unit = DataUnit, &
       file = DataFile, &
       action = 'read')
 
  do i1=1,N
   read(DataUnit,389) HHID(i1),date1(i1),shopid(i1),  &
                      qp(i1,:),nNonZero(i1),day(i1),err(i1,:)
  end do
389 format(3i10,<2*J>d25.16,2i10,<J>d25.16)

  close(DataUnit)
  deallocate(HHID,date1,shopid)
  deallocate(qp,err)
  deallocate(nNonzero,day) 

end program TestLoadData
