program TestReadCSV
! test reading of CSV file
! row 1:   integers
! row 2:   characters
! row 3:   characters
! row 4:   real numbers
! ...
! row 3+J: real numbers
  implicit none
  character(len=20)   :: csvfile
  integer             :: csvunit
  character(len=1024) :: buffer1,buffer2

  integer,            allocatable :: itemp(:),taxid(:)
  character(len=100), allocatable :: taxlabel(:)
  character(len=10),  allocatable :: taxtype(:)
  real(8),            allocatable :: tax(:,:)
  integer  :: i1,n1,icomma,J

  csvfile = "TestReadCSV.csv"
  J = 2

  open(unit = csvunit, &
       file = csvfile, &
       action = "read")

  ! read in row 1 (integers) and determine number of columns
  read(csvunit,'(a)') buffer1

  i1 = 1
  icomma = 1
  allocate(itemp(20))
  itemp = 0
  do while (icomma>0)
    icomma = index(buffer1,",")
    if (icomma>0) then
      read(buffer1(1:(icomma-1)),'(i2)') itemp(i1) 
      i1=i1+1
      buffer1=buffer1((icomma+1):len_trim(buffer1))
    else
      read(buffer1,'(i2)') itemp(i1) 
    end if
  end do
  n1 = i1
  allocate(taxid(n1))
  taxid = itemp(1:i1)

  ! Read in rows 2-3
  allocate(taxlabel(n1),taxtype(n1))

  read(csvunit,'(a)') buffer1
  read(csvunit,'(a)') buffer2

  do i1=1,n1
    if (i1<n1) then
      icomma       = index(buffer1,",")
      taxlabel(i1) = buffer1(1:(icomma-1))
      buffer1      = buffer1((icomma+1):len_trim(buffer1))
      icomma       = index(buffer2,",")
      taxtype(i1)  = buffer2(1:(icomma-1))
      buffer2      = buffer2((icomma+1):len_trim(buffer2))
    else
      taxlabel(i1) = buffer1
      taxtype(i1)  = buffer2
    end if
  end do

  ! read in (J x n1) tax rates
  allocate(tax(J,n1))
  tax = 1.0d0
  do i1=1,J
    read(csvunit,69) tax(i1,:)
  end do
 69 format(<n1>g12.4)  

close(csvunit)

print 56, "tax_id",taxid
56 format(a10,<n1>i4) 
print 57, "taxlabel",taxlabel
57 format(a10,<n1>a20)
print 58, "taxtype",taxtype
58 format(a10,<n1>a20)
do i1=1,J
  print 59, "tax",i1,tax(i1,:)
end do
59 format(a10,i2.2,<n1>g12.4)
  deallocate(itemp)
  deallocate(taxid,taxlabel,taxtype,tax)
end program
