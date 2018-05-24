program TestDate
  implicit none
  integer(4) :: date(5),year(5),month(5),day(5)

  date  = (/20080321,20090401,20100131,20080612,20081205/)
  year  = date/10000
  month = (date - 10000*year)/100
  day       = (date-10000*year-100*month)

  print *,date
  print *,year
  print *,month
  print *,day
end program
