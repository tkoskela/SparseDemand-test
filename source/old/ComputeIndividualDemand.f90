subroutine ComputeIndividualDemand(SimData1)
  use nrtype
  use GlobalModule, only : DataStructure,parms
  use OutputModule, only : MakeFullFileName
  implicit none
  type(DataStructure), intent(inout) :: SimData1

  integer(i4b)              :: nHH,np,SimUnit
  integer(i4b)              :: i1,i2,ihh
  real(dp), allocatable     :: p(:)
  character(len=99)         :: SimFile

  nHH = SimData1%n
  np  = 30
  allocate(p(np))
  p = 0.0d0

  SimUnit = 45
  SimFile = MakeFullFileName('IndividualDemand.txt')
  open(Unit= SimUNit, &
       File= SimFile, &
       action = 'write')


  do i1=1,parms%J
    p = SimData1%p(i1,1) * (0.5d0+(2.0d0-0.5d0)*real((/0:np-1/),dp)/real(np-1,dp))
    do i2=1,np
      SimData1%p(i1,:) = p(i2)
      call ComputeDemand(SimData1)

      do ihh=1,nHH
        ! iHH  i1 i2  eta(1)  eta(2)  nNonZero i1-i5 q1-q5 p1-p24
        write(SimUnit,30) ihh,i1,i2,SimData1%eta(:,ihh),SimData1%nNonZero(ihh), &
                          SimData1%iNonZero(:,ihh),SimData1%q(:,ihh), &
                          SimData1%p(:,ihh)
      end do
    end do
  end do
30 format(3(i2,","), &
          <parms%dim_eta>(g12.4,","), &
          i2,",",              &
          <parms%K>(i3,","),   &
          <parms%K>(g12.4,",") &
          <parms%J>(g12.4,:,","))

  close(SimUnit)

  deallocate(p)
end subroutine ComputeIndividualDemand
