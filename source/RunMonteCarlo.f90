subroutine MaxOneAtTime(pid)
  use nrtype
  use GlobalModule, only : MasterID,iFree,SelectFreeType,parms, &
                           ReadWriteParameters,DeallocateIFree

  implicit none
  integer(i4b), intent(in) :: pid
  type(SelectFreeType)     :: iFree0
  integer(i4b)             :: i1,ifail
  real(dp), allocatable    :: x(:),Grad(:),Hess(:,:)
  real(dp)                 :: LValue

  call CopyIFree(iFree,iFree0)

  do i1=1,iFree0%NALL
    call DeallocateIFree(iFree)
    call UpdateIFree(i1,iFree0,iFree)
    if (pid==MasterID) then
      allocate(x(iFree%NALL))
      allocate(Grad(iFree%NALL),Hess(iFree%NALL,iFree%NALL))
      call ComputeInitialGuess(parms,iFree,x)
      call MaximizeLikelihood(x,LValue,Grad,Hess,ifail)
      ! update parms
      call UpdateParms2(x,iFree,parms)
      ! b) save parms to disk
      call ReadWriteParameters(parms,'write')
      deallocate(x,Grad,Hess)
    elseif (pid .ne. MasterID) then
#if USE_MPI>0
      call WorkerTask(parms%model,pid)
#endif
    end if
  end do

  call DeallocateIFree(iFree)
  call CopyIFree(iFree0,iFree)
  call DeallocateIFree(iFree0)

end subroutine MaxOneAtTime

subroutine CopyIFree(iFree,iFreeCopy)
  use GlobalModule, only : SelectFreeType
  implicit none
  type(SelectFreeType), intent(in) :: iFree
  type(SelectFreeType), intent(inout) :: iFreeCopy

  iFreeCopy%ND           = iFree%ND
  iFreeCopy%NBC          = iFree%NBC
  iFreeCopy%NMUE         = iFree%NMUE
  iFreeCopy%NINVCDIAG    = iFree%NINVCDIAG
  iFreeCopy%NINVCOFFDIAG = iFree%NINVCOffDiag

  iFreeCopy%NBD_beta     = iFree%NBD_beta
  iFreeCopy%NBD_CDIAG    = iFree%NBD_CDIAG
  iFreeCopy%NBD_COFFDIAG = iFree%NBD_COFFDIAG

  iFreeCopy%NBC_CDiag     = iFree%NBC_beta
  iFreeCopy%NBC_CDIAG    = iFree%NBC_CDIAG
  iFreeCopy%NBC_COFFDIAG = iFree%NBC_COFFDIAG

  iFreeCopy%NALL         = iFREE%NALL

  iFreeCopy%flagD           = iFree%flagD
  iFreeCopy%flagBC          = iFree%flagBC
  iFreeCopy%flagMUE         = iFree%flagMUE
  iFreeCopy%flagInvCDiag    = iFree%flagInvCDiag
  iFreeCopy%flagInvCOffDiag = iFree%flagInvCOffDiag

  iFreeCopy%flagBC_beta     = iFree%flagBC_beta
  iFreeCopy%flagBC_CDiag    = iFree%flagBC_CDiag
  iFreeCopy%flagBC_COffDiag = iFree%flagBC_COffDiag

  iFreeCopy%flagBD_beta     = iFree%flagBD_beta
  iFreeCopy%flagBD_CDiag    = iFree%flagBD_CDiag
  iFreeCopy%flagBD_COffDiag = iFree%flagBD_COffDiag
  iFreeCopy%OneAtATime      = iFree%OneAtATime

  if (iFree%flagD>0) then
    allocate(iFreeCopy%D(iFree%ND))
    allocate(iFreeCopy%xD(iFree%ND))
    iFreeCopy%D = iFree%D
    iFreeCopy%xD = iFree%xD
  end if

  if (iFree%flagBC>0) then
    allocate(iFreeCopy%BC(iFree%NBC))
    allocate(iFreeCopy%xBC(iFree%NBC))
    iFreeCopy%BC = iFree%BC
    iFreeCopy%xBC = iFree%xBC
  end if

  if (iFree%flagMUE>0) then
    allocate(iFreeCopy%MUE(iFree%NMUE))
    allocate(iFreeCopy%xD(iFree%NMUE))
    iFreeCopy%MUE = iFree%MUE
    iFreeCopy%xMUE = iFree%xMUE
  end if

  if (iFree%flagInvCDiag>0) then
    allocate(iFreeCopy%InvCDiag(iFree%NInvCDiag))
    allocate(iFreeCopy%xInvCDiag(iFree%NInvCDiag))
    iFreeCopy%InvCDiag = iFree%InvCDiag
    iFreeCopy%xInvCDiag = iFree%xInvCDiag
  end if

  if (iFree%flagInvCOffDiag>0) then
    allocate(iFreeCopy%InvCOffDiag(iFree%NInvCOffDiag))
    allocate(iFreeCopy%xInvCOffDiag(iFree%NInvCOffDiag))
    iFreeCopy%InvCOffDiag = iFree%InvCOffDiag
    iFreeCopy%xInvCOffDiag = iFree%xInvCOffDiag
  end if

  if (iFree%flagBC_beta>0) then
    allocate(iFreeCopy%BC_beta(iFree%NBC_beta))
    allocate(iFreeCopy%xBC_beta(iFree%NBC_beta))
    iFreeCopy%BC_beta = iFree%BC_beta
    iFreeCopy%xBC_beta = iFree%xBC_beta
  end if

  if (iFree%flagBC_CDiag>0) then
    allocate(iFreeCopy%BC_CDiag(iFree%NBC_CDiag))
    allocate(iFreeCopy%xBC_CDiag(iFree%NBC_CDiag))
    iFreeCopy%BC_CDiag = iFree%BC_CDiag
    iFreeCopy%xBC_CDiag = iFree%xBC_CDiag
  end if

  if (iFree%flagBC_COffDiag>0) then
    allocate(iFreeCopy%BC_COffDiag(iFree%NBC_COffDiag))
    allocate(iFreeCopy%xBC_COffDiag(iFree%NBC_COffDiag))
    iFreeCopy%BC_COffDiag = iFree%BC_COffDiag
    iFreeCopy%xBC_COffDiag = iFree%xBC_COffDiag
  end if

  if (iFree%flagBD_beta>0) then
    allocate(iFreeCopy%BD_beta(iFree%NBD_beta))
    allocate(iFreeCopy%xBD_beta(iFree%NBD_beta))
    iFreeCopy%BD_beta = iFree%BD_beta
    iFreeCopy%xBD_beta = iFree%xBD_beta
  end if

  if (iFree%flagBD_CDiag>0) then
    allocate(iFreeCopy%BD_CDiag(iFree%NBD_CDiag))
    allocate(iFreeCopy%xBD_CDiag(iFree%NBD_CDiag))
    iFreeCopy%BD_CDiag = iFree%BD_CDiag
    iFreeCopy%xBD_CDiag = iFree%xBD_CDiag
  end if

  if (iFree%flagBD_COffDiag>0) then
    allocate(iFreeCopy%BD_COffDiag(iFree%NBD_COffDiag))
    allocate(iFreeCopy%xBD_COffDiag(iFree%NBD_COffDiag))
    iFreeCopy%BD_COffDiag = iFree%BD_COffDiag
    iFreeCopy%xBD_COffDiag = iFree%xBD_COffDiag
  end if

end subroutine CopyIFree

subroutine UpdateIFree(i1,iFree0,iFreeNew)
  use GlobalModule, only : SelectFreeType
  implicit none
  integer(i4b), intent(in) :: i1
  type(SelectFreeType), intent(in)    :: iFree0
  type(SelectFreeType), intent(inout) :: iFreeNew

  iFreeNew%nall         = 1
  iFreeNew%nD           = 0
  iFreeNew%nBC          = 0
  iFreeNew%nMUE         = 0
  iFreeNew%nInvCDiag    = 0
  iFreeNew%nInvCOffDiag = 0
  iFreeNew%nBC_beta     = 0
  iFreeNew%nBC_CDiag    = 0
  iFreeNew%nBC_COffDiag = 0
  iFreeNew%nBD_beta     = 0
  iFreeNew%nBD_CDiag    = 0
  iFreeNew%nBD_COffDiag = 0

  if (iFree0%nD>0) then
    if (any(iFree0%xD==i1)) then
      allocate(iFreeNew%D(1))
      allocate(IFreeNew%xD(1))
      iFreeNew%xD = 1
      iFreeNew%D  = pack(iFree0%D,iFree0%xD==i1)
      iFreeNew%nD = 1
    end if
  end if

  if (iFree0%nBC>0) then
    if (any(iFree0%xBC==i1)) then
      allocate(iFreeNew%BC(1))
      allocate(IFreeNew%xBC(1))
      iFreeNew%xBC = 1
      iFreeNew%BC  = pack(iFree0%BC,iFree0%xBC==i1)
      iFreeNew%nBC = 1
    end if
  end if

  if (iFree0%nMUE>0) then
    if (any(iFree0%xMUE==i1)) then
      allocate(iFreeNew%MUE(1))
      allocate(IFreeNew%xMUE(1))
      iFreeNew%xMUE = 1
      iFreeNew%MUE  = pack(iFree0%MUE,iFree0%xMUE==i1)
      iFreeNew%nMUE = 1
    end if
  end if

  if (iFree0%nInvCDiag>0) then
    if (any(iFree0%xInvCDiag==i1)) then
      allocate(iFreeNew%InvCDiag(1))
      allocate(IFreeNew%xInvCDiag(1))
      iFreeNew%xInvCDiag = 1
      iFreeNew%InvCDiag  = pack(iFree0%InvCDiag,iFree0%xInvCDiag==i1)
      iFreeNew%nInvCDiag = 1
    end if
  end if

  if (iFree0%nInvCOffDiag>0) then
    if (any(iFree0%xInvCOffDiag==i1)) then
      allocate(iFreeNew%InvCOffDiag(1))
      allocate(IFreeNew%xInvCOffDiag(1))
      iFreeNew%xInvCOffDiag = 1
      iFreeNew%InvCOffDiag  = pack(iFree0%InvCOffDiag,iFree0%xInvCOffDiag==i1)
      iFreeNew%nInvCOffDiag = 1
    end if
  end if

  if (iFree0%nBC_beta>0) then
    if (any(iFree0%xBC_beta==i1)) then
      allocate(iFreeNew%BC_beta(1))
      allocate(IFreeNew%xBC_beta(1))
      iFreeNew%xBC_beta = 1
      iFreeNew%BC_beta  = pack(iFree0%BC_beta,iFree0%xBC_beta==i1)
      iFreeNew%nBC_beta = 1
    end if
  end if
  if (iFree0%nBD_beta>0) then
    if (any(iFree0%xBD_beta==i1)) then
      allocate(iFreeNew%BD_beta(1))
      allocate(IFreeNew%xBD_beta(1))
      iFreeNew%xBD_beta = 1
      iFreeNew%BD_beta  = pack(iFree0%BD_beta,iFree0%xBD_beta==i1)
      iFreeNew%nBD_beta = 1
    end if
  end if

  if (iFree0%nBC_CDiag>0) then
    if (any(iFree0%xBC_CDiag==i1)) then
      allocate(iFreeNew%BC_CDiag(1))
      allocate(IFreeNew%xBC_CDiag(1))
      iFreeNew%xBC_CDiag = 1
      iFreeNew%BC_CDiag  = pack(iFree0%BC_CDiag,iFree0%xBC_CDiag==i1)
      iFreeNew%nBC_CDiag = 1
    end if
  end if

  if (iFree0%nBD_CDiag>0) then
    if (any(iFree0%xBD_CDiag==i1)) then
      allocate(iFreeNew%BD_CDiag(1))
      allocate(IFreeNew%xBD_CDiag(1))
      iFreeNew%xBD_CDiag = 1
      iFreeNew%BD_CDiag  = pack(iFree0%BD_CDiag,iFree0%xBD_CDiag==i1)
      iFreeNew%nBD_CDiag = 1
    end if
  end if

  if (iFree0%nBC_COffDiag>0) then
    if (any(iFree0%xBC_COffDiag==i1)) then
      allocate(iFreeNew%BC_COffDiag(1))
      allocate(IFreeNew%xBC_COffDiag(1))
      iFreeNew%xBC_COffDiag = 1
      iFreeNew%BC_COffDiag  = pack(iFree0%BC_COffDiag,iFree0%xBC_COffDiag==i1)
      iFreeNew%nBC_COffDiag = 1
    end if
  end if

  if (iFree0%nBD_COffDiag>0) then
    if (any(iFree0%xBD_COffDiag==i1)) then
      allocate(iFreeNew%BD_COffDiag(1))
      allocate(IFreeNew%xBD_COffDiag(1))
      iFreeNew%xBD_COffDiag = 1
      iFreeNew%BD_COffDiag  = pack(iFree0%BD_COffDiag,iFree0%xBD_COffDiag==i1)
      iFreeNew%nBD_COffDiag = 1
    end if
  end if

end subroutine UpdateIFree
