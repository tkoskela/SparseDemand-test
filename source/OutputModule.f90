module OutputModule
  use nrtype
  use GlobalModule, only : OutDir
  implicit none
  
  ! Output filenames
  character(len=200)          :: Results_FILE
  character(len=200)          :: BayesResults_FILE
  character(len=200)          :: Results1P_FILE
  character(len=200)          :: Hess_FILE
  character(len=200)          :: GradLHH_FILE
  character(len=200)          :: SaveDataFile_q
  character(len=200)          :: SaveDataFile_p
  character(len=200)          :: SaveDataFile_e
  character(len=200)          :: SaveDataFile_eta
  character(len=200)          :: SaveDataFile_iNonZero
  character(len=200)          :: SaveDataFile_iZero
  character(len=200)          :: SaveDataFile_nNonZero
  character(len=200)          :: SaveDataFile_market
  character(len=200)          :: BIC1_File
  character(len=200)          :: MC_FILE1
  character(len=200)          :: MC_FILE2
  character(len=200)          :: ElasFile
  character(len=200)          :: DemandFile

  ! Output file: unit numbers
  integer(i4b), parameter :: Results_UNIT           = 12
  integer(i4b), parameter :: BayesResults_UNIT      = 13
  integer(i4b), parameter :: Results1P_UNIT         = 14 
  integer(i4b), parameter :: Hess_UNIT              = 15
  integer(i4b), parameter :: SaveData_UNIT_q        = 16
  integer(i4b), parameter :: SaveData_UNIT_p        = 17
  integer(i4b), parameter :: SaveData_UNIT_e        = 18
  integer(i4b), parameter :: SaveData_UNIT_iNonZero = 19
  integer(i4b), parameter :: SaveData_UNIT_iZero    = 20
  integer(i4b), parameter :: SaveData_UNIT_nNonZero = 21
  integer(i4b), parameter :: SaveData_UNIT_market   = 22
  integer(i4b), parameter :: BIC1_UNIT              = 23
  integer(i4b), parameter :: MC_UNIT1               = 24
  integer(i4b), parameter :: SaveData_UNIT_eta      = 25
  integer(i4b), parameter :: MC_UNIT2               = 26
  integer(i4b), parameter :: Elas_UNIT              = 27
  integer(i4b), parameter :: Demand_UNIT            = 28
  integer(i4b), parameter :: GradLHH_unit           = 29
  integer(i4b), parameter :: taxresults_unit        = 30

contains

subroutine DefineFileNames(pid)
  implicit none
  integer(i4b), intent(in) :: pid
  character(len=30)        :: TempFileName

!  OutDir        = 'output'
  Results_File          = MakeFullFileName('results.txt')
  BayesResults_File     = MakeFullFileName('BayesResults.txt')
  Hess_FILE             = MakeFullFileName('hess.txt')
  GradLHH_FILE          = MakeFullFileName('gradLHH.txt')
  SaveDataFile_q        = MakeFullFileName('q.csv')
  SaveDataFile_p        = MakeFullFileName('p.csv')
  SaveDataFile_e        = MakeFullFileName('e.csv')
  SaveDataFile_eta      = MakeFullFileName('eta.csv')
  SaveDataFile_iNonZero = MakeFullFileName('iNonZero.csv')
  SaveDataFile_iZero    = MakeFullFileName('iZero.csv')
  SaveDataFile_nNonZero = MakeFullFileName('nNonZero.csv')
  SaveDataFile_market   = MakeFullFileName('market.csv')
  ElasFile              = MakeFullFileName('elas.csv')

  write(TempFileName,'(A9,i4.4,a4)') 'MCResults',pid,'.txt'
  MC_FILE1      = MakeFullFileName(TempFileName)
  write(TempFileName,'(A8,i4.4,a4)') 'MCLambda',pid,'.txt'
  MC_FILE2      = MakeFullFileName(TempFileName)

end subroutine DefineFileNames

subroutine DefinePenaltyFileNames(iter)
  implicit none
  integer(i4b), intent(in) :: iter
  character(len=50)        :: results1,results2
  write(results1,'(A10,I2.2,A4)') 'Results1_p',iter,'.txt'
  Results1P_File = MakeFullFileName(results1)
  BIC1_FILE      = MakeFullFileName('BIC1.txt')
end subroutine DefinePenaltyFileNames

!------------------------------------------------------------------------------
!
! subroutine SaveOutputs(xFree,LValue,Grad,Hess)
!
!  Hess   : estimate of Hessian
!------------------------------------------------------------------------------
subroutine SaveOutputs(xFree,LValue,Grad,Hess,Stats)
  use NewTools,     only : MatrixInverse
  use GlobalModule, only : iFree,           &
                           HHData,          &
                           ControlOptions,  &
                           ResultStructure, &
                           Penalty

  implicit none
  real(dp),     intent(in)  :: xFree(:)
  real(dp),     intent(in)  :: LValue
  real(dp),     intent(in)  :: Grad(:)
  real(dp),     intent(in)  :: Hess(:,:)
  type(ResultStructure), intent(in) :: stats
  
  ! variables used to compute standard errors  
  integer(i4b)          :: i1
  real(dp), allocatable :: StandardErrors(:)
  real(dp), allocatable :: InvHess(:,:)
  
  integer(i4b)          :: n 
  character(20)         :: cTemp

  ! Define filenames
  !call DefineFileNames

  n = size(xFree,1)
  allocate(InvHess(n,n))
  allocate(StandardErrors(n))
 
  call MatrixInverse(Hess,InvHess,'Symmetric')
  do i1=1,n
    StandardErrors(i1) = dsqrt(InvHess(i1,i1))
  end do
  
  !--------------------------------------------------------------------------
  ! Write parameters and std. err.
  !--------------------------------------------------------------------------
  open(UNIT = Results_UNIT, &
       FILE = Results_FILE, &
       ACTION = 'WRITE')

  write(Results_UNIT,2) 'Variable','Coef','Gradient','s.e.'

  ! write parms%B_D
  if (iFree%flagD>0) then
    do i1=1,iFree%nD
      ! parms%D(iFree%D) = xFree(iFree%xD)
      write(Results_UNIT,1) HHData%ColumnLabels(iFree%D(i1)), &
                            xFree(iFree%xD(i1)),  &
                            Grad(iFree%xD(i1)),   &
                            StandardErrors(iFree%xD(i1))
      end do
  end if

  ! write parms%CSig
  if (iFree%flagBC>0) then
    do i1=1,iFree%nBC
      ! parms%CSig(iFree%BC) = xFree(iFree%xBC)
      write(cTemp,'(a4,i2.2)') 'phi_',i1
      write(Results_UNIT,1) cTemp,                &
                            xFree(iFree%xBC(i1)),  &
                            Grad(iFree%xBC(i1)),   &
                            StandardErrors(iFree%xBC(i1))
    end do
  end if

  ! write parms%MUE
  if (iFree%flagMUE>0) then
    do i1=1,iFree%nMUE
      ! parms%MUE(iFree%mUE) = xFree(iFree%xMUE)
      write(cTemp,'(a4,i2.2)') 'MUE_',i1
      write(Results_UNIT,1) cTemp,                &
                            xFree(iFree%xMUE(i1)),  &
                            Grad(iFree%XMUE(i1)),   &
                            StandardErrors(iFree%XMUE(i1))
    end do
  end if

  ! write parms%INVCDiag
  if (iFree%flagInvCDiag>0) then
    do i1=1,iFree%nINVCDiag
      ! parms%INVCDiag(iFree%INVCDiag) = xFree(iFree%xInvCDiag)
      write(cTemp,'(a7,i2.2)') 'INVCDiag',i1
      write(Results_UNIT,1) cTemp,                &
                            xFree(iFree%xInvCDiag(i1)),  &
                            Grad(iFree%xInvCDiag(i1)),   &
                            StandardErrors(iFree%xInvCDiag(i1))
    end do
  end if

! write parms%INVCOffDiag
  if (iFree%flagInvCOffDiag>0) then
    do i1=1,iFree%nINVCOffDiag
      ! parms%INVCOffDiag(iFree%INVCOffDiag) = xFree(iFree%xInvCOffDiag)
      write(cTemp,'(a9,i2.2)') 'INVCOffDiag',i1
      write(Results_UNIT,1) cTemp,                &
                            xFree(iFree%xInvCOffDiag(i1)),  &
                            Grad(iFree%xInvCOffDiag(i1)),   &
                            StandardErrors(iFree%xInvCOffDiag(i1))
    end do
  end if

  if (iFree%flagBC_beta>0) then
    do i1=1,iFree%nBC_beta
      ! parms%BC_beta(iFree%BC_beta) = xFree(iFree%xBC_beta)
      write(cTemp,'(a9,i2.2)') 'BC_beta',i1
      write(Results_UNIT,1) cTemp,                &
                            xFree(iFree%xBC_beta(i1)),  &
                            Grad(iFree%xBC_beta(i1)),   &
                            StandardErrors(iFree%xBC_beta(i1))
    end do
  end if
  if (iFree%flagBC_CDiag>0) then
    do i1=1,iFree%nBC_CDiag
      ! parms%BC_beta(iFree%BC_CDiag) = xFree(iFree%xBC_CDiag)
      write(cTemp,'(a9,i2.2)') 'BC_CDiag',i1
      write(Results_UNIT,1) cTemp,                &
                            xFree(iFree%xBC_CDiag(i1)),  &
                            Grad(iFree%xBC_CDiag(i1)),   &
                            StandardErrors(iFree%xBC_CDiag(i1))
    end do
  end if

  if (iFree%flagBC_COffDiag>0) then
    do i1=1,iFree%nBC_COffDiag
      ! parms%BC_COffDiag(iFree%BC_COffDiag) = xFree(iFree%xBC_COffDiag)
      write(cTemp,'(a9,i2.2)') 'BC_COffDiag',i1
      write(Results_UNIT,1) cTemp,                &
                            xFree(iFree%xBC_COffDiag(i1)),  &
                            Grad(iFree%xBC_COffDiag(i1)),   &
                            StandardErrors(iFree%xBC_COffDiag(i1))
    end do
  end if

  if (iFree%flagBD_beta>0) then
    do i1=1,iFree%nBD_beta
      ! parms%BD_beta(iFree%BD_beta) = xFree(iFree%xBD_beta)
      write(cTemp,'(a9,i2.2)') 'BD_beta',i1
      write(Results_UNIT,1) cTemp,                &
                            xFree(iFree%xBD_beta(i1)),  &
                            Grad(iFree%xBD_beta(i1)),   &
                            StandardErrors(iFree%xBD_beta(i1))
    end do
  end if

  if (iFree%flagBD_CDiag>0) then
    do i1=1,iFree%nBD_CDiag
      ! parms%BD_CDiag(iFree%BD_CDiag) = xFree(iFree%xBD_CDiag)
      write(cTemp,'(a9,i2.2)') 'BD_CDiag',i1
      write(Results_UNIT,1) cTemp,                &
                            xFree(iFree%xBD_CDiag(i1)),  &
                            Grad(iFree%xBD_CDiag(i1)),   &
                            StandardErrors(iFree%xBD_CDiag(i1))
    end do
  end if

  if (iFree%flagBD_COffDiag>0) then
    do i1=1,iFree%nBD_COffDiag
      ! parms%BD_COffDiag(iFree%BD_COffDiag) = xFree(iFree%xBD_COffDiag)
      write(cTemp,'(a9,i2.2)') 'BD_COffDiag',i1
      write(Results_UNIT,1) cTemp,                &
                            xFree(iFree%xBD_COffDiag(i1)),  &
                            Grad(iFree%xBD_COffDiag(i1)),   &
                            StandardErrors(iFree%xBD_COffDiag(i1))
    end do
  end if

  1 format( a20,3es25.16)
  2 format( a20,3a25   )
  ! write stats: (likelihood,nobs,maxgrad,MinEigHess,BIC,nNonZero)
  write(Results_UNIT,'(A12,g11.4)') 'Likelihood',LValue
  write(Results_UNIT,'(A12,I11)')   'nobs',stats%N
  write(Results_UNIT,'(A12,I11)')   'No. Parms',stats%nNonZero
  write(Results_UNIT,'(A12,g11.4)') 'BIC',stats%BIC
  write(Results_UNIT,'(A12,g11.4)') 'Max. Grad.',stats%MaxGrad
  write(Results_UNIT,'(A12,g11.4)') 'Min. Eig.',stats%MinEigHess
  if (Penalty%method>0) then
    write(Results_UNIT,'(A12,g11.4)') 'Penalty',Penalty%lambda
  end if
  close(UNIT = Results_UNIT)
  
  !--------------------------------------------------------------------------
  ! Write Hessian for model 1
  !--------------------------------------------------------------------------
  open(UNIT = Hess_UNIT, &
       FILE = Hess_FILE, &
       ACTION = 'WRITE')
  do i1=1,n
    write(Hess_UNIT,3)  Hess(i1,:)
  end do
  3 format(<n>es25.16)
  close(UNIT = Hess_UNIT) 
  deallocate(InvHess,StandardErrors)
end subroutine SaveOutputs

subroutine WriteHess(hess)
  implicit none
  real(dp), intent(in) :: hess(:,:)
  integer(i4b) :: i1,n

  n = size(hess,1)

  open(UNIT = Hess_UNIT, &
       FILE = Hess_FILE, &
       ACTION = 'WRITE')
  do i1=1,n
    write(Hess_UNIT,3)  Hess(i1,:)
  end do
  3 format(<n>es25.16)
  close(UNIT = Hess_UNIT) 
  
end subroutine WriteHess

subroutine ReadWriteGradLHH(GradLHH,action)
  implicit none
  real(dp), intent(inout)      :: GradLHH(:,:)
  character(len=*), intent(in) :: action

  integer(i4b)                 :: n,nx,i1
  logical                      :: existflag

  n = size(GradLHH,1)
  nx = size(GradLHH,2)

  select case (action)
    case ('read')
      GradLHH = 0.0d0
      inquire(file=GradLHH_file,exist=existflag)
      if (existflag) then
        open(unit = GradLHH_unit, &
             file = GradLHH_file, &
             action = 'read')
        do i1=1,n
          read(GradLHH_unit,329) GradLHH(i1,:)
        end do
        close(GradLHH_unit)
      end if
    case('write')
      open(unit = GradLHH_unit, &
           file = GradLHH_file, &
           action = 'write')
      do i1=1,n
        write(GradLHH_unit,329) GradLHH(i1,:)
      end do
      close(GradLHH_unit)
  end select
329 format(<nx>g25.16)

end subroutine ReadWriteGradLHH
!------------------------------------------------------------------------------
!
! subroutine SavePenaltyOutputs(iter,model,xFree,LValue,Grad,Hess)
!
!------------------------------------------------------------------------------
subroutine SavePenaltyOutputs(iter,model,xFree,LValue,Grad,Hess,Stats)
  use NewTools,     only : MatrixInverse
  use GlobalModule, only : iFree,HHData,   &
                           ControlOptions, &
                           ResultStructure, &
                           Penalty

  implicit none
  integer(i4b), intent(in)  :: iter
  integer(i4b), intent(in)  :: model
  real(dp),     intent(in)  :: xFree(:)
  real(dp),     intent(in)  :: LValue
  real(dp),     intent(in)  :: Grad(:)
  real(dp),     intent(in)  :: Hess(:,:)
  type(ResultStructure), intent(in) :: stats
  
  ! variables used to compute standard errors  
  integer(i4b)          :: i1
  real(dp), allocatable :: StandardErrors(:)
  real(dp), allocatable :: InvHess(:,:)
  
  integer(i4b)          :: n 
  character(20)         :: cTemp
  ! Define filenames
  call DefinePenaltyFileNames(iter)

  n = size(xFree,1)
  allocate(InvHess(n,n))
  allocate(StandardErrors(n))

  call MatrixInverse(Hess,InvHess,'Symmetric')
  do i1=1,n
    StandardErrors(i1) = dsqrt(InvHess(i1,i1))
  end do

  !----------------------------------------------------------------------------
  !  Write Model 1 output
  !---------------------------------------------------------------------------- 
  if (model==1) then

    !--------------------------------------------------------------------------
    ! Write parameters and std. err.
    !--------------------------------------------------------------------------
    open(UNIT = Results1P_UNIT, &
         FILE = Results1P_FILE, &
         ACTION = 'WRITE')

    write(Results1P_UNIT,2) 'Variable','Coef','Gradient','s.e.' 
    do i1=1,iFree%nD
      ! parms%D(iFree%D) = xFree(iFree%xD)
      write(Results1P_UNIT,1) HHData%ColumnLabels(iFree%xD(i1)), &
                              xFree(iFree%xD(i1)),  &
                              Grad(iFree%xD(i1)),   &
                              StandardErrors(iFree%xD(i1))
    end do
    do i1=1,iFree%nBC
      ! parms%CSig(iFree%BC) = xFree(iFree%xBC)
      write(cTemp,'(a4,i4)') 'C',i1
      write(Results1P_UNIT,1) cTemp,          &
                              xFree(iFree%xBC(i1)),  &
                              Grad(iFree%xBC(i1)),   &
                              StandardErrors(iFree%xBC(i1))
    end do
    do i1=1,iFree%nMUE
      ! parms%MUE(iFree%MUE) = xFree(iFree%xMUE)
      write(cTemp,'(a4,i4)') 'MUE_',i1
      write(Results1P_UNIT,1) cTemp,          &
                              xFree(iFree%xmUE(i1)),  &
                              Grad(iFree%xMUE(i1)),   &
                              StandardErrors(iFree%xMUE(i1))
    end do
    do i1=1,iFree%nINVCDiag
      ! parms%INVCDiag(iFree%INVCDiag) = xFree(iFree%xINVCDiag)
      write(cTemp,'(a7,i4)') 'INVCDiag',i1
      write(Results1P_UNIT,1) cTemp,          &
                              xFree(iFree%XINVCDiag(i1)),  &
                              Grad(iFree%XINVCDiag(i1)),   &
                              StandardErrors(iFree%XINVCDiag(i1))
    end do
    do i1=1,iFree%nINVCOffDiag
      ! parms%INVCOffDiag(iFree%INVCOffDiag) = xFree(iFree%XINVCOffDiag)
      write(cTemp,'(a9,i4)') 'INVCOffDiag',i1
      write(Results1P_UNIT,1) cTemp,          &
                              xFree(iFree%XINVCOffDiag(i1)),  &
                              Grad(iFree%XINVCOffDiag(i1)),   &
                              StandardErrors(iFree%XINVCOffDiag(i1))
    end do

    1 format( a20,3es25.16)
    2 format( a20,3a25   )
    ! write stats: (likelihood,nobs,maxgrad,MinEigHess,BIC,nNonZero)
    write(Results1P_UNIT,'(A12,g11.4)') 'Likelihood',LValue
    write(Results1P_UNIT,'(A12,I11)')   'nobs',stats%N
    write(Results1P_UNIT,'(A12,I11)')   'No. Parms',stats%nNonZero
    write(Results1P_UNIT,'(A12,g11.4)') 'BIC',stats%BIC
    write(Results1P_UNIT,'(A12,g11.4)') 'Max. Grad.',stats%MaxGrad
    write(Results1P_UNIT,'(A12,g11.4)') 'Min. Eig.',stats%MinEigHess
    if (Penalty%method>0) then
      write(Results1P_UNIT,'(A12,g11.4)') 'Penalty',Penalty%lambda
    end if
    close(UNIT = Results1P_UNIT)
  end if
  deallocate(InvHess,StandardErrors)
end subroutine SavePenaltyOutputs

! Extract short input filename and append it to OutDir
character(len=200) function CopyInputFilename(InputFile)
  implicit none
  character(len=*), intent(in) :: InputFile
  integer(i4b)                  :: i1,n1

  i1 = index(InputFile,"/",back=.TRUE.)
  n1 = len_trim(InputFile)
  CopyInputFilename = MakeFullFileName(InputFile(i1+1:n1))

end function CopyInputFilename

character(LEN=1048) function MakeFullFileName(ShortFileName)
  implicit none
  character(len=*), intent(in) :: ShortFileName

  MakeFullFileName = trim(OutDir) // '/' // trim(ShortFileName)
end function MakeFullFileName

subroutine ComputeStats(x,L,Grad,Hess,N,Stats)
  use GlobalModule, only : ResultStructure
  implicit none
  real(dp), intent(in)     :: x(:)
  real(dp), intent(in)     :: L
  real(dp), intent(in)     :: Grad(:)
  real(dp), intent(in)     :: Hess(:,:)
  integer(i4b), intent(in) :: N
  type(ResultStructure), intent(inout) :: stats

  ! variables needed to compute eigenvalues of Hess
  integer(i4b) :: nH,lwork,ifail
  character(1) :: JOBZ,UPLO
  real(dp), allocatable :: TempH(:,:),eig(:),work(:)
  real(dp)              :: dummy(1)

  stats%nNonZero   = count(abs(x)>=1e-6)
  stats%BIC        = 2.0*L + stats%nNonZero*dlog(dble(N))
  stats%MaxGrad    = maxval(abs(grad))
  stats%N          = N

  ! Compute eigenvalues of Hess
  nH =size(Hess,1)
  allocate(TempH(nH,nH),eig(nH))
  JOBZ = 'N'  ! 'N' = compute eigenvalues only and not vectors
  UPLO = 'L'  ! 'L' = ignore lower triangular part of matrix
  TempH = Hess 
  
  ! determine optimal workspsace size fo F08FAF
  lwork=-1
  call F08FAF(JOBZ,UPLO,nH,TempH,nH,eig,dummy,lwork,ifail)
  lwork = max(66*nH,nint(dummy(1)))
  allocate(work(lwork))

  ! compute eigenvalues
  call F08FAF(JOBZ,UPLO,nH,TempH,nH,eig,work,lwork,ifail)

  ! minimum eigenvalue
  stats%MinEigHess = minval(abs(eig))
  deallocate(TempH,eig,work)
end subroutine ComputeStats

subroutine SaveData
  use GlobalModule, only : HHData,parms
  implicit none
  integer(i4b) :: i1,i2

  open(UNIT = SaveData_UNIT_q,    &
       File = SaveDataFile_q, &
       Action = 'write')
  do i1=1,HHData%N
    write(SaveData_UNIT_q,161) HHData%q(:,i1)
    161 format(<parms%K>(g25.16,','))
  end do
  close(SaveData_UNIT_q)

  open(UNIT = SaveData_UNIT_p,    &
       File = SaveDataFile_p, &
       Action = 'write')
  do i1=1,HHData%N
    write(SaveData_UNIT_p,374) HHData%p(:,i1)
    374 format(<parms%J>(g25.16,','))
  end do
  close(SaveData_UNIT_p)

  open(UNIT = SaveData_UNIT_e,    &
       File = SaveDataFile_e, &
       Action = 'write')
  do i1=1,HHData%N
    write(SaveData_UNIT_e,382) HHData%e(:,i1)
    382 format(<parms%K>(g25.16,','))
  end do
  close(SaveData_UNIT_e)

  open(UNIT = SaveData_UNIT_eta,    &
       File = SaveDataFile_eta, &
       Action = 'write')
  do i1=1,HHData%N
    write(SaveData_UNIT_eta,396) HHData%eta(:,i1)
    396 format(<parms%dim_eta>(g25.16,','))
  end do
  close(SaveData_UNIT_e)

  open(UNIT = SaveData_UNIT_iNonZero,    &
       File = SaveDataFile_iNonZero, &
       Action = 'write')
  do i1=1,HHData%N
    write(SaveData_UNIT_iNonZero,392) HHData%iNonZero(:,i1)
    392 format(<parms%K>(g25.16,','))
  end do
  close(SaveData_UNIT_iNonZero)

  open(UNIT = SaveData_UNIT_iZero,    &
       File = SaveDataFile_iZero, &
       Action = 'write')
  do i1=1,HHData%N
    write(SaveData_UNIT_iZero,400) HHData%iZero(:,i1)
    400 format(<parms%J>(g25.16,','))
  end do
  close(SaveData_UNIT_iZero)
 
  open(UNIT = SaveData_UNIT_nNonZero,    &
       File = SaveDataFile_nNonZero, &
       Action = 'write')
  do i1=1,HHData%N
    write(SaveData_UNIT_nNonZero,409) HHData%nNonZero(i1)
    409 format(g25.16)
  end do
  close(SaveData_UNIT_nNonZero)

  open(UNIT = SaveData_UNIT_market,    &
       File = SaveDataFile_market, &
       Action = 'write')
  do i1=1,HHData%N
    write(SaveData_UNIT_market,'(g25.16)') HHData%market(i1)
  end do
  close(SaveData_UNIT_market)
end subroutine SaveData

! print file with (lambda,BIC,nNonZero)
subroutine SaveBICResults(model,stats)
  use GlobalModule, only : ResultStructure,Penalty
  implicit none
  integer(i4b),          intent(in) :: model
  type(ResultStructure), intent(in) :: stats(:)
  integer(i4b)                      :: i1
  integer(i4b)                      :: Current_UNIT
  character(len=200)                :: Current_FILE
  if (model==1) then
    Current_UNIT = BIC1_UNIT
    Current_FILE = BIC1_FILE
  end if

  open(UNIT = Current_UNIT, &
       FILE = Current_FILE, &
       action = 'write')
  write(Current_UNIT,439) 'Lambda','BIC','No. Parms'
  439 format(3A8)
  do i1=1,Penalty%nLambda
    write(Current_UNIT,442) Penalty%VectorLambda(i1),stats(i1)%BIC,stats(i1)%nNonZero
    442 format(2es12.4,i12)
  end do
  close(Current_UNIT)
end subroutine SaveBICResults

subroutine SaveMCOutputs(model,MCX,MCLambda,IMC)
  implicit none
  integer(i4b), intent(in) :: model
  real(dp),     intent(in) :: MCX(:,:),MCLambda(:)
  integer(i4b), intent(in) :: IMC
  integer(i4b)             :: i1,NX,NMC

  ! write values of MCX
  open(UNIT = MC_UNIT1, &
       FILE = MC_FILE1, &
       ACTION = 'write')

  NX = size(MCX,1)
  NMC = size(MCX,2)
  do i1=1,NX
    write(MC_UNIT1,471) MCX(i1,1:IMC),MCX(i1,NMC)
    471 format(<IMC+1>es25.16)
  end do

  close(MC_UNIT1)

  ! write values of MCLambda
  open(UNIT = MC_UNIT2, &
       FILE = MC_FILE2, &
       action = 'write')
  do i1=1,IMC
    write(MC_UNIT2,'(es25.16)') MCLambda(i1)
  end do
  close(MC_UNIT2)

  ! What outputs for MC exercise
  !    1) MSE of estimator =  mean( (beta_hat - beta).*(beta_hat-beta) )
  !    2) size of model selected
  !    3) probability selected model contains all large coefficients
  !    4) size of lambda
  ! 3) 
end subroutine SaveMCOutputs

subroutine WriteBayes(DINEST,ERREST,IVALID)
  implicit none
  real(dp), intent(in) :: DINEST(:),ERREST(:)
  integer(i4b), intent(in) :: IVALID(:)
  integer(i4b)             :: n,i1

  ! write bayes results
  open(UNIT = BayesResults_UNIT, &
       FILE = BayesResults_FILE, &
       ACTION = 'write')

  N = size(DINEST,1)
  write(BayesResults_UNIT,'(a4,2a25,a8)') 'i1','Integral','Error','status'
  do i1=1,N
    write(BayesResults_UNIT,514) i1,DINEST(i1),ERREST(i1),IVALID(i1)
    514 format(i4,2es25.16,i4)
  end do

  close(BayesResults_UNIT)
  
end subroutine WriteBayes

subroutine WriteElasticities(elas)
  implicit none
  real(dp), intent(in) :: elas(:,:)

  integer(i4b)         :: J,j1

  open(UNIT = Elas_UNIT, &
       FILE = ElasFILE, &
       ACTION = 'WRITE')
  J = size(elas,1)
  do j1=1,J
    write(Elas_UNIT,534) elas(:,j1)
  end do
534 format(<J>(g25.16,','))

  close(Elas_UNIT)
end subroutine WriteElasticities

subroutine WriteDemandResults1(qdata,qhat,qaverage)
  implicit none
  real(dp), intent(in) :: qdata(:),qhat(:),qaverage(:)

!  character(len=200)       :: ShortFileName
  integer(i4b)             :: J,i1

!  write(ShortFileName,*) 'demand_data.csv'
  DemandFile = MakeFullFileName(trim('demand_data.csv'))

  open(unit = Demand_UNIT, &
       File = DemandFile,  &
       Action = 'WRITE')
  J = size(qdata)
  do i1=1,J
    write(Demand_UNIT,562) qdata(i1),qhat(i1),qaverage(i1)
  end do
562 format(3(g25.16,','))
  close(Demand_Unit)

end subroutine WriteDemandResults1


subroutine WriteDemandResults2(p,NewQ,j1)
  implicit none
  real(dp), intent(in) :: p(:),NewQ(:,:)
  integer(i4b), intent(in) :: j1

  character(len=200)       :: ShortFileName
  integer(i4b)             :: i1,np

  write(ShortFileName,535) 'demand',j1,'.csv'
535 format(a6,i2.2,a4)
  DemandFile = MakeFullFileName(trim(ShortFileName))

  open(unit = Demand_UNIT, &
       File = DemandFile,  &
       Action = 'WRITE')
  np = size(p)
  do i1=1,np
    write(Demand_UNIT,536) p(i1),NewQ(i1,:)
  end do
536 format(7(g25.16,','))
  close(Demand_Unit)

end subroutine WriteDemandResults2

! Write tax results: aggregate results:  baseline (q0,p0) and alternative (qtax,ptax) (J x ...)
!  (q0,p0)     baseline quantities and prices
!  (qtax,ptax) counterfactual quantitites and prices
subroutine WriteTaxResults1(q0,p0,qtax,ptax,filename)
  implicit none
  real(dp),         intent(in) :: q0(:),p0(:),qtax(:,:),ptax(:,:)
  character(len=*), intent(in) :: filename
  character(len=200)           :: resultsfile
  integer(i4b)                 :: ntax,i1,J
  character(len=4), allocatable :: qstring(:),pstring(:)

  resultsfile = MakeFullFileName(trim(filename))

  open(unit = taxresults_UNIT, &
       File = resultsfile,    &
       Action = 'WRITE')

  J    = size(q0,1)
  ntax = size(qtax,2)

  ! Create variable names
  allocate(qstring(ntax),pstring(ntax))
  do i1=1,ntax
    write(qstring(i1),'(a1,i2.2)') "q",i1
    write(pstring(i1),'(a1,i2.2)') "p",i1
  end do

  ! write variable names 
  write(taxresults_unit,780) "q0",qstring,"p0",pstring
780 format(<2*ntax+2>(a25,:,","))

  ! write (quantity,price)
  do i1=1,J
    write(taxresults_unit,781) q0(i1),qtax(i1,:),p0(i1),ptax(i1,:)
  end do
781 format(<2*ntax+2>(g25.16,:,","))

  close(taxresults_unit)
  deallocate(qstring,pstring) 
end subroutine WriteTaxResults1

! Write tax results: household (expenditure,utility)
!    (e0,u0)     = baseline results
!    (etax,utax) = counterfactual results
subroutine WriteTaxResults2(e0,u0,etax,utax,filename)
  implicit none
  real(dp),         intent(in) :: e0(:),u0(:),etax(:,:),utax(:,:)
  character(len=*), intent(in) :: filename
  character(len=200)           :: resultsfile
  integer(i4b)                 :: i1,n,ntax
  character(len=4), allocatable :: estring(:),ustring(:)

  resultsfile = MakeFullFileName(trim(filename))

  n    = size(e0,1)
  ntax = size(etax,2)
  ! create variable labels
  allocate(estring(ntax),ustring(ntax))
  do i1=1,ntax
    write(estring(i1),'(a1,i2.2)') "e",i1
    write(ustring(i1),'(a1,i2.2)') "u",i1
  end do

  open(unit = taxresults_UNIT, &
       File = resultsfile,    &
       Action = 'WRITE')

  ! Write variable names
  write(taxresults_unit,826) "e0",estring,"u0",ustring
826 format(<2*ntax+2>(a25,:,","))

  ! write (expenditure,utility)
  do i1=1,n
    write(taxresults_unit,830) e0(i1),etax(i1,:),u0(i1),utax(i1,:)
  end do
830 format(<2*ntax+2>(g25.16,:,","))
  close(taxresults_unit)
  deallocate(ustring,estring)
end subroutine WriteTaxResults2

end module OutputModule

