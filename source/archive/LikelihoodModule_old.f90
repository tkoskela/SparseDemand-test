module LikelihoodModule

! 1)  subroutine DefineFreeParameters(K,RCDIM,model)
! 2)  subroutine ComputeInitialGuess(model,b,CDiag,COffDiag,x)
! 3)  subroutine ComputeIntegral(xData,b,C,J,RCDim,mode,Prob,GradB,GradC)
! 4)  subroutine ComputeLogitProbability(xData,b,J,mode,Prob,GradB)
! 5)  subroutine MapToStructural(x,RCDim,b,CDiag,COffDiag,C)
! 6)  subroutine SetBounds(model,lower,upper)
! 7)  subroutine MaximizeLikelihood(model,x,LValue,Grad,Hess,ierr)
! 8)  subroutine MaximizePenalizedLikelihood(model,x,LValue,Grad,Hess,ierr)
! 9)  subroutine LikeFunc(mode,nFree,xFree,L,GradL,nstate,iuser,ruser)
! 10) subroutine MixedLogit(mode,nx,x,L,GradL,nstate,iuser,ruser)
! 11) subroutine PenalizedLikelihood(mode,nxP,xP,LP,GradLP,nstate,iuser,ruser)
! 12) subroutine PenaltyFunction(xPlus,xMinus,P,GradP)
! 13) subroutine Logit(mode,nx,x,L,GradL,nstate,iuser,ruser)

contains

subroutine DefineFreeParameters(model)
  use nrtype
  use GlobalModule, only : K,RCDIM,iFree,Penalty
  implicit none

  integer(i4b), intent(in)  :: model

  ! Define number of free parameters
  iFree%nb        = K
  allocate(iFree%b(iFree%nb),iFree%bX(iFree%nb))
  !  b(iFree%b)=x(iFree%bX)
  iFree%b    = (/1:K/)
  iFree%bX   = iFree%b
  iFree%nAll = iFree%nb

  if (model==2) then
    iFree%nCDiag = RCDim
    allocate(iFree%CDiag(iFree%nCDiag))
    allocate(iFree%CDiagX(iFree%nCDiag))
    iFree%CDiag  = (/1:RCDim/)
    iFree%CDiagX = iFree%nb + (/1:iFree%nCDiag/)

    if (Penalty%method==0) then
      iFree%nCOffDiag = RCDim*(RCDim-1)/2
      allocate(iFree%COffDiag(iFree%nCOffDiag))
      allocate(iFree%COffDiagX(iFree%nCOffDiag))
      iFree%COffDiag  = (/1:iFree%nCOffDiag/)
      iFree%COffDiagX = iFree%nb + iFree%nCDiag + (/1:iFree%nCOffDiag/)
    else if (Penalty%method==1) then
      iFree%nCOffDiag = 0
    end if
    iFree%nAll = iFree%nAll + iFree%nCDiag + iFree%nCOffDiag
  end if

end subroutine DefineFreeParameters

subroutine ComputeInitialGuess(model,x)
  use nrtype
  use GlobalModule, only : iFree,b,CDiag,COffDiag,Penalty
  implicit none

  integer(i4b), intent(in)  :: model
  real(dp),     intent(out) :: x(:)

  ! initial guess
  x = 0.0d0
  x(iFree%bX)        = b(iFree%b)
  !x(iFree%bX) = 1.0d0

  if (model==2) then
    x(iFree%CDiagX)    = CDiag(iFree%CDiag)
    if (iFree%nCOffDiag>0) then
      x(iFree%COffDiagX) = COffDiag(iFree%COffDiag)
      if (Penalty%method==1) then
        COffDiag = 0.0d0
      end if
    end if
  end if
end subroutine ComputeInitialGuess

subroutine ComputeIntegral(xData,b,C,J,RCDim,mode,Prob,GradB,GradC)
  use nrtype
  use GlobalModule, only : iXData,IntRule
  implicit none

  ! xData  (K x J x N)
  real(dp),               intent(in)  :: xData(:,:)
  real(dp),               intent(in)  :: b(:)
  real(dp),               intent(in)  :: C(:,:)
  integer(i4b),           intent(in)  :: J
  integer(i4b),           intent(in)  :: RCDim
  integer(i4b),           intent(in)  :: mode
  real(dp),               intent(out) :: Prob
  real(dp),               intent(out) :: GradB(:)
  real(dp),               intent(out) :: GradC(:)

  ! Prob  = probability of choosing option J
  ! GradB = gradient of P w.r.t. B
  ! GradC = gradient of P w.r.t. C
  !
  ! xData = (K x J)  data matrix
  ! 
  ! Revision history
  ! 06jul2012 LN   translate from ComputeIntegral.m
  !
  ! Remaining jobs
  ! 2) code everything for flag==1
  ! 3) code everything for flag==3
  ! 4) code everything for flag==4
  ! 5) code everything for flag==5
  !
  ! IntRule.flag = 0  Gauss-Hermite
  !              = 1  Sparse Grid Hermite
  !              = 2  pseudo-Montecarlo
  !              = 3  quasi-montecarlo
  !              = 4  Taylor expansion around 0
  !              = 5  Laplace
  integer(i4b) :: nb
  real(dp), allocatable :: v(:),v1(:),ev1(:)
  real(dp), allocatable :: e(:,:)
  integer(i4b) :: i1,i2,i3,ib,iC,row,col
  real(dp)     :: vMax
  real(dp)     :: TempProb(J)
  real(dp)     :: m1_b(size(b,1))
  real(dp)     :: m1_e(RCDim)
  real(dp)     :: DeltaX(RCDim,J)
  real(dp)     :: m2(RCDim)
  real(dp)     :: t0,p2
  real(dp)     :: DeltaX_b(size(b,1),J)
  real(dp)     :: m2_b(RCDim,size(b,1))
  real(dp)     :: m3_b(RCDim,size(b,1))
  real(dp)     :: t1(RCDim,J)
  real(dp)     :: p2_b

  allocate(v(J),v1(J),ev1(J))
  allocate(e(IntRule%nAll,RCDIM))

  nb = size(b,1)
  GradB = 0.0d0
  GradC = 0.0d0

  if (IntRule%flag<4) then
    ! size(v) = (J x 1)
    v = matmul(b,xData(iXData%b,:))
  
    e = matmul(IntRule%nodes,transpose(C))
    Prob = 0.0d0

    do i1=1,IntRule%nAll
      v1 = v +  matmul(transpose(xData(iXData%e,:)),e(i1,:))
      vMax = maxval(v1)
      ev1 = dexp(v1 - vMax)
      if (mode>0) then
        TempProb = ev1/sum(ev1)
        Prob = Prob + IntRule%weights(i1)*TempProb(J)

        m1_b  = matmul(xData(iXData.b,:),TempProb)
        GradB = GradB + IntRule%weights(i1)*            &
                TempProb(J)*(xData(iXData.b,J)-m1_b)
       
        ! m1_e = (RCDim x 1)
        m1_e  = matmul(xData(iXData.e,:),TempProb)
      
        ! GradC(iC) = Gradient P w.r.t C(row,col) 
        do row=1,RCDim
        do col=1,row
          iC = row*(row-1)/2+col
          GradC(iC) = GradC(iC) + IntRule%weights(i1)* &
                TempProb(J)*IntRule%nodes(i1,col)*(xData(iXData.e(row),J)-m1_e(row))
        end do
        end do
      
      elseif (mode==0) then
        Prob = Prob + IntRule%weights(i1)*ev1(J)/sum(ev1) 
      end if
    end do  ! do i1=1,nAll

  elseif (IntRule%flag==4) then
    ! size(v) = (J x 1)
    v = matmul(b,xData(iXData%b,:))
    vMax = maxval(v)
    ev1  = dexp(v-vMax)

    ! ProbVec = (J x 1)  base prob of each event
    TempProb = ev1/sum(ev1)

    ! m1 = (RCDim x 1)    =  weighted average of x
    m1_e = matmul(xData(iXData%e,:),TempProb)  

    ! m2 = (RCDim x 1) =  weighted second moment of x
    !                     only diagonal elements
    ! DeltaX = (RCDim x J)

    DeltaX = xData(iXData%e,:) - spread(m1_e,2,J)
    m2 = matmul(DeltaX*DeltaX,TempProb)
  
    t0 = 1.0d0
    do i1=1,RCDim
      ! D2_P w.r.t. C(i1,i2)
      p2 = DeltaX(i1,J)*DeltaX(i1,J)-m2(i1)
      do i2=1,i1
        ! expansion w.r.t. C(i1,i2)
        t0 = t0 + 0.5d0*p2*C(i1,i2)**2
        iC = i1*(i1-1)/2+i2
        GradC(iC) = GradC(iC) + p2*C(i1,i2)
      end do
    end do
  
    Prob = TempProb(J)*t0
  
    GradC = GradC*TempProb(J)
 

    if (mode>0) then
      ! Gradient w.r.t. b
      ! DeltaX_b  =  (nb x J)
      ! m1_b      =  (nb x 1)
      m1_b     = matmul(xData(iXData%b,:),TempProb)
      DeltaX_b = xData(iXData%b,:) - spread(m1_b,2,J)
    
      m2_b = 0.0d0
      m3_b = 0.0d0
    
      GradB = DeltaX_b(:,J)
      do ib=1,nb
        t1 = DeltaX * spread(DeltaX_b(ib,:),1,RCDim)
        m2_b(:,ib) = matmul(t1,TempProb)
        m3_b(:,ib) = matmul(DeltaX*t1,TempProb)
      
        do i2=1,RCDim
          p2_b = DeltaX(i2,J)*DeltaX(i2,J)*DeltaX_b(ib,J) &
                 - 2.0d0*m2_b(i2,ib)*DeltaX(i2,J)         &
                 - m2(i2)*DeltaX_b(ib,J)                  &
                 - m3_b(i2,ib)
          do i3=1,i2
            GradB(ib) = GradB(ib) + 0.5d0*p2_b*C(i2,i3)**2    
          end do
        end do     
      end do
      GradB = TempProb(J)*GradB
    end if  
  end if
  deallocate(e,v,v1,ev1)
end subroutine ComputeIntegral

! Compute logit probability
subroutine ComputeLogitProbability(xData,b,J,mode,Prob,GradB)
  use nrtype
  use GlobalModule, only : iXData
  implicit none

  ! xData  (K x J x N)
  real(dp),               intent(in)  :: xData(:,:)
  real(dp),               intent(in)  :: b(:)
  integer(i4b),           intent(in)  :: J
  integer(i4b),           intent(in)  :: mode
  real(dp),               intent(out) :: Prob
  real(dp),               intent(out) :: GradB(:)

  ! Prob  = probability of choosing option J
  ! GradB = gradient of P w.r.t. B
  !
  ! xData = (K x J)  data matrix
  ! 
  ! Revision history
  ! 19sep2012 LN adapt from ComputeIntegral
  ! 06jul2012 LN   translate from ComputeIntegral.m
  !
  integer(i4b) :: nb
  real(dp)     :: v(J),v1(J),ev1(J)
  integer(i4b) :: i1,i2,i3,ib,iC,row,col
  real(dp)     :: vMax
  real(dp)     :: TempProb(J)
  real(dp)     :: m1_b(size(b,1))
  real(dp)     :: t0,p2
  real(dp)     :: DeltaX_b(size(b,1))
  real(dp)     :: p2_b

  nb = size(b,1)
  GradB = 0.0d0

  ! size(v) = (J x 1)
  v = matmul(b,xData(iXData%b,:))
  vMax = maxval(v)
  ev1  = dexp(v-vMax)

  ! ProbVec = (J x 1)  base prob of each event
  TempProb = ev1/sum(ev1)
  Prob = TempProb(J)

  if (mode>0) then
    ! Gradient w.r.t. b
    ! DeltaX_b  =  (nb x 1)
    ! m1_b      =  (nb x 1)
    m1_b     = matmul(xData(iXData%b,:),TempProb)
    DeltaX_b = xData(iXData%b,J) - m1_b
  
    GradB = TempProb(J)*GradB
  end if  
end subroutine ComputeLogitProbability

subroutine MapToStructural(x,RCDim,b,CDiag,COffDiag,C)
  use nrtype
  use GlobalModule, only : iFree
  implicit none

  real(dp),             intent(in)    :: x(:)
  integer(i4b),         intent(in)    :: RCDim
  real(dp),             intent(inout) :: b(:)
  real(dp),             intent(inout) :: CDiag(:)
  real(dp),             intent(inout) :: COffDiag(:)
  real(dp),             intent(out)   :: C(:,:)

  integer(i4b)  :: i1

  ! Map free parameters to structural parameters
  !
  !   [b,C]=MapToStructural(x,iFree,b,C)
  !
  ! Revision History
  ! 06jul2012 LN  translate from MapToStructural.m

  ! Vector of utility coefficients
  !    x    free parameters
  !    b     structural utility parameters
  !    iFree    indexes used to map xFree2 to beta2 and C

  if (size(iFree%b,1)>0) then
    b(iFree%b)=x(iFree%bX)		
  end if

  ! Diagonal elements of C = chol(sig) where sig is covariance matrix of 
  ! random coefficients
  if (size(iFree%CDiag)>0) then
    CDiag(iFree%CDiag)=x(iFree%CDiagX)
  end if

  ! OffDiagonal elements of C
  if (size(iFree%COffDiag,1)>0) then
    COffDiag(iFree%COffDiag) = x(iFree%COffDiagX)   
  end if

  C = 0.0d0
  do i1=1,RCDim
    C(i1,i1) = CDiag(i1)
    if (i1>1) then
      C(i1,1:i1-1) = COffDiag((i1-1)*(i1-2)/2+(/1:i1-1/))
    end if
  end do
end subroutine MapToStructural

subroutine SetBounds(model,lower,upper)
  use nrtype
  use GlobalModule, only : iFree
  implicit none

  integer(i4b),         intent(in)  :: model
  real(dp),             intent(out) :: lower(:)
  real(dp),             intent(out) :: upper(:)

  ! Set lower and upper bounds for parameters
  lower = 0.0d0
  upper = 0.0d0

  lower(iFree%bX) = -40.0d0
  upper(iFree%bX) = 40.0d0
  if (model==2) then
    lower(iFree%CDiagX)    = 0.0d0
    upper(iFree%CDiagX)    = 2.0d0
    lower(iFree%COffDiagX) = -2.0d0
    upper(iFree%COffDiagX) = 2.0d0
  end if
end subroutine SetBounds

subroutine MaximizeLikelihood(model,x,LValue,Grad,Hess,ierr)
  use nrtype
  use GlobalModule, only : MaxOptions,Penalty
  implicit none
  integer(i4b), intent(in)    :: model
  real(dp),     intent(inout) :: x(:)
  real(dp),     intent(out)   :: LValue
  real(dp),     intent(out)   :: Grad(:)
  real(dp),     intent(out)   :: Hess(:,:)
  integer(i4b), intent(out)   :: ierr

  if (MaxOptions%Algorithm(model)==1) then
    ! E04WDF:  Dense optimization: constrained non-linear optimization for dense problem
    if (Penalty%method==0) then
      call MaximizeLikelihood1(model,x,LValue,Grad,Hess,ierr)
    else if (Penalty%method==1) then
      call MaximizePenalizedLikelihood(model,x,LValue,Grad,Hess,ierr)
    end if
  elseif (MaxOptions%Algorithm(model)==2) then
    ! E04VHF:  Sparse optimization: constrained non-linear optimization for sparse problem
!    call MaximizeLikelihood2(model,x,LValue,Grad,Hess,ierr)
  end if
end subroutine MaximizeLikelihood

! Maximise using Constrained maximization: E04WDF
subroutine MaximizeLikelihood1(model,x,LValue,Grad,Hess,ierr)
  use nrtype
  use GlobalModule, only : ControlOptions, &
                           MaxOptions
  implicit none
  integer(i4b), intent(in) :: model
  real(dp), intent(inout)  :: x(:)
  real(dp), intent(out)    :: LValue
  real(dp), intent(out)    :: Grad(:)
  real(dp), intent(out)    :: Hess(:,:)
  integer(i4b), intent(out) :: ierr

  ! Default values of parameters required for NAG routine e04WDF
  integer(i4b), parameter   :: leniw=600,lenrw=600
  integer(i4b), parameter   :: eflag=1,eunit=-1

  ! used to transmit completion message to slave in MPI version
  integer(i4b)              :: mode

  ! Dimensions of optimization problem
  integer(i4b)              :: nx,nclin,ncnln,ldcj,ldh,lda,nctotal

  ! linear constraints on optimization problem
  real(dp), allocatable     :: A(:,:)

  ! LOWER AND UPPER BOUNDS for optimisation
  real(dp), allocatable     :: BL(:),BU(:)            

  ! MULTIPLIERS ON CONSTRAINTS
  REAL(DP), allocatable     :: CLAMBDA(:)                          

  ! PARAMETERS FOR E04WCF and E04WDF (constrained optimisation subroutines from NAG)
  INTEGER(I4B)              :: IFAIL,iter
  integer(i4b), allocatable :: ISTATE(:)
  real(dp), allocatable     :: ccon(:)
  real(dp), allocatable     :: cjac(:,:)
  
  ! Workspace for E04WCF and E04WDF (constrained optimisation subroutines from NAG)
  real(dp)                  :: RW(lenrw)
  integer(i4b)              :: IW(leniw)
  
  ! E04WDF outputs
  real(DP)                  :: OBJF

  ! User INPUTS TO OBJECTIVE FUNCTION
  REAL(DP)                  :: RUSER(1)
  INTEGER(I4B)              :: iuser(1)

  ! Control parameter for saving intermediate results
  !   0 DO NOT SAVE intermediate results
  !   1 SAVE intermediate results
  integer(i4b), parameter   :: NAG_SAVE=0    
    
  external E04WCF		! initialisation routine for E04WDF
  external E04WDP		! dummy subroutine to evaluate nonlinear constraints
  external E04WDF  	! constrained optimisation
  external E04WFF		! set options for E04WDF
  external X04AAF		! suppress and/or redirect NAG error messages

  ! Dimensions of maximization problem
  nx=size(x,1)
  nclin=0
  ncnln=0
  lda=max(1,NCLIN)
  ldcj=max(1,NCNLN)
  ldh=size(HESS,1)
  nctotal=nx+nclin+ncnln

  ! Allocate memory for linear constraints, bound constraints, and some aux. outputs
  allocate(A(lda,1))
  allocate(BL(nctotal),BU(nctotal))
  allocate(CLAMBDA(nctotal),ISTATE(nctotal))
  allocate(CJAC(LDCJ,1))
  allocate(CCON(max(1,NCNLN)))

  ! Initialize E04WDF by calling E04WCF
  ifail=-1
  call e04wcf(IW,LENIW,RW,LENRW,ifail)

  ! Suppress NAG error messages
  call x04aaf(eflag,eunit)

  ! Constraints
  !   NONE

  ! constraints:   BL <= A*x <= BU
  !  A is not referenced when nclin==0
		 
  ! lower and upper bounds on xFree 
  call SetBounds(model,BL,BU)
   
  ! Option to verify gradient calculations
  call E04WFF('Verify Level = -1',iw,rw,ifail)     ! -1   do not check derivatives
                                                  ! 0    default  cheap derivative check
                                                  ! 3 comprehensive check
  !open(1,File=*,'WRITE')
  !call E04WFF('Major Optimality Tolerance = 2.0D-8',iw,rw,ifail)
  !call E04WFF('Minor Feasibility Tolerance = 2.0D-8',iw,rw,ifail) 
  !MaxIter = 'Major Iterations Limit = 5'

  !call E04Wff('Major Iterations Limit = 5',iw,rw,ifail)            
  call E04WFF(trim(MaxOptions%MaxIter(model)),iw,rw,ifail)            
  call E04WFF('Major Print Level = 00001',iw,rw,ifail)  ! default = 00001
  call E04Wff('Print File = 6',iw,rw,ifail)            
  !call E04WFF('Minor Print Level = 0',iw,rw,ifail)      ! default = 1
  !call E04WFF('Print Frequency = 1',iw,rw,ifail)         ! default = 100
 ! call E04WFF('Summary File = 6',iw,rw,ifail) 
 ! call E04WFF('Summary Frequency = 1',iw,rw,ifail)
  ! call E04WFF('Derivative Linesearch')    ! default
  !call E04WFF('Nonderivative Linesearch',iw,rw,ifail)   ! possibly faster, but less accurate
  call E04WFF('Hessian Full Memory',iw,rw,ifail)      ! more accurate, more meory
  ! call E04WFF('Hessian Limited Memory',iw,rw,ifail)   ! default when n>75
  !call E04WFF('Check Frequency = 600',iw,rw,ifail)      ! default = 60
  !call E04WFF('Linesearch Tolerance = 0.99',iw,rw,ifail) ! default=0.9
                                                         ! 0.99 faster, less accurate
                                                         ! 0.1 or 0.01 slower, more accurate
  if (NAG_SAVE==1) then  
    if (model==1) then 
      open(101,FILE = 'output/NewBasis1.txt',ACTION='WRITE')
      open(103,FILE = 'output/BackupBasis1.txt',ACTION='WRITE')
      call E04WFF('New Basis File = 101',iw,rw,ifail)   ! unit number to save intermediate results
      call E04WFF('Backup Basis File = 103',iw,rw,ifail)   ! unit number to save intermediate results
    elseif (model>1) then
      open(102,FILE = 'output/NewBasis2.txt',ACTION='WRITE')
      open(104,FILE = 'output/BackupBasis2.txt',ACTION='WRITE')
      call E04WFF('New Basis File = 102',iw,rw,ifail)   ! unit number to save intermediate results
      call E04WFF('Backup Basis File = 104',iw,rw,ifail)   ! unit number to save intermediate results
    end if
    call E04WFF('Save Frequency = 100',iw,rw,ifail)
  end if

  ifail=-1
  iuser(1) = model
  ruser    = 0.0d0
 
  LValue = 0.0D0
  Grad   = 0.0d0
  Hess   = 0.0d0
  do iter=1,nx
    Hess(iter,iter)=1.0d0
  end do
 
 ! Maximise likelihood 
  if (ControlOptions%TestLikeFlag==1) then
    call TestGrad(LikeFunc,nx,x,1,iuser,ruser)
  else
    print *,'Begin maximization.'
    call E04WDF(nx,NCLIN,NCNLN,LDA,LDCJ,LDH,A,BL,BU,              &
                'E04WDP',LikeFunc,iter,ISTATE,CCON,CJAC,CLAMBDA,  &
                LValue,GRAD,HESS,x,IW,LENIW,RW,LENRW,iuser,RUSER, &
                IFAIL)
    call ComputeHess(x,nx,model,LValue,GRAD,Hess,iuser,ruser)
  end if

  deallocate(A,BL,BU,CLAMBDA,ISTATE,CJAC)
  deallocate(CCON)

  if (NAG_SAVE==1) then
    if (model==1) then
      close(101)
      close(103)
    elseif (model>=2) then
      close(102)
      close(104)
    end if
  end if

  ierr = 0   ! no error in subroutine
end subroutine MaximizeLikelihood1

! Maximise using Constrained maximization: E04WDF
subroutine MaximizePenalizedLikelihood(model,x,LValue,Grad,Hess,ierr)
  use nrtype
  use GlobalModule, only : Penalty,ControlOptions,MaxOptions
  implicit none
  integer(i4b), intent(in) :: model
  real(dp), intent(inout)  :: x(:)
  real(dp), intent(out)    :: LValue
  real(dp), intent(out)    :: Grad(:)
  real(dp), intent(out)    :: Hess(:,:)
  integer(i4b), intent(out) :: ierr

  ! xP  = [x;xPlus;xMinus]
  real(dp), allocatable :: xP(:)
  real(dp), allocatable :: GradP(:)
  real(dp), allocatable :: HessP(:,:)

  ! Default values of parameters required for NAG routine e04WDF
  integer(i4b), parameter   :: leniw=600,lenrw=600
  integer(i4b), parameter   :: eflag=1,eunit=-1

  ! used to transmit completion message to slave in MPI version
  integer(i4b)              :: mode

  integer(i4b)              :: i1
  ! Dimensions of optimization problem
  integer(i4b)              :: nxP,nclin,ncnln,ldcj,ldh,lda,nctotal

  ! linear constraints on optimization problem
  real(dp), allocatable     :: A(:,:)

  ! LOWER AND UPPER BOUNDS for optimisation
  real(dp), allocatable     :: BL(:),BU(:)            

  ! MULTIPLIERS ON CONSTRAINTS
  REAL(DP), allocatable     :: CLAMBDA(:)                          

  ! PARAMETERS FOR E04WCF and E04WDF (constrained optimisation subroutines from NAG)
  INTEGER(I4B)              :: IFAIL,iter
  integer(i4b), allocatable :: ISTATE(:)
  real(dp), allocatable     :: ccon(:)
  real(dp), allocatable     :: cjac(:,:)
  
  ! Workspace for E04WCF and E04WDF (constrained optimisation subroutines from NAG)
  real(dp)                  :: RW(lenrw)
  integer(i4b)              :: IW(leniw)
  
  ! E04WDF outputs
  real(DP)                  :: OBJF

  ! User INPUTS TO OBJECTIVE FUNCTION
  REAL(DP)                  :: RUSER(1)
  INTEGER(I4B)              :: iuser(1)

  ! Control parameter for saving intermediate results
  !   0 DO NOT SAVE intermediate results
  !   1 SAVE intermediate results
  integer(i4b), parameter   :: NAG_SAVE=0    
  
  
  external E04WCF		! initialisation routine for E04WDF
  external E04WDP		! dummy subroutine to evaluate nonlinear constraints
  external E04WDF        	! constrained optimisation
  external E04WFF		! set options for E04WDF
  external X04AAF		! suppress and/or redirect NAG error messages

  ! Dimensions of maximization problem
  nxP   = Penalty%nxP
  nclin = Penalty%nx2
  ncnln = 0
  lda   = max(1,NCLIN)
  ldcj  = max(1,NCNLN)
  ldh   = nxP
  nctotal=nxP+nclin+ncnln

  ! Allocate memory for linear constraints, bound constraints, and some aux. outputs
  allocate(A(lda,nxP))
  allocate(BL(nctotal),BU(nctotal))
  allocate(CLAMBDA(nctotal),ISTATE(nctotal))
  allocate(CJAC(LDCJ,1))
  allocate(CCON(max(1,NCNLN)))
  allocate(xP(nxP))
  allocate(GradP(nxP))
  allocate(HessP(nxP,nxP))

  ! Initialize E04WDF by calling E04WCF
  ifail=-1
  call e04wcf(IW,LENIW,RW,LENRW,ifail)

  ! Suppress NAG error messages
  call x04aaf(eflag,eunit)

  ! constraints:   BL <= A*xP <= BU
  !  A is not referenced when nclin==0
  !    xP = [x1;x2;xPlus;xMinus]
  !    0 <= x2 - xPlus+xMinus <=0
  ! size(A) = nx2 x nxP
  A = 0.0d0
  do i1=1,nclin
    A(i1,Penalty%xP_index2(i1)) = 1.0d0
    A(i1,Penalty%xP_index_xPlus(i1)) = -1.0d0
    A(i1,Penalty%xP_index_xMinus(i1)) = 1.0d0
  end do  
		 
  ! lower and upper bounds on xP
  ! x1L  <= x1              <= x1H
  ! x2L  <= x2              <= x2H
  ! 0    <= xPlus           <= x2H
  ! 0    <= xMinus          <= -x2L
  ! 0    <= x2-xPlus+xMinus <= 0
  BL = 0.0d0
  BU = 0.0d0

  ! set bounds on x = (x1,x2)
  call SetBounds(model,BL(1:Penalty%nx),BU(1:Penalty%nx))

  ! set bounds on (xPlus,xMinus)
  BU(Penalty%xP_index_xPlus)  = BU(Penalty%xP_index2)
  BU(Penalty%xP_index_xMinus) = -BL(Penalty%xP_index2)
 
  ! Option to verify gradient calculations
  call E04WFF('Verify Level = -1',iw,rw,ifail)     ! -1   do not check derivatives
                                                  ! 0    default  cheap derivative check
                                                  ! 3 comprehensive check
  !open(1,File=*,'WRITE')
  !call E04WFF('Major Optimality Tolerance = 2.0D-8',iw,rw,ifail)
  !call E04WFF('Minor Feasibility Tolerance = 2.0D-8',iw,rw,ifail) 
  !MaxIter = 'Major Iterations Limit = 5'

  !call E04Wff('Major Iterations Limit = 5',iw,rw,ifail)            
  call E04WFF(trim(MaxOptions%MaxIter(model)),iw,rw,ifail)            
  call E04WFF('Major Print Level = 00001',iw,rw,ifail)  ! default = 00001
  call E04Wff('Print File = 6',iw,rw,ifail)            
  !call E04WFF('Minor Print Level = 0',iw,rw,ifail)      ! default = 1
  !call E04WFF('Print Frequency = 1',iw,rw,ifail)         ! default = 100
 ! call E04WFF('Summary File = 6',iw,rw,ifail) 
 ! call E04WFF('Summary Frequency = 1',iw,rw,ifail)
  ! call E04WFF('Derivative Linesearch')    ! default
  !call E04WFF('Nonderivative Linesearch',iw,rw,ifail)   ! possibly faster, but less accurate
  call E04WFF('Hessian Full Memory',iw,rw,ifail)      ! more accurate, more meory
  ! call E04WFF('Hessian Limited Memory',iw,rw,ifail)   ! default when n>75
  !call E04WFF('Check Frequency = 600',iw,rw,ifail)      ! default = 60
  !call E04WFF('Linesearch Tolerance = 0.99',iw,rw,ifail) ! default=0.9
                                                         ! 0.99 faster, less accurate
                                                         ! 0.1 or 0.01 slower, more accurate
  if (NAG_SAVE==1) then  
    if (model==1) then 
      open(101,FILE = 'output/NewBasis1.txt',ACTION='WRITE')
      open(103,FILE = 'output/BackupBasis1.txt',ACTION='WRITE')
      call E04WFF('New Basis File = 101',iw,rw,ifail)   ! unit number to save intermediate results
      call E04WFF('Backup Basis File = 103',iw,rw,ifail)   ! unit number to save intermediate results
    elseif (model>1) then
      open(102,FILE = 'output/NewBasis2.txt',ACTION='WRITE')
      open(104,FILE = 'output/BackupBasis2.txt',ACTION='WRITE')
      call E04WFF('New Basis File = 102',iw,rw,ifail)   ! unit number to save intermediate results
      call E04WFF('Backup Basis File = 104',iw,rw,ifail)   ! unit number to save intermediate results
    end if
    call E04WFF('Save Frequency = 100',iw,rw,ifail)
  end if

  ifail=-1
  iuser(1) = model
  ruser    = 0.0d0
 
  LValue = 0.0D0
  Grad   = 0.0d0
  Hess   = 0.0d0
  GradP  = 0.0d0
  HessP  = 0.0d0
  do iter=1,nxP
    HessP(iter,iter)=1.0d0
  end do
  Hess = HessP(1:Penalty%nx,1:Penalty%nx)

  ! Initial guess for xP
  xP=0.0d0
  if (Penalty%nx1>0) then
    xP(Penalty%xP_index1) = x(Penalty%x_index1)
  end if
  xP(Penalty%xP_index2) = x(Penalty%x_index2)
  do i1=1,Penalty%nx2
    xP(Penalty%xP_index_xPlus(i1)) = max(x(Penalty%x_index2(i1)),0.0d0)
    xP(Penalty%xP_index_xMinus(i1)) = max(-x(Penalty%x_index2(i1)),0.0d0)
  end do

 ! Maximise likelihood 
  if (ControlOptions%TestLikeFlag==1) then
    call TestGrad(LikeFunc,nxP,xP,1,iuser,ruser)
  else
    print *,'Begin maximization.'
    call E04WDF(nxP,NCLIN,NCNLN,LDA,LDCJ,LDH,A,BL,BU,                &
                'E04WDP',LikeFunc,iter,ISTATE,CCON,CJAC,CLAMBDA,     &
                LValue,GRADP,HESSP,xP,IW,LENIW,RW,LENRW,iuser,RUSER, &
                IFAIL)
    !call ComputeHess(xP(1:Penalty%nx),Penalty%nx,model,LValue,GRAD,Hess,iuser,ruser)
    if (Penalty%nx1>0) then
      x(Penalty%x_index1)    = xP(Penalty%xP_index1)
      Grad(Penalty%x_index1) = GradP(Penalty%xP_index1)
      Hess(Penalty%x_index1,Penalty%x_index1) = HessP(Penalty%xP_index1,Penalty%xP_index1)
      Hess(Penalty%x_index1,Penalty%x_index2) = HessP(Penalty%xP_index1,Penalty%xP_index2)
      Hess(Penalty%x_index2,Penalty%x_index1) = HessP(Penalty%xP_index2,Penalty%xP_index1)
    end if
    x(Penalty%x_index2) = xP(Penalty%xP_index2)
    Grad(Penalty%x_index2) = GradP(Penalty%xP_index2)
    Hess(Penalty%x_index2,Penalty%x_index2) = HessP(Penalty%xP_index2,Penalty%xP_index2)

  end if

  deallocate(A,BL,BU,CLAMBDA,ISTATE,CJAC)
  deallocate(CCON)
  deallocate(xP)
  if (NAG_SAVE==1) then
    if (model==1) then
      close(101)
      close(103)
    elseif (model>=2) then
      close(102)
      close(104)
    end if
  end if

  ierr = 0   ! no error in subroutine
end subroutine MaximizePenalizedLikelihood

!------------------------------------------------------------------------------
! subroutine LikeFunc
!
!------------------------------------------------------------------------------
subroutine LikeFunc(mode,nx,x,L,GradL,nstate,iuser,ruser)
  use nrtype
  use GlobalModule, only : Penalty
  implicit none
  integer(i4b), intent(inout) :: mode 
  integer(i4b), intent(in)    :: nx,nstate,iuser(*)
  real(dp), intent(in)        :: x(nx),ruser(*)
  real(dp), intent(out)       :: L
  real(dp), intent(inout)     :: GradL(nx)
  integer(i4b)                :: model

  model = iuser(1)
#if USE_MPI<=1
  if (model==1) then
    ! Logit model
    if (Penalty%method==0) then
      call Logit(mode,nx,x,L,GradL,nstate,iuser,ruser)
    else
      ! Penalized logit
      call PenalizedLikelihood(mode,nx,x,L,GradL,nstate,iuser,ruser)
    end if
  else if (model==2) then
    if (Penalty%method==0) then
      ! Mixed logit
      call MixedLogit(mode,nx,x,L,GradL,nstate,iuser,ruser)
    else 
      ! Penalized mixed logit
      call PenalizedLikelihood(mode,nx,x,L,GradL,nstate,iuser,ruser)
    end if
  end if
#elif USE_MPI==2
  if (model==1) then
    call Logit_master(mode,nx,x,L,GradL,nstate,iuser,ruser)
  else if (model==2) then
    call MixedLogit_master(mode,nx,x,L,GradL,nstate,iuser,ruser)
  else if (model==3) then
    call PenalizedLikelihood_master(mode,nx,x,L,GradL,nstate,iuser,ruser)
  end if
#endif
end subroutine LikeFunc
 
! MixedLogit model
subroutine MixedLogit(mode,nx,x,L,GradL,nstate,iuser,ruser)
  use nrtype
  use GlobalModule, only : b,CDiag,COffDiag,xData, &
                           K,J,N,RCDIM,              &
                           iXData,iFree,IntRule
  implicit none
! x       =  (nx x 1) free parameters
! b       =  (K x 1)     b coefficients
!
! xData    = (K x J x N) data. data for option chosen is in xData(:,J,:)
! iXData%b = (K x 1)     indexes of xData that are multiplied by b
! iXData%e = (K x 1)     indexes of xData that are multiplied by random 
!                        coefficients
! IntRule.flag = 1   Taylor expand around 0
!              = 2   Laplace approximation
!              = 3   quadrature rule
! IntRule.nodes    = (nAll x RCDim)
! IntRule.weights  = (nAll x 1)
! IntRule.nAll = (1 x 1) number of quadrature nodes
!
! Revision history
! 19sep2012 LN   edit to work with E04WDF
! 06jul2012 LN   translate from MixedLogit.m

  integer(i4b), intent(inout)        :: mode
  integer(i4b), intent(in)           :: nx
  real(dp),     intent(in)           :: x(:)
  real(dp),     intent(out)          :: L
  real(dp),     intent(inout)        :: GradL(:)
  integer(i4b), intent(in)           :: nstate
  integer(i4b), intent(in)           :: iuser(*)
  real(dp),     intent(in)           :: ruser(*)

  integer(i4b)   :: h
  integer(i4b)   :: iVec(1),iTemp
  integer(i4b)   :: row,col,iC
  real(dp), allocatable :: C(:,:)
  real(dp)              :: Prob
  real(dp), allocatable :: GradB(:),GradC(:)

  allocate(C(RCDIM,RCDIM))
  allocate(GradB(K))
  allocate(GradC(RCDIM*(RCDIM+1)/2))
  ! Map x to structural (b,C)
  ! [b,C,CDiag,COffDiag] = MapToStructural
  call MapToStructural(x,RCDim,b,CDiag,COffDiag,C)

  L = 0.0d0
  GradL = 0.0d0

  do h=1,N
    ! P     = probability of choosing option J
    ! GradB = gradient of P w.r.t. b
    ! GradC = gradient of P w.r.t. C
    !% GradC = Grad( C(1,1) )
    !%       = Grad( C(2,1) )
    !%       = Grad( C(2,2) )
    !%       = Grad( C(3,1) ) ...
    call ComputeIntegral(reshape(xData(:,:,h),(/K,J/)),b,C,J,RCDim,mode,Prob,GradB,GradC)

    ! Gradient w.r.t. b
    if (mode>0) then
      GradL(iFree%bX) = GradL(iFree%bX) - GradB(iFree%b)/Prob
      do row=1,RCDim
      do col=1,row
        iC = row*(row-1)/2+col
      
        if (row==col) then
          if (any(iFree%CDiag==row)) then
            iVec = pack(iFree%CDiagX,iFree%CDiag==row)
            GradL(iVec) = GradL(iVec) - GradC(iC)/Prob
          end if
        else
          iTemp = (row-1)*(row-2)/2+col  
          if (any(iFree%COffDiag==iTemp)) then
            iVec = pack(iFree%COffDiagX,iFree%COffDiag==iTemp)
            GradL(iVec) = GradL(iVec) - GradC(iC)/Prob
          end if
        end if
        
      end do
      end do
    end if ! if (mode>0)

    L  = L - dlog(Prob)
  
  end do !  do h=1,N

  deallocate(C,GradB,GradC)
end subroutine MixedLogit

! PenalizedLikelihood:
!  compute objective function for penalized likelihood for model == iuser(1)
subroutine PenalizedLikelihood(mode,nxP,xP,LP,GradLP,nstate,iuser,ruser)
  use nrtype
  use GlobalModule, only : Penalty,N
  implicit none
! xP       =  (nxP x 1) free parameters
!
! Revision history
! 19sep2012 LN   created subroutine 

  integer(i4b), intent(inout)        :: mode
  integer(i4b), intent(in)           :: nxP
  real(dp),     intent(in)           :: xP(:)
  real(dp),     intent(out)          :: LP
  real(dp),     intent(inout)        :: GradLP(:)
  integer(i4b), intent(in)           :: nstate
  integer(i4b), intent(in)           :: iuser(*)
  real(dp),     intent(in)           :: ruser(*)

  ! Method 1 : max L(x1,x2) - P(x2)    subject to      x2 = x2Plus - x2Minus
  !                                                    x2Plus >= 0
  !                                                    x2Minus >= 0
  !                                                    P(x2) = LASSO penalty
  real(dp)                           :: x(Penalty%nx)
  real(dp)                           :: xPlus(Penalty%nx2)
  real(dp)                           :: xMinus(Penalty%nx2)
  real(dp)                           :: L,GradL(Penalty%nx)
  real(dp)                           :: P,GradP(Penalty%nxP)
  integer(i4b)                       :: model
  model = iuser(1)

  x = 0.0d0
  if (Penalty%nx1>0) then
    x(Penalty%x_index1) = xP(Penalty%xP_index1)
  end if

  if (Penalty%nx2>0) then
    x(Penalty%x_index2) = xP(Penalty%xP_index2)
  end if
  
  if (model==1) then 
    ! Logit model
    call Logit(mode,Penalty%nx,x,L,GradL,nstate,iuser,ruser)
  elseif (model==2) then
    ! MixedLogit
    call MixedLogit(mode,Penalty%nx,x,L,GradL,nstate,iuser,ruser)
  end if

  xPlus  = xP(Penalty%xP_index_xPlus)
  xMinus = xP(Penalty%xP_index_xMinus)

  call PenaltyFunction(xPlus,xMinus,P,GradP)

  LP = L + real(N)*P
  if (mode>0) then
    GradLP = 0.d0
    ! elements of GradLP corresponding to x1
    if (Penalty%nx1>0) then
      GradLP(Penalty%xP_index1) = GradL(Penalty%x_index1)
    end if

    ! elements of GradLP corresponding to x2
    GradLP(Penalty%xP_index2) = GradL(Penalty%x_index2)

    ! elements of GradLP corresponding to xPlus
    ! elements of GradLP corresponding to xMinus
    GradLP(Penalty%xP_index_xPlus)  = real(N)*GradP(Penalty%xPlus_index)
    GradLP(Penalty%xP_index_xMinus) = real(N)*GradP(Penalty%xMinus_index)
  end if
end subroutine PenalizedLikelihood

subroutine PenaltyFunction(xPlus,xMinus,P,GradP)
  use nrtype
  use GlobalModule, only : Penalty
  implicit none
  real(dp), intent(in) :: xPlus(:)
  real(dp), intent(in) :: xMinus(:)
  real(dp), intent(out) :: P
  real(dp), intent(out) :: GradP(:)

  if (Penalty%method==1) then
    ! LASSO penalty
    !P = Penalty%lambda*(sum(abs(xPlus)) + sum(abs(xMinus)))
    !GradP(Penalty%xPlus_index) = Penalty%lambda*sign(1.0d0,xPlus)
    !GradP(Penalty%xMinus_index) = Penalty%lambda*sign(1.0d0,xMinus)
    P = Penalty%lambda*(sum(xPlus) + sum(xMinus))
    GradP(Penalty%xPlus_index) = Penalty%lambda
    GradP(Penalty%xMinus_index) = Penalty%lambda
  end if

end subroutine PenaltyFunction

! Logit model
subroutine Logit(mode,nx,x,L,GradL,nstate,iuser,ruser)
  use nrtype
  use GlobalModule, only : b,CDiag,COffDiag,xData, &
                           K,J,N,RCDIM,            &
                           iXData,iFree
  implicit none
! x       =  (nx x 1) free parameters
! b       =  (K x 1)     b coefficients
!
! xData    = (K x J x N) data. data for option chosen is in xData(:,J,:)
! iXData%b = (K x 1)     indexes of xData that are multiplied by b
! iXData%e = (K x 1)     indexes of xData that are multiplied by random 
!                        coefficients
!
! Revision history
! 19sep2012 LN   adapt from MixedLogit
! 19sep2012 LN   edit to work with E04WDF
! 06jul2012 LN   translate from MixedLogit.m

  integer(i4b), intent(inout)        :: mode
  integer(i4b), intent(in)           :: nx
  real(dp),     intent(in)           :: x(:)
  real(dp),     intent(out)          :: L
  real(dp),     intent(inout)        :: GradL(:)
  integer(i4b), intent(in)           :: nstate
  integer(i4b), intent(in)           :: iuser(*)
  real(dp),     intent(in)           :: ruser(*)

  integer(i4b)   :: h
  integer(i4b)   :: iVec(1),iTemp
  integer(i4b)   :: row,col,iC
  real(dp)       :: C(RCDim,RCDim)
  real(dp)       :: Prob,GradB(K),GradC(RCDim*(RCDim+1)/2)

  ! Map x to structural b
  b(iFree%b)=x(iFree%bX)		

  L = 0.0d0
  GradL = 0.0d0

  do h=1,N
    ! P     = probability of choosing option J
    ! GradB = gradient of P w.r.t. b
    ! GradC = gradient of P w.r.t. C
    !% GradC = Grad( C(1,1) )
    !%       = Grad( C(2,1) )
    !%       = Grad( C(2,2) )
    !%       = Grad( C(3,1) ) ...
    call ComputeLogitProbability(reshape(xData(:,:,h),(/K,J/)),b,J,mode,Prob,GradB)

    ! Gradient w.r.t. b
    if (mode>0) then
      GradL(iFree%bX) = GradL(iFree%bX) - GradB(iFree%b)/Prob
    end if ! if (mode==2)

    L  = L - dlog(Prob)
  
  end do !  do h=1,N

end subroutine Logit

subroutine TestGrad(LikeFunc,nx,x,nstate,iuser,ruser)
  use nrtype
  use OutputModule, only : MakeFullFileName
  implicit none
  integer(i4b), intent(in) :: nx,nstate
  real(dp),     intent(in) :: x(:)
  integer(i4b), intent(in) :: iuser(:)
  real(dp),     intent(in) :: ruser(:)

  interface
    subroutine LikeFunc(mode,n,x,L,GradL,nstate,iuser,ruser)
      use nrtype
      integer(i4b), intent(inout) :: mode
      integer(i4b), intent(in)    :: n,nstate,iuser(*)
      real(dp),     intent(in)    :: x(n),ruser(*)
      real(dp),     intent(out)   :: L
      real(dp),     intent(inout) :: GradL(n)
    end subroutine LikeFunc
  end interface

  integer(i4b)  :: mode
  real(dp)      :: L0,GradL0(nx)
  real(dp)      :: L1,L2,GradL1(nx),DummyGrad(nx)
  integer(i4b)  :: i1
  real(dp)      :: x1(nx),x2(nx)

  mode=2
  call LikeFunc(mode,nx,x,L0,GradL0,nstate,iuser,ruser)
  
  open(UNIT   = 1033, &
       File   = MakeFullFileName('TestGrad.txt'), &
       action = 'write')
  GradL1 = 0.0d0
  do i1=1,nx
    x1=x
    x2=x
    x1(i1) = x(i1) + 1.0d-6
    x2(i1) = 2.0*x(i1) - x1(i1)
    mode=0
    call LikeFunc(mode,nx,x1,L1,DummyGrad,nstate,iuser,ruser)
    call LikeFunc(mode,nx,x2,L2,DummyGrad,nstate,iuser,ruser)
    GradL1(i1) = (L1-L2)/(x1(i1)-x2(i1))

    write(1033,1033) i1,x(i1),GradL0(i1),GradL1(i1)
    1033 format(i5,3d25.16)
  end do
  close(1033)
end subroutine TestGrad

!------------------------------------------------------------------------------
! subroutine ComputeHess(xFree,nFree,model,L,GradL,Hess,iuser,ruser)
!------------------------------------------------------------------------------
subroutine ComputeHess(xFree,nFree,model,L,GradL,Hess,iuser,ruser)
  use nrtype
  implicit none
  integer(i4b), intent(in)  :: nFree
  real(dp),     intent(in)  :: xFree(nFree)
  integer(i4b), intent(in)  :: model
  real(dp),     intent(out) :: L
  real(dp),     intent(out) :: GradL(nFree)
  real(dp),     intent(out) :: Hess(nFree,nFree)
  integer(i4b), intent(in)  :: iuser(:)
  real(dp),     intent(in)  :: ruser(:)

  integer(i4b)              :: i1,i2
  integer(i4b)              :: mode,nstate
  real(dp)                  :: x1(nFree),x2(nFree)
  real(dp)                  :: L1,L2
  real(dp)                  :: GradL1(nFree),GradL2(nFree)
  real(dp)                  :: h

  Hess=0.0d0
  mode=1
  nstate=1
  h=1e-6
  if (model==1) then
    call Logit(mode,nFree,xFree,L,GradL,nstate,iuser,ruser)
  else if (model==2) then
    call MixedLogit(mode,nFree,xFree,L,GradL,nstate,iuser,ruser)
  end if
  do i1=1,nFree
    x1=xFree
    x2=xFree
    x1(i1) = xFree(i1) + h
    x2(i1) = 2.0d0*xFree(i1) - x1(i1)
    if (model==1) then
      call Logit(mode,nFree,x1,L1,GradL1,nstate,iuser,ruser)
      call Logit(mode,nFree,x2,L2,GradL2,nstate,iuser,ruser)
    else if (model==2) then
      call MixedLogit(mode,nFree,x1,L1,GradL1,nstate,iuser,ruser)
      call MixedLogit(mode,nFree,x2,L2,GradL2,nstate,iuser,ruser)
    end if
    Hess(:,i1) = (GradL1 - GradL2)/(x1(i1)-x2(i1))
  end do
  Hess = 0.5d0 * (Hess+transpose(Hess))
end subroutine ComputeHess

subroutine MinimizeBIC(model,x,LValue,Grad,Hess,Stats)
  use IFPORT
  use nrtype
  use GlobalModule, only : ResultStructure,Penalty,N
  use OutputModule, only : ComputeStats,SavePenaltyOutputs,SaveBICResults
  implicit none
  integer(i4b), intent(in) :: model
  real(dp),     intent(inout) :: x(:)
  real(dp),     intent(out)   :: LValue
  real(dp),     intent(out)   :: Grad(:)
  real(dp),     intent(out)   :: Hess(:,:)
  type(ResultStructure), intent(inout) :: Stats
  integer(i4b)              :: ierr
  real(dp), allocatable :: x1(:,:),L1(:),Grad1(:,:),Hess1(:,:,:)
  type(ResultStructure), allocatable :: stats1(:)
  real(dp), allocatable :: BIC1(:)
  integer(i4b)          :: i1,MinIndex(1)
  
  ! User INPUTS TO ComputeHess
  REAL(DP)                  :: RUSER(1)
  INTEGER(I4B)              :: iuser(1)

  allocate(x1(size(x,1),Penalty%nLambda))
  allocate(L1(Penalty%nLambda))
  allocate(Grad1(size(x,1),Penalty%nLambda))
  allocate(Hess1(size(x,1),size(x,1),Penalty%nLambda))
  allocate(Stats1(Penalty%nLambda))
  allocate(BIC1(Penalty%nLambda))
  x1    = 0.0d0
  L1    = 0.0d0
  Grad1 = 0.0d0
  Hess1 = 0.0d0
  BIC1  = 0.0d0

  do i1=1,Penalty%nLambda
    Penalty%lambda = Penalty%VectorLambda(i1)
    ierr = 0
    stats1(i1)%lambda = Penalty%lambda
    stats1(i1)%EstimationTime = timef()
    
    call MaximizeLikelihood(model,x,LValue,Grad,Hess,ierr)
    stats1(i1)%EstimationTime = timef()
    call ComputeStats(x,LValue,Grad,Hess,N,Stats1(i1))
    call SavePenaltyOutputs(i1,model,x,LValue,Grad,Hess,Stats1(i1))
    x1(:,i1)      = x
    L1(i1)        = LValue
    Grad1(:,i1)   = Grad
    Hess1(:,:,i1) = Hess
    BIC1(i1)      = Stats1(i1)%BIC
  end do

  call SaveBICResults(model,stats1)
  MinIndex = minloc(BIC1)
  x        = x1(:,MinIndex(1))
  iuser(1) = model
  ruser    = 0.0d0
  call ComputeHess(x,Penalty%nx,model,LValue,Grad,Hess,iuser,ruser)

  !reshape(x1(:,MinIndex),(/size(x,1)/))
  !LValue   = L1(MinIndex(1))
  !Grad     = Grad1(:,MinIndex(1))
  !Hess     = Hess1(:,:,MinIndex(1))
  stats    = stats1(MinIndex(1))

  deallocate(x1,L1,Grad1,Hess1,Stats1,BIC1)
end subroutine MinimizeBIC

subroutine RunMonteCarlo(IMC1,IMC2,model)
  use nrtype
  use GlobalModule, only : K,N,J,RCDIM,iFree,           &
                           ResultStructure
  use DataModule, only       : CreateData
  use OutputModule, only     : SaveMCOutputs
  implicit none
  integer(i4b), intent(in) :: IMC1,IMC2,model
  integer(i4b)             :: IMC
  real(dp), allocatable    :: x(:),xTrue(:),MCX(:,:),MCLambda(:)
  real(dp), allocatable    :: Grad(:),Hess(:,:)
  real(dp)                 :: LVALUE
  type(ResultStructure)    :: stats
  integer(i4b)             :: NMC

  NMC = IMC2-IMC1+1
  call DefineFreeParameters(model)
  allocate(x(iFree%NALL),xTrue(iFree%NALL))
  allocate(MCX(iFree%NALL,NMC+1))
  allocate(MCLambda(NMC))

  call ComputeInitialGuess(model,x)
  xTrue = x
  MCX(:,NMC+1) = xTrue

  allocate(Grad(iFree%NALL),Hess(iFree%NALL,iFree%NALL))
  
  do iMC=IMC1,IMC2
    call CreateData(N,J,K,RCDIM,iMC)
    !call RescaleData

    ! x = parameters
    ! 
    call MinimizeBIC(model,x,LVALUE,Grad,Hess,Stats)
    MCX(:,IMC-IMC1+1) = x
    MCLambda(IMC-IMC1+1) = stats%lambda
    call SaveMCOutputs(model,MCX,MCLambda,IMC)
  end do
  deallocate(x,Grad,Hess)
  deallocate(MCX,MCLambda)
end subroutine RunMonteCarlo
    
end module LikelihoodModule

