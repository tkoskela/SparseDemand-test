subroutine LogDensityFunc(e,F,GradF_e,GradF_MuE,GradF_InvC)
! e          = (K x 1)
! F          = (1 x 1) log density evaluated at e
! GradF_e    = (K x 1) gradient of log density w.r.t. e
! GradF_MuE  = (K x 1) gradient of log density w.r.t. MuE
! GradF_InvC = (K x K) gradient of log density w.r.t. InvC 
!                           
! Revision history
! 09DEC2012 LN translated to Fortran from matlab LogDensity.m
  use nrtype
  use GlobalModule, only : parms
  implicit none
  real(dp), intent(in)     :: e(:)
  real(dp), intent(out)    :: F
  real(dp), intent(out)    :: GradF_e(:)
  real(dp), intent(out)    :: GradF_MuE(:)
  real(dp), intent(out)    :: GradF_InvC(:,:)

  real(dp)                 :: z(parms%K)
  integer(i4b)             :: i1,i2
  real(dp)                 :: temp2(parms%K,parms%K)
  real(dp)                 :: temp1(parms%K)

  ! workspace to compute determinant and matrix inverse
  real(dp)                 :: D
  integer(i4b)             :: ID
  integer(i4b)             :: ifail
  integer(i4b)             :: iPivot(parms%K)
  real(dp)                 :: rCond,errBound

  external F03BFF  ! determinant of the cholesky c
  external F04BAF  ! compute inverse of matrix using LU decomp

  z = matmul(parms%InvCSig,e - parms%MuE)

  ! F03BFF: Compute determinant of InvC
  ! log(det(InvC)) = 0.5* (log(D) + ID*log(2.0))
  call F03BFF(parms%K,parms%InvC,parms%K,D,ID,ifail)

  ! Ln(L) = -0.5*z'*z -0.5*K*log(2*pi) + log(det(InvC))
  F     = -0.5d0*dot_product(z,z) - 0.5d0*real(parms%K)*log(2.0d0*pi) &
          + 0.5d0*log(D) + 0.5d0*real(ID)*log(2.0d0)

  GradF_e    = -transpose(parms%InvC)*z
  GradF_MuE  = -GradF_e
  GradF_InvC = 0.0d0
  do i1=1,parms%K
  do i2=1,i1
    temp2 = 0.0d0
    temp2(:,i2) = temp2(:,i2) + transpose(parms%InvC(i1,:));
    temp2(i1,:) = temp2(i1,:) + parms%InvC(i2,:);
    temp1       = matmul(temp2,e-parms%MuE)
    GradF%InvC(i2,i1) = -0.5*matmul(transpose(e-parms%MuE),temp1)
  end
  end
  
  temp2 = 0.0d0
  do i1=1,parms%K
    temp2(i1,i1) = 1.0d0
  end do
  ! compute temp2 = inv(parms%InvC)
  call F04BAF(parms%K,parms%K,parms%InvC,parms%K,iPivot,temp2, &
              parms%K,rCond,errBound,ifail)

  GradF%InvC = GradF%InvC + transpose(temp2)

end subroutine LogDensityFunc