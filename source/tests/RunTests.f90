program RunTests

  use GaussianQuadratureTests
  use LinearAlgebraTests

  implicit none

  print "(A)", "--------------------------------------"
  print "(A)", "Running Gaussian quadrature rule tests"
  print "(A)", "--------------------------------------"
  call TestLegendreRule()
  call TestChebyshev1Rule()
  call TestChebyshev2Rule()
  call TestHermiteRule()
  call TestLaguerreRule()
  call TestJacobiRule()
  call TestLegendreRadauRule
  call TestLegendreLobattoRule()


  print "(A)", "--------------------------------------"
  print "(A)", "Running linear algebra tests"
  print "(A)", "--------------------------------------"
  call TestCopyTriangular()
  call TestComputeSVD()
  call TestComputeLQ()
  call TestComputeQR()
  call TestComputeCholesky()
  call TestComputeEigenSymmetric()
  call TestInvertTriangular()
  call TestInvertSymmetric()
  call TestSolveMultiple()
  call TestSolveSingle()
  call TestSolveTriangularMultiple()
  call TestSolveTriangularSingle()
  call TestLogAbsDet()

end program