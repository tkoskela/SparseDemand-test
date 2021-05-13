program RunTests

    use GaussianQuadratureTests

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

end program