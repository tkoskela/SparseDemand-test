! Tests for Gaussian-type quadrature rules

module GaussianQuadratureTests

  use ConstantsModule
  use GaussianQuadrature

  implicit none

  real(dp), parameter :: tolerance = 1e-14_dp

  contains

  function RuleNameString(kind) result(string)
    ! Convert rule kind code to corresponding string description

    integer(i4b), intent(in) :: kind
    character(len=25) :: string

    select case (kind)

      case (legendre)
        string = "Legendre"
      case (chebyshev1)
        string = "Chebyshev (first kind)"
      case (chebyshev2)
        string = "Chebyshev (second kind)"
      case (hermite)
        string = "Hermite"
      case (jacobi)
        string = "Jacobi"
      case (laguerre)
        string = "Laguerre"

    end select

  end function

  function MaxAbsDifference(x1, x2) result(diff)
    ! Compute L-infinity norm (maximum absolute difference) between two arrays

    real(dp), intent(in) :: x1(:), x2(:)
    real(dp) :: diff

    diff = maxval(abs(x1 - x2))

  end function

  subroutine TestGaussianQuadratureRule(&
      kind, true_nodes, true_weights, alpha, beta, endpts)
    
    ! Test that computed nodes and weights match reference values

    integer(i4b), intent(in) :: kind
    real(dp), intent(in) :: true_nodes(:), true_weights(:)
    real(dp), intent(in), optional :: alpha, beta, endpts(:)
    real(dp), allocatable :: nodes(:), weights(:)
    integer(i4b) :: num_points
    real(dp) :: true_total_weight

    num_points = size(true_nodes)
    allocate(nodes(num_points))
    allocate(weights(num_points))

    ! Only pass optional arguments if present
    if (present(endpts)) then
      if (present(alpha)) then
        if (present(beta)) then
          call GaussianQuadratureRule(kind, nodes, weights, alpha, beta, endpts)
        else
          call GaussianQuadratureRule(kind, nodes, weights, alpha, endpts=endpts)
        end if
      else
        call GaussianQuadratureRule(kind, nodes, weights, endpts=endpts)
      end if
    else
      if (present(alpha)) then
        if (present(beta)) then
          call GaussianQuadratureRule(kind, nodes, weights, alpha, beta)
        else
          call GaussianQuadratureRule(kind, nodes, weights, alpha)
        end if
      else
        call GaussianQuadratureRule(kind, nodes, weights)
      end if
    end if

    ! All weights should be non-negative
    if (any(weights < 0)) then
      print "(A,A)", "Negative weights in rule: ", RuleNameString(kind)
      print "(A)", "Computed weights"
      print *, weights
    end if

    ! Total weight depends on rule kind and alpha and beta for Jacobi + Laguerre
    select case (kind)
      case (legendre)
        true_total_weight = 2.0_dp
      case (chebyshev1)
        true_total_weight = pi_d
      case (chebyshev2)
        true_total_weight = pi_d / 2
      case (hermite)
        true_total_weight = dsqrt(pi_d)
      case (jacobi)
        if (present(alpha) .and. present(beta)) then
          if (alpha == 0.0_dp .and. beta == 0.0_dp) then
            ! equivalent to Legendre
            true_total_weight = 2.0_dp
          else if (alpha == -0.5_dp .and. beta == -0.5_dp) then
            ! equivalent to Chebyshev-1
            true_total_weight = pi_d
          else if (alpha == 0.5_dp .and. beta == 0.5_dp) then
            ! equivalent to Chebyshev-2
            true_total_weight = pi_d / 2
          else
            stop "Unknown total weight for Jacobi rule for alpha and beta used"
          end if
        else
          stop "alpha and beta must both be specified for Jacobi rule test"
        end if
      case (laguerre)
        if (present(alpha)) then
          if (alpha /= 0.0_dp) then
            stop "Unknown total weight for Laguerre rule for alpha used"
          end if
        end if
        true_total_weight = 1.0_dp
    end select

    ! Total of computed weights should match reference value within tolerance
    if (abs(sum(weights) - true_total_weight) > tolerance) then
      print "(A,A)", "Total weight wrong for rule: ", RuleNameString(kind)
      print "(A,F7.5)", "Computed total weight: ", sum(weights)
      print "(A,F7.5)", "True total weight: ", true_total_weight
    end if

    ! Node values should be within tolerance of reference values
    if (MaxAbsDifference(true_nodes, nodes) > tolerance) then
      print "(A,A)", "Mismatch in nodes for rule: ", RuleNameString(kind)
      print "(A)", "Computed nodes"
      print *, nodes
      print "(A)", "Reference nodes"
      print *, true_nodes
    end if

    ! Weight values should be within tolerance of reference values
    if (MaxAbsDifference(true_weights, weights) > tolerance) then
      print "(A,A)", "Mismatch in weights for rule: ", RuleNameString(kind)
      print "(A)", "Computed weights"
      print *, weights
      print "(A)", "Reference weights"
      print *, true_weights
    end if

  end subroutine

end module

program TestGaussianQuadrature

  use ConstantsModule
  use GaussianQuadrature
  use GaussianQuadratureTests

  implicit None

  ! References values for test cases generated using numpy.polynomial functions
  real(dp), parameter, dimension(1) :: legendre_nodes_1 = &
    [0.00000000000000000_dp]
  real(dp), parameter, dimension(1) :: legendre_weights_1 = &
    [2.00000000000000000_dp]
  real(dp), parameter, dimension(2) :: legendre_nodes_2 = &
    [-0.57735026918962573_dp, 0.57735026918962573_dp]
  real(dp), parameter, dimension(2) :: legendre_weights_2 = &
    [1.00000000000000000_dp, 1.00000000000000000_dp]
  real(dp), parameter, dimension(3) :: legendre_nodes_3 = &
    [-0.77459666924148340_dp, 0.00000000000000000_dp, 0.77459666924148340_dp]
  real(dp), parameter, dimension(3) :: legendre_weights_3 = &
    [0.55555555555555569_dp, 0.88888888888888884_dp, 0.55555555555555569_dp]
  real(dp), parameter, dimension(5) :: legendre_nodes_5 = &
    [-0.90617984593866396_dp, -0.53846931010568311_dp, 0.00000000000000000_dp, &
    0.53846931010568311_dp, 0.90617984593866396_dp]
  real(dp), parameter, dimension(5) :: legendre_weights_5 = &
    [0.23692688505618942_dp, 0.47862867049936619_dp, 0.56888888888888900_dp, &
    0.47862867049936619_dp, 0.23692688505618942_dp]
  real(dp), parameter, dimension(10) :: legendre_nodes_10 = &
    [-0.97390652851717174_dp, -0.86506336668898454_dp, -0.67940956829902444_dp, &
    -0.43339539412924721_dp, -0.14887433898163122_dp, 0.14887433898163122_dp, &
    0.43339539412924721_dp, 0.67940956829902444_dp, 0.86506336668898454_dp, &
    0.97390652851717174_dp]
  real(dp), parameter, dimension(10) :: legendre_weights_10 = &
    [0.06667134430868807_dp, 0.14945134915058036_dp, 0.21908636251598201_dp, &
    0.26926671930999652_dp, 0.29552422471475298_dp, 0.29552422471475298_dp, &
    0.26926671930999652_dp, 0.21908636251598201_dp, 0.14945134915058036_dp, &
    0.06667134430868807_dp]
  real(dp), parameter, dimension(1) :: chebyshev1_nodes_1 = &
    [0.00000000000000006_dp]
  real(dp), parameter, dimension(1) :: chebyshev1_weights_1 = &
    [3.14159265358979312_dp]
  real(dp), parameter, dimension(2) :: chebyshev1_nodes_2 = &
    [-0.70710678118654746_dp, 0.70710678118654757_dp]
  real(dp), parameter, dimension(2) :: chebyshev1_weights_2 = &
    [1.57079632679489656_dp, 1.57079632679489656_dp]
  real(dp), parameter, dimension(3) :: chebyshev1_nodes_3 = &
    [-0.86602540378443871_dp, 0.00000000000000006_dp, 0.86602540378443871_dp]
  real(dp), parameter, dimension(3) :: chebyshev1_weights_3 = &
    [1.04719755119659763_dp, 1.04719755119659763_dp, 1.04719755119659763_dp]
  real(dp), parameter, dimension(5) :: chebyshev1_nodes_5 = &
    [-0.95105651629515353_dp, -0.58778525229247303_dp, 0.00000000000000006_dp, &
    0.58778525229247314_dp, 0.95105651629515353_dp]
  real(dp), parameter, dimension(5) :: chebyshev1_weights_5 = &
    [0.62831853071795862_dp, 0.62831853071795862_dp, 0.62831853071795862_dp, &
    0.62831853071795862_dp, 0.62831853071795862_dp]
  real(dp), parameter, dimension(10) :: chebyshev1_nodes_10 = &
    [-0.98768834059513766_dp, -0.89100652418836779_dp, -0.70710678118654746_dp, &
    -0.45399049973954669_dp, -0.15643446504023059_dp, 0.15643446504023092_dp, &
    0.45399049973954680_dp, 0.70710678118654757_dp, 0.89100652418836790_dp, &
    0.98768834059513777_dp]
  real(dp), parameter, dimension(10) :: chebyshev1_weights_10 = &
    [0.31415926535897931_dp, 0.31415926535897931_dp, 0.31415926535897931_dp, &
    0.31415926535897931_dp, 0.31415926535897931_dp, 0.31415926535897931_dp, &
    0.31415926535897931_dp, 0.31415926535897931_dp, 0.31415926535897931_dp, &
    0.31415926535897931_dp]
  real(dp), parameter, dimension(1) :: chebyshev2_nodes_1 = &
    [0.00000000000000006_dp]
  real(dp), parameter, dimension(1) :: chebyshev2_weights_1 = &
    [1.57079632679489656_dp]
  real(dp), parameter, dimension(2) :: chebyshev2_nodes_2 = &
    [-0.49999999999999978_dp, 0.50000000000000011_dp]
  real(dp), parameter, dimension(2) :: chebyshev2_weights_2 = &
    [0.78539816339744839_dp, 0.78539816339744817_dp]
  real(dp), parameter, dimension(3) :: chebyshev2_nodes_3 = &
    [-0.70710678118654746_dp, 0.00000000000000006_dp, 0.70710678118654757_dp]
  real(dp), parameter, dimension(3) :: chebyshev2_weights_3 = &
    [0.39269908169872425_dp, 0.78539816339744828_dp, 0.39269908169872403_dp]
  real(dp), parameter, dimension(5) :: chebyshev2_nodes_5 = &
    [-0.86602540378443871_dp, -0.49999999999999978_dp, 0.00000000000000006_dp, &
    0.50000000000000011_dp, 0.86602540378443871_dp]
  real(dp), parameter, dimension(5) :: chebyshev2_weights_5 = &
    [0.13089969389957468_dp, 0.39269908169872420_dp, 0.52359877559829882_dp, &
    0.39269908169872408_dp, 0.13089969389957468_dp]
  real(dp), parameter, dimension(10) :: chebyshev2_nodes_10 = &
    [-0.95949297361449737_dp, -0.84125353283118109_dp, -0.65486073394528499_dp, &
    -0.41541501300188632_dp, -0.14231483827328500_dp, 0.14231483827328512_dp, &
    0.41541501300188644_dp, 0.65486073394528510_dp, 0.84125353283118121_dp, &
    0.95949297361449737_dp]
  real(dp), parameter, dimension(10) :: chebyshev2_weights_10 = &
    [0.02266894250185884_dp, 0.08347854093418908_dp, 0.16312217745481658_dp, &
    0.23631356020348732_dp, 0.27981494230309656_dp, 0.27981494230309650_dp, &
    0.23631356020348732_dp, 0.16312217745481658_dp, 0.08347854093418901_dp, &
    0.02266894250185884_dp]
  real(dp), parameter, dimension(1) :: hermite_nodes_1 = &
    [0.00000000000000000_dp]
  real(dp), parameter, dimension(1) :: hermite_weights_1 = &
    [1.77245385090551588_dp]
  real(dp), parameter, dimension(2) :: hermite_nodes_2 = &
    [-0.70710678118654746_dp, 0.70710678118654746_dp]
  real(dp), parameter, dimension(2) :: hermite_weights_2 = &
    [0.88622692545275794_dp, 0.88622692545275794_dp]
  real(dp), parameter, dimension(3) :: hermite_nodes_3 = &
    [-1.22474487139158894_dp, 0.00000000000000000_dp, 1.22474487139158894_dp]
  real(dp), parameter, dimension(3) :: hermite_weights_3 = &
    [0.29540897515091941_dp, 1.18163590060367718_dp, 0.29540897515091941_dp]
  real(dp), parameter, dimension(5) :: hermite_nodes_5 = &
    [-2.02018287045608558_dp, -0.95857246461381851_dp, 0.00000000000000000_dp, &
    0.95857246461381851_dp, 2.02018287045608558_dp]
  real(dp), parameter, dimension(5) :: hermite_weights_5 = &
    [0.01995324205904592_dp, 0.39361932315224107_dp, 0.94530872048294179_dp, &
    0.39361932315224107_dp, 0.01995324205904592_dp]
  real(dp), parameter, dimension(10) :: hermite_nodes_10 = &
    [-3.43615911883773739_dp, -2.53273167423278966_dp, -1.75668364929988163_dp, &
    -1.03661082978951358_dp, -0.34290132722370459_dp, 0.34290132722370459_dp, &
    1.03661082978951358_dp, 1.75668364929988163_dp, 2.53273167423278966_dp, &
    3.43615911883773739_dp]
  real(dp), parameter, dimension(10) :: hermite_weights_10 = &
    [0.00000764043285523_dp, 0.00134364574678123_dp, 0.03387439445548111_dp, &
    0.24013861108231471_dp, 0.61086263373532579_dp, 0.61086263373532579_dp, &
    0.24013861108231471_dp, 0.03387439445548111_dp, 0.00134364574678123_dp, &
    0.00000764043285523_dp]
  real(dp), parameter, dimension(1) :: laguerre_nodes_1 = &
    [1.00000000000000000_dp]
  real(dp), parameter, dimension(1) :: laguerre_weights_1 = &
    [1.00000000000000000_dp]
  real(dp), parameter, dimension(2) :: laguerre_nodes_2 = &
    [0.58578643762690497_dp, 3.41421356237309492_dp]
  real(dp), parameter, dimension(2) :: laguerre_weights_2 = &
    [0.85355339059327373_dp, 0.14644660940672624_dp]
  real(dp), parameter, dimension(3) :: laguerre_nodes_3 = &
    [0.41577455678347908_dp, 2.29428036027904181_dp, 6.28994508293747945_dp]
  real(dp), parameter, dimension(3) :: laguerre_weights_3 = &
    [0.71109300992917290_dp, 0.27851773356924098_dp, 0.01038925650158613_dp]
  real(dp), parameter, dimension(5) :: laguerre_nodes_5 = &
    [0.26356031971814087_dp, 1.41340305910651676_dp, 3.59642577104072192_dp, &
    7.08581000585883736_dp, 12.64080084427578221_dp]
  real(dp), parameter, dimension(5) :: laguerre_weights_5 = &
    [0.52175561058280850_dp, 0.39866681108317598_dp, 0.07594244968170769_dp, &
    0.00361175867992205_dp, 0.00002336997238578_dp]
  real(dp), parameter, dimension(10) :: laguerre_nodes_10 = &
    [0.13779347054049260_dp, 0.72945454950317101_dp, 1.80834290174031587_dp, &
    3.40143369785489957_dp, 5.55249614006380376_dp, 8.33015274676449735_dp, &
    11.84378583790006623_dp, 16.27925783137810356_dp, 21.99658581198076135_dp, &
    29.92069701227389089_dp]
  real(dp), parameter, dimension(10) :: laguerre_weights_10 = &
    [0.30844111576501732_dp, 0.40111992915527611_dp, 0.21806828761180960_dp, &
    0.06208745609867777_dp, 0.00950151697518110_dp, 0.00075300838858754_dp, &
    0.00002825923349600_dp, 0.00000042493139850_dp, 0.00000000183956482_dp, &
    0.00000000000099118_dp]

  ! Reference values for Radau (fixed left endpoint) Legendre rules taken from
  ! tabulations at https://mathworld.wolfram.com/RadauQuadrature.html
  real(dp), parameter, dimension(1) :: radau_nodes_1 = [-1.0_dp]
  real(dp), parameter, dimension(1) :: radau_weights_1 = [2.0_dp]
  real(dp), parameter, dimension(2) :: radau_nodes_2 = [-1.0_dp, 1.0_dp / 3]
  real(dp), parameter, dimension(2) :: radau_weights_2 = [0.5_dp, 1.5_dp]
  real(dp), parameter, dimension(3) :: radau_nodes_3 = &
    [-1.0_dp, (1 - sqrt(6.0_dp)) / 5, (1 + sqrt(6.0_dp)) / 5]
  real(dp), parameter, dimension(3) :: radau_weights_3 = &
    [2.0_dp / 9, (16 + sqrt(6.0_dp)) / 18, (16 - sqrt(6.0_dp)) / 18]

  ! Reference values for Lobatto (fixed endpoints) Legendre rules taken from
  ! tabulations at https://mathworld.wolfram.com/LobattoQuadrature.html
  real(dp), parameter, dimension(3) :: lobatto_nodes_3 = &
    [-1.0_dp, 0.0_dp, 1.0_dp]
  real(dp), parameter, dimension(3) :: lobatto_weights_3 = &
    [1.0_dp / 3, 4.0_dp / 3, 1.0_dp / 3]
  real(dp), parameter, dimension(4) :: lobatto_nodes_4 = &
    [-1.0_dp, -sqrt(5.0_dp) / 5, sqrt(5.0_dp) / 5, 1.0_dp]
  real(dp), parameter, dimension(4) :: lobatto_weights_4 = &
    [1.0_dp / 6, 5.0_dp / 6, 5.0_dp / 6, 1.0_dp / 6]
  real(dp), parameter, dimension(5) :: lobatto_nodes_5 = &
    [-1.0_dp, -sqrt(21.0_dp) / 7, 0.0_dp, sqrt(21.0_dp) / 7, 1.0_dp]
  real(dp), parameter, dimension(5) :: lobatto_weights_5 = &
    [1.0_dp / 10, 49.0_dp / 90, 32.0_dp / 45, 49.0_dp / 90, 1.0_dp / 10]

  call TestGaussianQuadratureRule(&
    legendre, legendre_nodes_1, legendre_weights_1)
  call TestGaussianQuadratureRule(&
    legendre, legendre_nodes_2, legendre_weights_2)
  call TestGaussianQuadratureRule(&
    legendre, legendre_nodes_3, legendre_weights_3)
  call TestGaussianQuadratureRule(&
    legendre, legendre_nodes_5, legendre_weights_5)
  call TestGaussianQuadratureRule(&
    legendre, legendre_nodes_10, legendre_weights_10)

  call TestGaussianQuadratureRule(&
    chebyshev1, chebyshev1_nodes_1, chebyshev1_weights_1)
  call TestGaussianQuadratureRule(&
    chebyshev1, chebyshev1_nodes_2, chebyshev1_weights_2)
  call TestGaussianQuadratureRule(&
    chebyshev1, chebyshev1_nodes_3, chebyshev1_weights_3)
  call TestGaussianQuadratureRule(&
    chebyshev1, chebyshev1_nodes_5, chebyshev1_weights_5)
  call TestGaussianQuadratureRule(&
    chebyshev1, chebyshev1_nodes_10, chebyshev1_weights_10)

  call TestGaussianQuadratureRule(&
    chebyshev2, chebyshev2_nodes_1, chebyshev2_weights_1)
  call TestGaussianQuadratureRule(&
    chebyshev2, chebyshev2_nodes_2, chebyshev2_weights_2)
  call TestGaussianQuadratureRule(&
    chebyshev2, chebyshev2_nodes_3, chebyshev2_weights_3)
  call TestGaussianQuadratureRule(&
    chebyshev2, chebyshev2_nodes_5, chebyshev2_weights_5)
  call TestGaussianQuadratureRule(&
    chebyshev2, chebyshev2_nodes_10, chebyshev2_weights_10)

  call TestGaussianQuadratureRule(&
    hermite, hermite_nodes_1, hermite_weights_1)
  call TestGaussianQuadratureRule(&
    hermite, hermite_nodes_2, hermite_weights_2)
  call TestGaussianQuadratureRule(&
    hermite, hermite_nodes_3, hermite_weights_3)
  call TestGaussianQuadratureRule(&
    hermite, hermite_nodes_5, hermite_weights_5)
  call TestGaussianQuadratureRule(&
    hermite, hermite_nodes_10, hermite_weights_10)

  call TestGaussianQuadratureRule(&
    laguerre, laguerre_nodes_1, laguerre_weights_1)
  call TestGaussianQuadratureRule(&
    laguerre, laguerre_nodes_2, laguerre_weights_2)
  call TestGaussianQuadratureRule(&
    laguerre, laguerre_nodes_3, laguerre_weights_3)
  call TestGaussianQuadratureRule(&
    laguerre, laguerre_nodes_5, laguerre_weights_5)
  call TestGaussianQuadratureRule(&
    laguerre, laguerre_nodes_10, laguerre_weights_10)


  ! Jacobi rule with alpha = 0, beta = 0 equivalent to Legendre rule
  call TestGaussianQuadratureRule(&
    jacobi, legendre_nodes_1, legendre_weights_1, 0.0_dp, 0.0_dp)
  call TestGaussianQuadratureRule(&
    jacobi, legendre_nodes_2, legendre_weights_2, 0.0_dp, 0.0_dp)
  call TestGaussianQuadratureRule(&
    jacobi, legendre_nodes_3, legendre_weights_3, 0.0_dp, 0.0_dp)
  call TestGaussianQuadratureRule(&
    jacobi, legendre_nodes_5, legendre_weights_5, 0.0_dp, 0.0_dp)
  call TestGaussianQuadratureRule(&
    jacobi, legendre_nodes_10, legendre_weights_10, 0.0_dp, 0.0_dp)

  ! Jacobi rule with alpha = -0.5, beta = -0.5 equivalent to Chebyshev-1 rule
  call TestGaussianQuadratureRule(&
    jacobi, chebyshev1_nodes_1, chebyshev1_weights_1, -0.5_dp, -0.5_dp)
  call TestGaussianQuadratureRule(&
    jacobi, chebyshev1_nodes_2, chebyshev1_weights_2, -0.5_dp, -0.5_dp)
  call TestGaussianQuadratureRule(&
    jacobi, chebyshev1_nodes_3, chebyshev1_weights_3, -0.5_dp, -0.5_dp)
  call TestGaussianQuadratureRule(&
    jacobi, chebyshev1_nodes_5, chebyshev1_weights_5, -0.5_dp, -0.5_dp)
  call TestGaussianQuadratureRule(&
    jacobi, chebyshev1_nodes_10, chebyshev1_weights_10, -0.5_dp, -0.5_dp)

  ! Jacobi rule with alpha = 0.5, beta = 0.5 equivalent to Chebyshev-2 rule
  call TestGaussianQuadratureRule(&
    jacobi, chebyshev2_nodes_1, chebyshev2_weights_1, 0.5_dp, 0.5_dp)
  call TestGaussianQuadratureRule(&
    jacobi, chebyshev2_nodes_2, chebyshev2_weights_2, 0.5_dp, 0.5_dp)
  call TestGaussianQuadratureRule(&
    jacobi, chebyshev2_nodes_3, chebyshev2_weights_3, 0.5_dp, 0.5_dp)
  call TestGaussianQuadratureRule(&
    jacobi, chebyshev2_nodes_5, chebyshev2_weights_5, 0.5_dp, 0.5_dp)
  call TestGaussianQuadratureRule(&
    jacobi, chebyshev2_nodes_10, chebyshev2_weights_10, 0.5_dp, 0.5_dp)

  call TestGaussianQuadratureRule(&
    legendre, radau_nodes_1, radau_weights_1, endpts=[-1.0_dp])
  call TestGaussianQuadratureRule(&
    legendre, radau_nodes_2, radau_weights_2, endpts=[-1.0_dp])
  call TestGaussianQuadratureRule(&
    legendre, radau_nodes_3, radau_weights_3, endpts=[-1.0_dp])

  call TestGaussianQuadratureRule(&
    legendre, lobatto_nodes_3, lobatto_weights_3, endpts=[-1.0_dp, 1.0_dp])
  call TestGaussianQuadratureRule(&
    legendre, lobatto_nodes_4, lobatto_weights_4, endpts=[-1.0_dp, 1.0_dp])
  call TestGaussianQuadratureRule(&
    legendre, lobatto_nodes_5, lobatto_weights_5, endpts=[-1.0_dp, 1.0_dp])

end program