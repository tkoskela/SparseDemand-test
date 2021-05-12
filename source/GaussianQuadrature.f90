! Gaussian-type quadrature rules with pre-assigned nodes.
!
! Based on `gaussq.f` from go ('Golden Oldies') netlib library
! (editor Eric Grosse, see https://www.netlib.org/go/ for more details).
!
! The (private) procedures `gaussq`, `solve`, `class` and `gausq2` are direct
! translations of corresponding procedures from netlib go library to a more
! 'modern' Fortran style.
!
! The public subroutine `GaussianQuadratureRule` wraps `gaussq` to provide a
! more convenient interface, with assumed-shape semantics for passing arrays
! and arguments relevant to only some rule types made optional. Further the
! module exposes integer constants `legendre`, `chebyshev1`, `hermite`,
! `jacobi` and `laguerre` for setting the `kind` type code argument - see
! documentation of `GaussianQuadratureRule` for more details.

module GaussianQuadrature

  use ConstantsModule

  implicit none

  ! Enumeration of `kind` parameter defining quadrature rule to use

  ! Legendre quadrature
  integer(i4b), parameter :: legendre = 1
  ! Chebyshev quadrature of the first kind 
  integer(i4b), parameter :: chebyshev1 = 2 
  ! Chebyshev quadrature of the second kind
  integer(i4b), parameter :: chebyshev2 = 3
  ! Hermite quadrature
  integer(i4b), parameter :: hermite = 4
  ! Jacobi quadrature
  integer(i4b), parameter :: jacobi = 5
  ! Generalized Laguerre quadrature
  integer(i4b), parameter :: laguerre = 6

  private
  public legendre, chebyshev1, chebyshev2, hermite, jacobi, laguerre, &
    GaussianQuadratureRule

  contains

    subroutine GaussianQuadratureRule(kind, nodes, weights, alpha, beta, endpts)

      ! Compute the nodes and weights for a Gaussian-type quadrature rule.
      !
      ! A Gaussian-type quadrature rule is used when one wishes to approximate 
      !
      !  integral (x = a to b)  f(x) * w(x) * dx
      !
      ! by
      !
      !  sum (j = 1 to n) weights(j) * f(nodes(j))
      !
      ! for a given integrand function `f`, weight function `w`, integration
      ! interval `[a, b]` and number of points `n` to use in approximation.
      !
      ! Arguments:
      !
      !   kind: An integer between 1 and 6 giving the type of quadrature rule:
      !     legendre = 1: 
      !       Legendre quadrature, `w(x) = 1`, `a = -1`, `b = 1`.
      !     chebyshev1 = 2:  
      !       Chebyshev quadrature of the first kind, `w(x) = 1/sqrt(1 - x*x)`,
      !       `a = -1`, `b = +1`.
      !     chebyshev2 = 3:  
      !       Chebyshev quadrature of the second kind, `w(x) = sqrt(1 - x*x)`,
      !       `a = -1`, `b = +1`.
      !     hermite = 4:  
      !       Hermite quadrature, `w(x) = exp(-x*x)`, `a = -infinity`, 
      !       `b = +infinity`.
      !     jacobi = 5:  
      !       Jacobi quadrature, `w(x) = (1 - x)**alpha * (1 + x)**beta`,
      !       `a = -1`, `b= +1`, `alpha` and `beta` should be greater than -1,
      !       `kind = chebyshev1` and `kind = chebyshev2` are special cases.
      !     laguerre = 6: 
      !       Generalized Laguerre quadrature, `w(x) = exp(-x) * x**alpha`,
      !       `a = 0`, `b = +infinity`, `alpha` should be greater than -1.
      !
      !  nodes: Real array to output nodes to evaluate integrand at. Size `n`
      !    should match size of `weights` and determines number of points.
      !
      !  weights: Real array to output weights for terms in summation. Size `n`
      !    should match size of `nodes` and determines number of points.
      !
      !  alpha: Optional, real. Parameter for weight function when using Jacobi
      !    and Laguerre type rules.
      !
      !  beta: Optional, real. Parameter for weight function when using Jacobi 
      !    type rule.
      !
      !  endpts: Optional, real array. Values for any fixed endpoint nodes.
      !    Normally not required unless the left or right end-point (or both)
      !    of the integration interval is required to be a node (this is called 
      !    Gauss-Radau or Gauss-Lobatto quadrature). If specified should either
      !    be a size 1 array specifying left endpoint, or a size 2 array
      !    specifying left and right endpoints in that order.

      integer(i4b), intent(in) :: kind
      real(dp), intent(out) :: nodes(:), weights(:)
      real(dp), intent(in), optional :: alpha, beta
      real(dp), intent(in), optional :: endpts(:)
      real(dp) :: a_alpha, a_beta
      real(dp) :: a_endpts(2)
      integer(i4b) :: n, kpts
      real(dp), allocatable :: b(:)

      if (present(alpha)) then
        if (kind /= jacobi .and. kind /= laguerre) then
          stop "'alpha' should only be specified for Jacobi and Laguerre rules"
        end if
        a_alpha = alpha
      else
        a_alpha = 0.0_dp
      end if

      if (present(beta)) then
        if (kind /= jacobi) then
          stop "'alpha' should only be specified for Jacobi rule"
        end if
        a_beta = beta
      else
        a_beta = 0.0_dp
      end if

      if (present(endpts)) then
        kpts = size(endpts)
        if (kpts == 1) then
          a_endpts = [endpts(1), 0.0_dp]
        else if (kpts == 2) then
          a_endpts = endpts
        else
          stop "'endpts' should be of size 1 or 2" 
        end if
      else
        kpts = 0
        a_endpts = 0.0_dp
      end if

      n = size(nodes)
      if (size(weights) /= n) stop "Sizes of 'nodes' and 'weights' do not match"

      allocate(b(n))

      call gaussq(kind, n, a_alpha, a_beta, kpts, a_endpts, b, nodes, weights)

    end subroutine

    subroutine gaussq(kind, n, alpha, beta, kpts, endpts, b, t, w)
    
      ! This set of routines computes the nodes `t(j)` and weights `w(j)` for
      ! Gaussian-type quadrature rules with pre-assigned nodes. These are used 
      ! when one wishes to approximate 
      !
      !  integral (x = a to b)  f(x) * w(x) * dx
      !
      ! by
      !
      !  sum (j = 1 to n) w(j) * f(t(j))
      !
      ! (note `w(x)` and `w(j)` have no connection with each other). Here 
      ! `w(x)` is one of six possible non-negative weight functions (listed 
      ! below), and `f(x)` is the function to be integrated. Gaussian 
      ! quadrature is particularly useful on infinite intervals (with 
      ! appropriate weight functions), since then other techniques often fail.
      !
      ! Associated with each weight function `w(x)` is a set of orthogonal 
      ! polynomials. The nodes `t(j)` are just the zeros of the proper n-th 
      ! degree polynomial.
      !
      ! Input parameters (all real numbers are in double precision)
      !
      ! kind: An integer between 1 and 6 giving the type of quadrature rule:
      !
      !   kind = 1:  Legendre quadrature, 
      !              w(x) = 1 on (-1, 1)
      !   kind = 2:  Chebyshev quadrature of the first kind
      !              w(x) = 1/sqrt(1 - x*x) on (-1, +1)
      !   kind = 3:  Chebyshev quadrature of the second kind
      !              w(x) = sqrt(1 - x*x) on (-1, 1)
      !   kind = 4:  Hermite quadrature, 
      !              w(x) = exp(-x*x) on (-infinity, +infinity)
      !   kind = 5:  Jacobi quadrature, 
      !              w(x) = (1-x)**alpha * (1+x)**beta on (-1, 1), 
      !              alpha, beta greater than -1
      !              note: kind = 2 and 3 are a special case of this.
      !   kind = 6:  generalized Laguerre quadrature, 
      !              w(x) = exp(-x) * x**alpha on (0, +infinity), 
      !              alpha greather than -1
      !
      ! n: The number of points used for the quadrature rule.
      !
      ! alpha: Real parameter used only for Gauss-Jacobi and Gauss-Laguerre 
      !   quadrature (otherwise use 0).
      !
      ! beta: Real parameter used only for Gauss-Jacobi quadrature (otherwise 
      !   use 0).
      !
      ! kpts: Integer, normally 0, unless the left or right end-point (or both)
      !  of the interval is required to be a node (this is called Gauss-Radau 
      !  or Gauss-Lobatto quadrature). Then kpts is the number of fixed 
      !  endpoints (1 or 2).
      !
      ! endpts: Real array of length 2. Contains the values of any fixed 
      !   endpoints, if kpts = 1 or 2.
      !
      ! b: Real scratch array of length n.
      !
      ! Output parameters (both double precision arrays of length n)
      !
      ! t: Will contain the desired nodes.
      ! w: Will contain the desired weights w(j).
      !
      ! Underflow may sometimes occur, but is harmless.
      !
      ! References 
      !
      ! 1.  Golub, G. H., and Welsch, J. H., "Calculation of Gaussian 
      !     quadrature rules", Mathematics of Computation 23 (April, 1969), 
      !     pp. 221-230.
      ! 2.  Golub, G. H., "Some modified matrix eigenvalue problems", SIAM 
      !     Review 15 (April, 1973), pp. 318-334 (section 7).
      ! 3.  Stroud and Secrest, Gaussian quadrature formulas, Prentice-Hall, 
      !     Englewood Cliffs, N.J., 1966.
      !
      ! Original version 20 Jan 1975 from Stanford
      ! Modified 21 Dec 1983 by Eric Grosse
      ! Modified to use more modern Fortran style on 12 May 2021 by Matt Graham

      integer(i4b), intent(in) :: kind
      integer(i4b), intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: beta
      integer (i4b), intent(in) :: kpts
      real(dp), intent(in) :: endpts(2)
      real(dp), intent(out) :: b(n)
      real(dp), intent(out) :: t(n)
      real(dp), intent(out) :: w(n)

      integer(i4b) :: i, ierr
      real(dp) :: muzero, t1, gam

      call class(kind, n, alpha, beta, b, t, muzero)

      ! The matrix of coefficients is assumed to be symmetric. The array `t` 
      ! contains the diagonal elements, the array `b` the off-diagonal elements. 
      ! Make appropriate changes in the lower right 2 by 2 submatrix.

      if (kpts == 1) then
        ! if kpts=1, only `t(n)` must be changed
        t(n) = solve(endpts(1), n, t, b) * b(n-1)**2 + endpts(1)
      else if (kpts == 2)  then
        ! if kpts=2, `t(n)` and `b(n-1)` must be recomputed
        gam = solve(endpts(1), n, t, b)
        t1 = (endpts(1) - endpts(2)) / (solve(endpts(2), n, t, b) - gam)
        b(n - 1) = dsqrt(t1)
        t(n) = endpts(1) + gam * t1
      end if 

      ! Note that the indices of the elements of b run from 1 to n-1 and thus 
      ! the value of b(n) is arbitrary. Now compute the eigenvalues of the 
      ! symmetric tridiagonal matrix, which has been modified as necessary. The 
      ! method used is a ql-type method with origin shifting.

      w(1) = 1.0_dp
      do i = 2, n
        w(i) = 0.0_dp
      end do

      call gausq2(n, t, b, w, ierr)

      do i = 1, n
        w(i) = muzero * w(i) * w(i)
      end do

    end subroutine

    function solve(shift, n, a, b) result(delta_n)

      ! This procedure performs elimination to solve for the `n`th component of 
      ! the solution `delta` to the equation
      !
      !   (j_n - shift * identity) * delta  = e_n,
      !
      ! where `e_n` is the vector of all zeroes except for 1 in the `n`th 
      ! position.
      !
      ! The matrix `j_n` is symmetric tridiagonal, with diagonal elements `a`,
      ! and off-diagonal elements `b`. This equation must be solved to obtain 
      ! the appropriate changes in the lower 2 by 2 submatrix of coefficients 
      ! for orthogonal polynomials.

      real(dp), intent(in) :: shift
      integer(i4b), intent(in) :: n
      real(dp), intent(in) :: a(n), b(n)
      real(dp) :: delta_n, alpha
      integer(i4b) :: i

      alpha = a(1) - shift
      do i = 2, n - 1
        alpha = a(i) - shift - b(i - 1)**2 / alpha
      end do
      delta_n = 1.0_dp / alpha

    end function

    subroutine class(kind, n, alpha, beta, b, a, muzero)
    
      ! This procedure supplies the coefficients `a(j)`, `b(j)` of the 
      ! recurrence relation
      !
      !   b(j) * p(j)(x) = (x - a(j)) * p(j - 1)(x) - b(j - 1) * p(j -2)(x)
      !
      ! for the various classical (normalized) orthogonal polynomials, and the 
      ! zero-th moment
      !
      !   muzero = integral(x = -infinity to infinity) w(x) * dx
      !
      ! of the given polynomial's weight function `w(x)`. Since the polynomials
      ! are orthonormalized, the tridiagonal matrix is guaranteed to be 
      ! symmetric.
      !
      ! The input parameter `alpha` is used only for Laguerre and Jacobi 
      ! polynomials, and the parameter `beta` is used only for Jacobi 
      ! polynomials. The Laguerre and Jacobi polynomials require the `gamma` 
      ! function.

      integer(i4b), intent(in) :: kind
      integer(i4b), intent(in) :: n 
      real(dp), intent(in) :: alpha
      real(dp), intent(in) :: beta
      real(dp), intent(out) :: b(n)
      real(dp), intent(out) :: a(n)
      real(dp), intent(out) :: muzero

      integer(i4b) :: i
      real(dp) :: abi, a2b2, ab

      select case (kind)

        case (legendre)
          ! Legendre polynomials p(x) on (-1, +1), 
          ! w(x) = 1
          muzero = 2.0_dp
          do i = 1, n - 1
            a(i) = 0.0_dp
            abi = i
            b(i) = abi / dsqrt(4 * abi * abi - 1.0_dp)
          end do
          a(n) = 0.0_dp

        case (chebyshev1)
          ! Chebyshev polynomials of the first kind t(x) on (-1, +1),
          ! w(x) = 1 / sqrt(1 - x*x)
          muzero = pi_d
          do i = 1, n - 1
            a(i) = 0.0_dp
            b(i) = 0.5_dp
          end do
          b(1) = dsqrt(0.5_dp)
          a(n) = 0.0_dp

        case (chebyshev2)
          ! Chebyshev polynomials of the second kind u(x) on (-1, +1),
          ! w(x) = sqrt(1 - x * x)
          muzero = pi_d / 2.0_dp
          do i = 1, n - 1
            a(i) = 0.0_dp
            b(i) = 0.5_dp
          end do
          a(n) = 0.0_dp

        case (hermite)
          ! Hermite polynomials h(x) on (-infinity, +infinity), 
          ! w(x) = exp(-x**2)
          muzero = dsqrt(pi_d)
          do i = 1, n - 1
            a(i) = 0.0_dp
            b(i) = dsqrt(i / 2.0_dp)
          end do
          a(n) = 0.0_dp

        case (jacobi)
          ! Jacobi polynomials p(alpha, beta)(x) on (-1, +1), 
          ! w(x) = (1 - x)**alpha + (1 + x)**beta, 
          ! alpha and beta greater than -1
          ab = alpha + beta
          abi = 2.0_dp + ab
          muzero = 2.0_dp**(ab + 1.0_dp) * dgamma(alpha + 1.0_dp) &
                   * dgamma(beta + 1.0_dp) / dgamma(abi)
          a(1) = (beta - alpha)/abi
          b(1) = dsqrt(4.0_dp * (1.0_dp + alpha) * (1.0_dp + beta) &
                 / ((abi + 1.0_dp) * abi * abi))
          a2b2 = beta * beta - alpha * alpha
          do i = 2, n - 1
            abi = 2.0_dp * i + ab
            a(i) = a2b2 / ((abi - 2.0_dp) * abi)
            b(i) = dsqrt(4.0_dp * i * (i + alpha) * (i + beta) * (i + ab) &
                   / ((abi * abi - 1) * abi * abi))
          end do
          abi = 2.0_dp * n + ab
          a(n) = a2b2 / ((abi - 2.0_dp) * abi)

        case (laguerre)
          ! Laguerre polynomials l(alpha)(x) on (0, +infinity), 
          ! w(x) = exp(-x) * x**alpha, 
          ! alpha greater than -1.
          muzero = dgamma(alpha + 1.0_dp)
          do i = 1, n - 1
            a(i) = 2.0_dp * i - 1.0_dp + alpha
            b(i) = dsqrt(i * (i + alpha))
          end do
          a(n) = 2.0_dp * n - 1 + alpha

      end select 

    end subroutine

    subroutine gausq2(n, d, e, z, ierr)

      ! This subroutine is a translation of an algol procedure, Num. Math. 12, 
      ! 377-383(1968) by Martin and Wilkinson, as modified in Num. Math. 15, 
      ! 450(1970) by Dubrulle. Handbook for Auto. Comp., Vol.II-Linear Algebra,
      ! 241-248(1971). This is a modified version of the 'eispack' routine 
      ! imtql2.
      !
      ! This subroutine finds the eigenvalues and first components of the
      ! eigenvectors of a symmetric tridiagonal matrix by the implicit QL
      ! method.
      !
      ! On input:
      !
      !   n: Order of the matrix.
      !   d: Contains the diagonal elements of the input matrix. Overwritten.
      !   e: Contains the subdiagonal elements of the input matrix in its first
      !      `n-1` positions. `e(n)` is arbitrary. Overwritten.
      !   z: Contains the first row of the identity matrix. Overwritten.
      !
      ! On outputs:
      !
      !   d: Contains the eigenvalues in ascending order. If an error exit is 
      !      made, the eigenvalues are correct but unordered for indices 
      !      `1, 2, ..., ierr - 1`.
      !   z: Contains the first components of the orthonormal eigenvectors of
      !      the symmetric tridiagonal matrix. If an error exit is made, `z` 
      !      contains the eigenvectors associated with the stored eigenvalues.
      !   ierr: Set to zero for normal return and an integer `j` if the `j`th 
      !      eigenvalue has not been determined after 30 iterations.

      integer(i4b), intent(in) :: n
      real(dp), intent(inout) :: d(n)
      real(dp), intent(inout) :: e(n)
      real(dp), intent(inout) :: z(n)
      integer(i4b), intent(out) :: ierr

      integer(i4b) :: i, j, k, l, m, ii, mml
      real(dp) :: b, c, f, g, p, r, s, machep

      ! machep = d1mach(4)
      ! Based on implementation at
      ! https://github.com/certik/fortran-utils/blob/
      !   b43bd24cd421509a5bc6d3b9c3eeae8ce856ed88/src/legacy/amos/d1mach.f90 
      machep = radix(1.0_dp)
      machep = machep**(1 - digits(1.0_dp))

      ierr = 0

      if (n == 1) return

      e(n) = 0.0_dp

      do l = 1, n
        
        j = 0

        do while (j <= 30)
        
          ! look for small sub-diagonal element
          do m = l, n
            if (m == n) exit
            if (dabs(e(m)) <= machep * (dabs(d(m)) + dabs(d(m+1)))) exit
          end do

          p = d(l)
          if (m == l) exit
          if (j == 30) then
            ! set error - no convergence to an eigenvalue after 30 iterations
            ierr = l
            return
          end if
          j = j + 1
          
          ! form shift
          g = (d(l+1) - p) / (2.0_dp * e(l))
          r = dsqrt(g * g + 1.0_dp)
          g = d(m) - p + e(l) / (g + dsign(r, g))
          s = 1.0_dp
          c = 1.0_dp
          p = 0.0_dp
          mml = m - l

          ! for i = m-1 step -1 until l do
          do ii = 1, mml
            i = m - ii
            f = s * e(i)
            b = c * e(i)
            if (dabs(f) < dabs(g)) then
              s = f / g
              r = dsqrt(s * s + 1.0_dp)
              e(i + 1) = g * r
              c = 1.0_dp / r
              s = s * c
            else
              c = g / f
              r = dsqrt(c * c + 1.0_dp)
              e(i + 1) = f * r
              s = 1.0_dp / r
              c = c * s
            end if
            g = d(i + 1) - p
            r = (d(i) - g) * s + 2.0_dp * c * b
            p = s * r
            d(i + 1) = g + p
            g = c * r - b
            ! form first component of vector
            f = z(i + 1)
            z(i + 1) = s * z(i) + c * f
            z(i) = c * z(i) - s * f
          end do

          d(l) = d(l) - p
          e(l) = g
          e(m) = 0.0_dp
        
        end do 

      end do
      
      ! order eigenvalues and eigenvectors
      do ii = 2, n
        i = ii - 1
        k = i
        p = d(i)
        do j = ii, n
          if (d(j) >= p) cycle
          k = j
          p = d(j)
        end do
        if (k == i) cycle
        d(k) = d(i)
        d(i) = p
        p = z(i)
        z(i) = z(k)
        z(k) = p
      end do
  
    end subroutine

end module GaussianQuadrature