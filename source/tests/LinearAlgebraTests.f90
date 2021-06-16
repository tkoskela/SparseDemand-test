! Tests for linear algebra routines

module LinearAlgebraTests

  use ConstantsModule
  use LinearAlgebra

  implicit none

  abstract interface

    subroutine matrix_check_subroutine(mtx)
      import :: dp
      real(dp), intent(in) :: mtx(:, :)
    end subroutine

    subroutine matrix_matrix_check_subroutine(mtx_1, mtx_2)
      import :: dp
      real(dp), intent(in) :: mtx_1(:, :)
      real(dp), intent(in) :: mtx_2(:, :)
    end subroutine

    subroutine matrix_vector_check_subroutine(mtx, vct)
      import :: dp
      real(dp), intent(in) :: mtx(:, :)
      real(dp), intent(in) :: vct(:)
    end subroutine

  end interface

  real(dp), parameter :: abs_tolerance = 1e-12_dp
  real(dp), parameter :: rel_tolerance = 1e-10_dp
  integer(i4b), parameter :: seed_val = 20210525
  integer(i4b), parameter, dimension(4) :: size_1 = [1, 2, 5, 10]
  integer(i4b), parameter, dimension(4) :: size_2 = [1, 2, 5, 10]

  contains

  subroutine SeedRNG(val)
    integer(i4b), intent(in) :: val
    integer :: n
    integer, allocatable :: seed(:)

    call random_seed(size=n)
    allocate(seed(n))
    seed = val
    call random_seed(put=seed)
  end

  elemental function AreClose(x, y) result(are_close)

    real(dp), intent(in) :: x, y
    logical are_close

    are_close = abs(x - y) <= abs_tolerance + rel_tolerance * abs(y)

  end function

  function IsCloseToIdentity(a) result(is_close)

    real(dp), intent(in) :: a(:, :)
    logical :: is_close
    integer(i4b) :: i, j, m, n

    is_close = .True.

    m = size(a, 1)
    n = size(a, 2)

    do i = 1, m
      do j = 1, n
        if (i == j) then
          if (.not. AreClose(a(i, j), 1.0_dp)) then
            is_close = .False.
          end if
        else
          if (.not. AreClose(a(i, j), 0.0_dp)) then
            is_close = .False.
          end if
        end if
      end do
    end do
  end function

  function IsLowerTriangular(a) result(is_lower)
    real(dp), intent(in) :: a(:, :)
    logical :: is_lower
    integer(i4b) :: i, j, m, n
    m = size(a, 1)
    n = size(a, 2)
    is_lower = .True.
    do i = 1, m
      do j = i + 1, n
        if (a(i, j) > 0) then
          is_lower = .False.
        end if
      end do
    end do
  end function

  function IsUpperTriangular(a) result(is_upper)
    real(dp), intent(in) :: a(:, :)
    logical :: is_upper
    integer(i4b) :: i, j, m, n
    m = size(a, 1)
    n = size(a, 2)
    is_upper = .True.
    do i = 1, m
      do j = 1, min(i - 1, n)
        if (a(i, j) > 0) then
          is_upper = .False.
        end if
      end do
    end do
  end function

  subroutine CheckCopyTriangular(a)

    real(dp), intent(in) :: a(:, :)

    integer(i4b) :: m, n
    real(dp), allocatable :: b(:, :)

    m = size(a, 1)
    n = size(a, 2)

    allocate(b(m, n))

    call CopyTriangular(a, b, .True.)

    if (.not. IsLowerTriangular(b)) then
      print "(A)", "Destination array not lower triangular"
    end if

    call CopyTriangular(a, b, .False.)

    if (.not. IsUpperTriangular(b)) then
      print "(A)", "Destination array not upper triangular"
    end if
  end subroutine

  subroutine CheckComputeSVD(a)

    real(dp), intent(in) :: a(:, :)

    integer(i4b) :: m, n, k, i
    real(dp), allocatable :: original_a(:, :), u(:, :), s(:), vt(:, :), s_full(:, :)

    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)

    allocate(original_a(m, n))
    allocate(u(m, m))
    allocate(s(k))
    allocate(vt(n, n))
    allocate(s_full(m, n))

    original_a = a

    call ComputeSVD(a, u, s, vt)

    if (any(a /= original_a)) then
      print "(A)", "Input array changed"
    end if

    if (any(s < 0)) then
      print "(A)", "One or more negative singular values"
    end if

    if (any((s(1:k-1) - s(2:k)) < 0)) then
      print "(A)", "Singular values not in descending order"
    end if

    if (.not. IsCloseToIdentity(matmul(u, transpose(u)))) then
      print "(A)", "Left singular vectors not orthogonal"
    end if

    if (.not. IsCloseToIdentity(matmul(transpose(vt), vt))) then
      print "(A)", "Right singular vectors not orthogonal"
    end if

    s_full = 0.0_dp
    do i = 1, k
      s_full(i, i) = s(i)
    end do

    if (.not. all(AreClose(matmul(matmul(u, s_full), vt), a))) then
      print "(A)", "Product of factors does not equal original matrix"
    end if

  end subroutine

  subroutine CheckComputeLQ(a)

    real(dp), intent(in) :: a(:, :)

    integer(i4b) :: m, n, k
    real(dp), allocatable :: original_a(:, :), l(:, :), q(:, :)

    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)

    allocate(original_a(m, n))
    allocate(l(m, n))
    allocate(q(n, n))

    original_a = a

    call ComputeLQ(a, l, q)

    if (any(a /= original_a)) then
      print "(A)", "Input array changed"
    end if

    if (.not. IsLowerTriangular(l)) then
      print "(A)", "l factor is not lower triangular"
    end if

    if (.not. IsCloseToIdentity(matmul(q, transpose(q)))) then
      print "(A)", "q factor not orthogonal"
    end if

    if (.not. all(AreClose(matmul(l, q), a))) then
      print "(A)", "Product of factors does not equal original matrix"
    end if

  end subroutine

  subroutine CheckComputeQR(a)

    real(dp), intent(in) :: a(:, :)

    integer(i4b) :: m, n, k
    real(dp), allocatable :: original_a(:, :), q(:, :), r(:, :)

    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)

    allocate(original_a(m, n))
    allocate(r(m, n))
    allocate(q(m, m))

    original_a = a

    call ComputeQR(a, q, r)

    if (any(a /= original_a)) then
      print "(A)", "Input array changed"
    end if

    if (.not. IsUpperTriangular(r)) then
      print "(A)", "r factor is not upper triangular"
    end if

    if (.not. IsCloseToIdentity(matmul(q, transpose(q)))) then
      print "(A)", "q factor not orthogonal"
    end if

    if (.not. all(AreClose(matmul(q, r), a))) then
      print "(A)", "Product of factors does not equal original matrix"
    end if

  end subroutine

  subroutine CheckComputeCholesky(a, compute_lower_triangular)

    real(dp), intent(in) :: a(:, :)
    logical, intent(in) :: compute_lower_triangular

    integer(i4b) :: n
    real(dp), allocatable :: original_a(:, :), chol_a(:, :)

    n = size(a, 1)

    allocate(original_a(n, n))
    allocate(chol_a(n, n))

    original_a = a

    call ComputeCholesky(a, chol_a, compute_lower_triangular)

    if (any(a /= original_a)) then
      print "(A)", "Input array changed"
    end if

    ! if computing upper triangular factor transpose prior to doing checks
    if (.not. compute_lower_triangular) then
      chol_a = transpose(chol_a)
    end if

    if (.not. IsLowerTriangular(chol_a)) then
      print "(A)", "Cholesky factor is not triangular"
    end if
    if (.not. all(AreClose(matmul(chol_a, transpose(chol_a)), a))) then
      print "(A)", "Product of factors does not equal original matrix"
    end if

  end subroutine

  subroutine CheckComputeEigenSymmetric(a)

    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable :: original_a(:, :), w1(:), w2(:), w_full(:, :), v(:, :)

    integer(i4b) :: i, n

    n = size(a, 1)
    allocate(original_a(n, n), v(n, n), w1(n), w2(n), w_full(n, n))
    original_a = a

    ! Compute full eigendecomposition
    call ComputeEigenSymmetric(a, w1, v)

    if (any(a /= original_a)) then
      print "(A)", "Input array changed"
    end if

    if (any((w1(2:n) - w1(1:n-1)) < 0)) then
      print "(A)", "Eigenvalues not in ascending order"
    end if

    if (.not. IsCloseToIdentity(matmul(v, transpose(v)))) then
      print "(A)", "Eigenvectors not orthonormal"
    end if

    w_full = 0.0_dp
    do i = 1, n
      w_full(i, i) = w1(i)
    end do

    if (.not. all(AreClose(matmul(matmul(v, w_full), transpose(v)), a))) then
      print "(A)", "Product of factors does not equal original matrix"
    end if

    ! Compute just eigenvalues
    call ComputeEigenSymmetric(a, w2)

    if (any(a /= original_a)) then
      print "(A)", "Input array changed"
    end if

    if (.not. all(AreClose(w1, w2))) then
      print "(A)", "Eigenvalues do not match values from full eigendecomposition"
    end if

  end subroutine

  subroutine CheckInvertTriangular(a, a_is_lower)

    real(dp), intent(in) :: a(:, :)
    logical, intent(in) :: a_is_lower
    real(dp), allocatable :: original_a(:, :), inv_a(:, :)

    integer(i4b) :: n

    n = size(a, 1)
    allocate(original_a(n, n))
    original_a = a
    allocate(inv_a(n, n))

    call InvertTriangular(a, inv_a, a_is_lower)

    if (a_is_lower) then
      if (.not. IsLowerTriangular(inv_a)) then
        print "(A)", "Non-lower-triangular inverse of lower-triangular matrix"
      end if
    else
      if (.not. IsUpperTriangular(inv_a)) then
        print "(A)", "Non-upper-triangular inverse of upper-triangular matrix"
      end if
    end if

    if (any(a /= original_a)) then
      print "(A)", "Input array changed"
    end if

    if (.not. IsCloseToIdentity(matmul(a, inv_a))) then
      print "(A)", "Right multiplication by inverse does not equal identity"
    end if

    if (.not. IsCloseToIdentity(matmul(inv_a, a))) then
      print "(A)", "Left multiplication by inverse does not equal identity"
    end if

  end subroutine

  subroutine CheckInvertSymmetric(a)

    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable :: original_a(:, :), inv_a(:, :)

    integer(i4b) :: n

    n = size(a, 1)
    allocate(original_a(n, n))
    original_a = a
    allocate(inv_a(n, n))

    call InvertSymmetric(a, inv_a)

    if (any(a /= original_a)) then
      print "(A)", "Input array changed"
    end if

    if (.not. all(AreClose(inv_a, (inv_a + transpose(inv_a)) / 2))) then
      print "(A)", "Inverse not symmetric"
    end if

    if (.not. IsCloseToIdentity(matmul(a, inv_a))) then
      print "(A)", "Right multiplication by inverse does not equal identity"
    end if

    if (.not. IsCloseToIdentity(matmul(inv_a, a))) then
      print "(A)", "Left multiplication by inverse does not equal identity"
    end if

  end subroutine

  subroutine CheckSolveMultiple(a, b)

    real(dp), intent(in) :: a(:, :), b(:, :)
    real(dp), allocatable :: original_a(:, :), original_b(:, :), x(:, :)

    integer(i4b) :: m, n

    m = size(a, 1)
    n = size(b, 2)
    allocate(original_a(m, m), original_b(m, n), x(m, n))
    original_a = a
    original_b = b

    call Solve(a, b, x)

    if (any(a /= original_a)) then
      print "(A)", "Input array `a` changed"
    end if
    if (any(b /= original_b)) then
      print "(A)", "Input array `b` changed"
    end if

    if (.not. all(AreClose(matmul(a, x), b))) then
      print "(A)", "Solution `x` does not satisfy  `matmul(a, x) = b`"
    end if

  end subroutine

  subroutine CheckSolveSingle(a, b)

    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), allocatable :: original_a(:, :), original_b(:), x(:)

    integer(i4b) :: n

    n = size(a, 1)
    allocate(original_a(n, n), original_b(n), x(n))
    original_a = a
    original_b = b

    call Solve(a, b, x)

    if (any(a /= original_a)) then
      print "(A)", "Input array `a` changed"
    end if
    if (any(b /= original_b)) then
      print "(A)", "Input array `b` changed"
    end if

    if (.not. all(AreClose(matmul(a, x), b))) then
      print "(A)", "Solution `x` does not satisfy  `matmul(a, x) = b`"
    end if

  end subroutine

  subroutine CheckSolveTriangularMultiple(a, b, a_is_lower, solve_transposed)

    real(dp), intent(in) :: a(:, :), b(:, :)
    logical, intent(in) :: a_is_lower, solve_transposed
    real(dp), allocatable :: original_a(:, :), original_b(:, :), x(:, :)

    integer(i4b) :: m, n

    m = size(a, 1)
    n = size(b, 2)
    allocate(original_a(m, m), original_b(m, n), x(m, n))
    original_a = a
    original_b = b

    call SolveTriangular(a, b, x, a_is_lower, solve_transposed)

    if (any(a /= original_a)) then
      print "(A)", "Input array `a` changed"
    end if
    if (any(b /= original_b)) then
      print "(A)", "Input array `b` changed"
    end if

    if (solve_transposed) then
      if (.not. all(AreClose(matmul(transpose(a), x), b))) then
        print "(A)", "Solution `x` does not satisfy  `matmul(transpose(a), x) = b`"
      end if
    else
      if (.not. all(AreClose(matmul(a, x), b))) then
        print "(A)", "Solution `x` does not satisfy  `matmul(a, x) = b`"
      end if
    end if

  end subroutine

  subroutine CheckSolveTriangularSingle(a, b, a_is_lower, solve_transposed)

    real(dp), intent(in) :: a(:, :), b(:)
    logical, intent(in) :: a_is_lower, solve_transposed
    real(dp), allocatable :: original_a(:, :), original_b(:), x(:)

    integer(i4b) :: n

    n = size(a, 1)
    allocate(original_a(n, n), original_b(n), x(n))
    original_a = a
    original_b = b

    call SolveTriangular(a, b, x, a_is_lower, solve_transposed)

    if (any(a /= original_a)) then
      print "(A)", "Input array `a` changed"
    end if
    if (any(b /= original_b)) then
      print "(A)", "Input array `b` changed"
    end if

    if (solve_transposed) then
      if (.not. all(AreClose(matmul(transpose(a), x), b))) then
        print "(A)", "Solution `x` does not satisfy  `matmul(transpose(a), x) = b`"
      end if
    else
      if (.not. all(AreClose(matmul(a, x), b))) then
        print "(A)", "Solution `x` does not satisfy  `matmul(a, x) = b`"
      end if
    end if

  end subroutine

  subroutine CheckLogAbsDet(a)

    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable :: original_a(:, :), a_sq(:, :)

    real(dp) :: log_abs_det_a, log_abs_det_a_t, log_abs_det_a_sq, trace_a_sq
    integer(i4b) :: i, n

    n = size(a, 1)
    allocate(original_a(n, n), a_sq(n, n))
    original_a = a

    log_abs_det_a = LogAbsDet(a)

    if (any(a /= original_a)) then
      print "(A)", "Input array changed"
    end if

    log_abs_det_a_t = LogAbsDet(transpose(a))

    if (.not. AreClose(log_abs_det_a, log_abs_det_a_t)) then
      print "(A)", "Log absolute determinant not invariant under transposition"
    end if

    a_sq = matmul(a, transpose(a))
    log_abs_det_a_sq = LogAbsDet(a_sq)

    if (.not. AreClose(2 * log_abs_det_a, log_abs_det_a_sq)) then
      print "(A)", "Log absolute determinant not scaling correctly"
    end if

    trace_a_sq = 0.0_dp
    do i = 1, n
      trace_a_sq = trace_a_sq + a_sq(i, i)
    end do

    ! For the positive definite matrix b = matmul(a, transpose(a)) then
    ! LogAbsDet(b) <= Trace(b - identity(n)) = Trace(b) - n
    if (log_abs_det_a_sq > (trace_a_sq - n)) then
      print "(A)", "Log absolute determinant violating trace inequality"
    end if

    if (exp(LogAbsDet(a * 0)) /= 0) then
      print "(A)", "Log absolute determinant of zero matrix not negative infinity"
    end if

  end subroutine

  subroutine DoTestsOnGeneratedRectangularMatrices(check_subroutine, subroutine_name)

    procedure(matrix_check_subroutine) :: check_subroutine

    character(*), intent(in) :: subroutine_name

    integer(i4b) :: i, j, m, n
    real(dp), allocatable :: a(:, :)

    print "(A,A)", "Running tests for ", subroutine_name
    call SeedRNG(seed_val)

    do i = 1, size(size_1)
      m = size_1(i)
      do j = 1, size(size_2)
          n = size_2(j)
          allocate(a(m, n))
          call random_number(a)
          call check_subroutine(a)
          deallocate(a)
      end do
    end do

  end subroutine

  subroutine DoTestsOnGeneratedSquareMatrices(check_subroutine, subroutine_name)

    procedure(matrix_check_subroutine) :: check_subroutine

    character(*), intent(in) :: subroutine_name

    integer(i4b) :: i, n
    real(dp), allocatable :: a(:, :)

    print "(A,A)", "Running tests for ", subroutine_name
    call SeedRNG(seed_val)

    do i = 1, size(size_1)
      n = size_1(i)
      allocate(a(n, n))
      call random_number(a)
      call check_subroutine(a)
      deallocate(a)
    end do

  end subroutine

  subroutine DoTestsOnGeneratedSquareMatrixVectorPair(check_subroutine, subroutine_name)

    procedure(matrix_vector_check_subroutine) :: check_subroutine

    character(*), intent(in) :: subroutine_name

    integer(i4b) :: i, n
    real(dp), allocatable :: a(:, :), b(:)

    print "(A,A)", "Running tests for ", subroutine_name
    call SeedRNG(seed_val)

    do i = 1, size(size_1)
      n = size_1(i)
      allocate(a(n, n), b(n))
      call random_number(a)
      call random_number(b)
      call check_subroutine(a, b)
      deallocate(a, b)
    end do

  end subroutine

  subroutine DoTestsOnGeneratedSquareMatrixMatrixPair(check_subroutine, subroutine_name)

    procedure(matrix_matrix_check_subroutine) :: check_subroutine

    character(*), intent(in) :: subroutine_name

    integer(i4b) :: i, j, m, n
    real(dp), allocatable :: a(:, :), b(:, :)

    print "(A,A)", "Running tests for ", subroutine_name
    call SeedRNG(seed_val)

    do i = 1, size(size_1)
      m = size_1(i)
      do j = 1, size(size_2)
          n = size_2(j)
          allocate(a(m, m), b(m, n))
          call random_number(a)
          call random_number(b)
          call check_subroutine(a, b)
          deallocate(a, b)
      end do
    end do

  end subroutine

  subroutine TestCopyTriangular()

    call DoTestsOnGeneratedRectangularMatrices(CheckCopyTriangular, "CopyTriangular")

  end subroutine

  subroutine TestComputeSVD()

    call DoTestsOnGeneratedRectangularMatrices(CheckComputeSVD, "ComputeSVD")

  end subroutine

  subroutine TestComputeLQ()

    call DoTestsOnGeneratedRectangularMatrices(CheckComputeLQ, "ComputeLQ")

  end subroutine

  subroutine TestComputeQR()

    call DoTestsOnGeneratedRectangularMatrices(CheckComputeQR, "ComputeQR")

  end subroutine

  subroutine TestComputeCholesky()

    call DoTestsOnGeneratedSquareMatrices(check_subroutine, "ComputeCholesky")

    contains

    subroutine check_subroutine(r)

      real(dp), intent(in) :: r(:, :)
      real(dp), allocatable :: a(:, :)

      integer(i4b) :: n

      n = size(r, 1)
      allocate(a(n, n))

      a = matmul(r, transpose(r))
      call CheckComputeCholesky(a, .True.)
      call CheckComputeCholesky(a, .False.)

    end subroutine

  end subroutine

  subroutine TestComputeEigenSymmetric()

    call DoTestsOnGeneratedSquareMatrices(check_subroutine, "ComputeEigenSymmetric")

    contains

    subroutine check_subroutine(r)

      real(dp), intent(in) :: r(:, :)
      real(dp), allocatable :: a(:, :)

      integer(i4b) :: n

      n = size(r, 1)
      allocate(a(n, n))

      a = (r + transpose(r)) / 2
      call CheckComputeEigenSymmetric(a)

    end subroutine

  end subroutine

  subroutine TestInvertTriangular()

    call DoTestsOnGeneratedSquareMatrices(check_subroutine, "InvertTriangular")

    contains

    subroutine check_subroutine(r)

      real(dp), intent(in) :: r(:, :)
      real(dp), allocatable :: a(:, :)

      integer(i4b) :: n

      n = size(r, 1)
      allocate(a(n, n))

      call CopyTriangular(r, a, .True.)
      call CheckInvertTriangular(a, .True.)
      call CopyTriangular(r, a, .False.)
      call CheckInvertTriangular(a, .False.)

    end subroutine

  end subroutine

  subroutine TestInvertSymmetric()

    call DoTestsOnGeneratedSquareMatrices(check_subroutine, "InvertSymmetric")

    contains

    subroutine check_subroutine(r)

      real(dp), intent(in) :: r(:, :)
      real(dp), allocatable :: a(:, :)

      integer(i4b) :: n

      n = size(r, 1)
      allocate(a(n, n))

      a = (r + transpose(r)) / 2
      call CheckInvertSymmetric(a)

    end subroutine

  end subroutine

  subroutine TestSolveMultiple()

    call DoTestsOnGeneratedSquareMatrixMatrixPair(CheckSolveMultiple, "SolveMultiple")

  end subroutine

  subroutine TestSolveSingle()

    call DoTestsOnGeneratedSquareMatrixVectorPair(CheckSolveSingle, "SolveSingle")

  end subroutine

  subroutine TestSolveTriangularMultiple()

    call DoTestsOnGeneratedSquareMatrixMatrixPair(&
      check_subroutine, "SolveTriangularMultiple")

    contains

    subroutine check_subroutine(r, b)

      real(dp), intent(in) :: r(:, :), b(:, :)
      real(dp), allocatable :: a(:, :)

      integer(i4b) :: n

      n = size(r, 1)
      allocate(a(n, n))

      call CopyTriangular(r, a, .True.)
      call CheckSolveTriangularMultiple(a, b, .True., .True.)
      call CheckSolveTriangularMultiple(a, b, .True., .False.)
      call CopyTriangular(r, a, .False.)
      call CheckSolveTriangularMultiple(a, b, .False., .True.)
      call CheckSolveTriangularMultiple(a, b, .False., .False.)

    end subroutine

  end subroutine

  subroutine TestSolveTriangularSingle()

    call DoTestsOnGeneratedSquareMatrixVectorPair(&
      check_subroutine, "SolveTriangularSingle")

    contains

    subroutine check_subroutine(r, b)

      real(dp), intent(in) :: r(:, :), b(:)
      real(dp), allocatable :: a(:, :)

      integer(i4b) :: n

      n = size(r, 1)
      allocate(a(n, n))

      call CopyTriangular(r, a, .True.)
      call CheckSolveTriangularSingle(a, b, .True., .True.)
      call CheckSolveTriangularSingle(a, b, .True., .False.)
      call CopyTriangular(r, a, .False.)
      call CheckSolveTriangularSingle(a, b, .False., .True.)
      call CheckSolveTriangularSingle(a, b, .False., .False.)

    end subroutine

  end subroutine

  subroutine TestLogAbsDet()

    call DoTestsOnGeneratedSquareMatrices(CheckLogAbsDet, "LogAbsDet")

  end subroutine

end module