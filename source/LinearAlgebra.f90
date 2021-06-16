! Linear algebra routines for computing matrix factorizations, inverses, solving systems
!
! Most procedures simply wrap LAPACK routines with a more user-friendly interface.

module LinearAlgebra

  use ConstantsModule
  use ieee_arithmetic

  implicit none

  private
  public &
    ComputeLQ, ComputeQR, ComputeSVD, ComputeCholesky, ComputeEigenSymmetric, &
    InvertTriangular, InvertSymmetric, Solve, SolveTriangular, LogAbsDet, CopyTriangular

  interface Solve
    ! Solve a general linear system with a vector or matrix right-hand-side
    module procedure SolveSingle, SolveMultiple
  end interface

  interface SolveTriangular
    ! Solve a triangular linear system with a vector or matrix right-hand-side
    module procedure SolveTriangularSingle, SolveTriangularMultiple
  end interface

  contains

  subroutine CopyTriangular(source_array, dest_array, copy_lower)
    ! Copy the elements above or below and including diagonal from one array to another.
    !
    ! If `copy_lower = .True.` copies the elements below and including the diagonal of
    ! `source_array` to `dest_array`, which should be of the same size, with all other
    ! elements of `dest_array` set to zero.
    !
    ! Else if `copy_lower = .False.` copies the elements above and including the
    ! diagonal of `source_array` to `dest_array`, which should be of the same size, with
    ! all other elements of `dest_array` set to zero.
    !
    ! Arguments:
    !
    ! source_array: Array containing values to copy from.
    ! dest_array: Array which on exit will contain copied values.
    ! copy_lower: Whether to copy lower (True) or upper (False) triangular elements.

    real(dp), intent(in) :: source_array(:, :)
    logical, intent(in) :: copy_lower
    real(dp), intent(out) :: dest_array(:, :)
    integer (i4b) :: i, j, m, n
    m = size(source_array, 1)
    n = size(source_array, 2)
    if ((size(dest_array, 1) /= m) .or. (size(dest_array, 2) /= n)) then
        stop "Source and destination arrays must be the same size"
    end if
    dest_array = 0.0_dp
    do i = 1, m
      if (copy_lower) then
        do j = 1, min(i, n)
          dest_array(i, j) = source_array(i, j)
        end do
      else
        do j = i, n
          dest_array(i, j) = source_array(i, j)
        end do
      end if
    end do
  end subroutine

  subroutine ComputeLQ(a, l, q)
    ! Compute LQ factorization of a square or rectangular matrix.
    !
    ! For a 2D input array `a` of size `(m, n)` computes the `(m, n)` lower triangular /
    ! trapezoidal matrix `l` and `(n, n)` orthogonal matrix `q` such that
    ! `a = matmul(l, q)`.
    !
    ! Uses the LAPACK subroutines dgelqf and dorglq.
    !
    ! Arguments:
    !
    ! a: Array containing matrix to be factorized, of size `(m, n)`.
    ! l: Array which on exit will contain the size `(m, n)` lower triangular /
    !  trapezoidal matrix factor.
    ! q: Array which on exit will contain the size `(n, n)` orthogonal matrix factor.

    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: l(:, :),  q(:, :)

    integer(i4b) :: k, m, n, lwork, info

    real(dp), allocatable :: a_copy(:, :), tau(:), work(:)

    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)

    if (size(l, 1) /= m .or. size(l, 2) /= n) then
        stop "`l` must be of size `(m, n)` where `m = size(a, 1)` and `n = size(a, 2)`"
    end if

    if (size(q, 1) /= n .or. size(q, 2) /= n) then
        stop "`q` must be of size `(n, n)` where `n = size(a, 2)`"
    end if

    ! `a` argument overwritten by dgelqf therefore create copy
    allocate(a_copy(m, n))
    a_copy = a

    allocate(tau(k))

    ! Workspace query to get optimal lwork
    allocate(work(1))
    call dgelqf(m, n, a_copy, m, tau, work, -1, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))

    call dgelqf(m, n, a_copy, m, tau, work, lwork, info)

    if (info < 0) then
      stop "An argument to dgelqf has an illegal value"
    end if

    ! On exit from dgelqf elements on and below diagonal of a_copy correspond to
    ! non-zero elements of l
    call CopyTriangular(a_copy, l, .True.)

    deallocate(work)
    allocate(work(1))
    if (n > m) then
      ! a is 'wide' with n > m and smaller than q matrix with size (n, n) therefore copy
      ! values written in dgelqf to first m rows of q and call dorglq on q
      q(1:m, 1:n) = a_copy
      ! dorgqr workspace query to get optimal lwork
      call dorglq(n, n, k, q, n, tau, work, -1, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
      ! compute orthogonal q matrix from reflectors outputted in a_copy by dgelqf
      call dorglq(n, n, k, q, n, tau, work, lwork, info)
    else
      ! a is 'tall' with m >= n and larger than q matrix with size (n, n) therefore
      ! call dorgqr directly on values written in dgelqf to a_copy
      ! dorglq workspace query to get optimal lwork
      call dorglq(n, n, k, a_copy, m, tau, work, -1, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
      ! compute orthogonal q matrix from reflectors outputted in a_copy by dgelqf
      call dorglq(n, n, k, a_copy, m, tau, work, lwork, info)
      q = a_copy(1:n, 1:n)
    end if

    if (info < 0) then
      stop "An argument to dorglq has an illegal value"
    end if

  end subroutine

  subroutine ComputeQR(a, q, r)
    ! Compute QR factorization of a square or rectangular matrix.
    !
    ! For a 2D input array `a` of size `(m, n)` computes the `(m, m)` orthogonal matrix
    ! `q` and `(m, n)` upper triangular / trapezoidal matrix `r` such that
    ! `a = matmul(q, r)`.
    !
    ! Uses the LAPACK subroutines dgeqrf and dorgqr.
    !
    ! Arguments:
    !
    ! a: Array containing matrix to be factorized, of size `(m, n)`.
    ! q: Array which on exit will contain the size `(m, m)` orthogonal matrix factor.
    ! r: Array which on exit will contain the size `(m, n)` upper triangular /
    !  trapezoidal matrix factor.

    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: q(:, :),  r(:, :)

    integer(i4b) :: k, m, n, lwork, info

    real(dp), allocatable :: a_copy(:, :), tau(:), work(:)

    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)

    if (size(q, 1) /= m .or. size(q, 2) /= m) then
      stop "`q` must be of size `(m, m)` where `m = size(a, 1)`"
    end if

    if (size(r, 1) /= m .or. size(r, 2) /= n) then
      stop "`r` must be of size `(m, n)` where `m = size(a, 1)` and `n = size(a, 2)`"
    end if

    ! `a` argument overwritten by dgelqf therefore create copy
    allocate(a_copy(m, n))
    a_copy = a

    allocate(tau(k))

    ! dgeqrf workspace query to get optimal lwork
    allocate(work(1))
    call dgeqrf(m, n, a_copy, m, tau, work, -1, info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))

    ! Compute upper triangular factor r and reflectors defining orthogonal matrix q
    call dgeqrf(m, n, a_copy, m, tau, work, lwork, info)

    if (info < 0) then
      stop "An argument to dgeqrf has an illegal value"
    end if

    ! On exit from dgeqrf elements on and above diagonal of a_copy correspond to
    ! non-zero elements of r
    call CopyTriangular(a_copy, r, .False.)

    deallocate(work)
    allocate(work(1))
    if (m > n) then
      ! a is 'tall' with m > n and smaller than q matrix with size (m, m) therefore copy
      ! values written in dgeqrf to first n columns of q and call dorgqr on q
      q(1:m, 1:n) = a_copy
      ! dorgqr workspace query to get optimal lwork
      call dorgqr(m, m, k, q, m, tau, work, -1, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
      ! compute orthogonal q matrix from reflectors outputted in a_copy by dorgqr
      call dorgqr(m, m, k, q, m, tau, work, lwork, info)
    else
      ! a is 'wide' with n >= m and larger than q matrix with size (m, m) therefore
      ! call dorgqr directly on values written in dgeqrf to a_copy
      ! dorgqr workspace query to get optimal lwork
      call dorgqr(m, m, k, a_copy, m, tau, work, -1, info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
      ! compute orthogonal q matrix from reflectors outputted in a_copy by dorgqr
      call dorgqr(m, m, k, a_copy, m, tau, work, lwork, info)
      q = a_copy(1:m, 1:m)
    end if

    if (info < 0) then
      stop "An argument to dorgqr has an illegal value"
    end if

  end subroutine

  subroutine ComputeSVD(a, u, s, vt, info)
    ! Compute singular value decomposition (SVD) of a square or rectangular matrix.
    !
    ! For a matrix `a` of size `(m, n)` computes the `(m, m)` orthogonal matrix `u`,
    ! `min(m, n)` array `s` and `(n, n)` orthogonal matrix `vt` such that `s` contains
    ! the non-zero singular values of `a`, the columns of `u` the left singular vectors
    ! of `a`, the rows of `vt` the right singular vectors of `a` and such that
    ! `a = matmul(matmul(u, s_full), vt)` where `s_full` is a `(m, n)` rectangular
    ! matrix with the elements of `s` along the first `min(m, n)` elements of its
    ! diagonal and zeros elsewhere.
    !
    ! Uses the LAPACK subroutine dgesvd.
    !
    ! Arguments:
    !
    ! a: Array containing matrix to be factorized, of size `(m, n)`.
    ! u: Array which on exit will contain a size `(m, m)` orthogonal matrix with
    !   columns corresponding to the left singular vectors of `a`.
    ! s: Array which on exit will contain the `min(m, n)` non-zero singular values of
    !   `a`.
    ! vt: Array which on exit will contain a size `(n, n)` orthogonal matrix with
    !   rows corresponding to the right singular vectors of `a`.
    ! info: Optional. If present, any non-zero output from the `info` diagnostic
    !  argument (corresponding to failure conditions) to the LAPACK routine dgesvd is
    !  outputted here to allow appropriate handling in downstream code. Otherwise if
    !  not present, execution is halted and an error message printed on failure
    !  conditions.

    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: u(:, :), s(:), vt(:, :)
    integer(i4b), intent(out), optional :: info

    integer(i4b) :: k, m, n, lwork, local_info
    real(dp), allocatable :: a_copy(:, :), work(:)

    m = size(a, 1)
    n = size(a, 2)
    k = min(m, n)

    if (size(u, 1) /= m .or. size(u, 2) /= m) then
        stop "`u` must be of size (m, m) where m = size(a, 1)"
    end if

    if (size(s, 1) /= k) then
        stop "`s` must be of size `min(m, n)` where `[m, n] = shape(a)`"
    end if

    if (size(vt, 1) /= n .or. size(vt, 2) /= n) then
        stop "`vt` must be of size `(n, n)` where `n = size(a, 2)`"
    end if

    ! `a` argument overwritten by dgesvd therefore create copy
    allocate(a_copy(m, n))
    a_copy = a

    ! dgesvd workspace query to get optimal lwork
    allocate(work(1))
    call dgesvd("A", "A", m, n, a_copy, m, s, u, m, vt, n, work, -1, local_info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))

    ! compute SVD of a
    call dgesvd("A", "A", m, n, a_copy, m, s, u, m, vt, n, work, lwork, local_info)

    if (local_info < 0) then
      stop "An argument to dgesvd has an illegal value"
    else if (present(info)) then
      info = local_info
    else if (local_info > 0) then
      stop "Iteration in dgesvd did not converge"
    end if

  end subroutine

  subroutine ComputeCholesky(a, chol_a, compute_lower_triangular, info)
    ! Compute the Cholesky factorization of a positive definite matrix.
    !
    ! If `compute_lower_triangular = True` then computes a lower-triangular matrix
    ! `chol_a` such that
    !
    !   a = matmul(chol_a, transpose(chol_a))
    !
    ! otherwise computes an upper-triangular matrix `chol_a` such that
    !
    !   a = matmul(transpose(chol_a), chol_a)
    !
    ! Uses the LAPACK routine dpotrf.
    !
    ! Arguments:
    !
    ! a: Array containing positive definite matrix to be factorized, of size `(n, n)`.
    ! u: Array which on exit will contain a size `(n, n)` triangular matrix
    !   corresponding to the lower- or upper-triangular Choleksy factor of `a`.
    ! compute_lower_triangular: Whether to compute the lower-triangular Cholesky factor
    !   (True) or the upper-triangular Cholesky factor (False).
    ! info: Optional. If present, any non-zero output from the `info` diagnostic
    !  argument (corresponding to failure conditions) to the LAPACK routine dpotrf is
    !  outputted here to allow appropriate handling in downstream code. Otherwise if
    !  not present, execution is halted and an error message printed on failure
    !  conditions.

    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: chol_a(:, :)
    logical, intent(in) :: compute_lower_triangular
    integer(i4b), intent(out), optional :: info

    integer(i4b) :: n, local_info

    n = size(a, 1)
    if (size(a, 2) /= n) then
      stop "Matrix to be factorized `a` must be square"
    end if
    if (size(chol_a, 1) /= n .or. size(chol_a, 2) /= n) then
      stop "`chol_a` argument must be same size as `a`"
    end if

    ! dpotrf computes Cholesky factor in place on `a` argument therefore copy relevant
    ! triangle of `a` to `chol_a`
    if (compute_lower_triangular) then
      call CopyTriangular(a, chol_a, .True.)
      call dpotrf("l", n, chol_a, n, local_info)
    else
      call CopyTriangular(a, chol_a, .False.)
      call dpotrf("u", n, chol_a, n, local_info)
    end if

    if (local_info < 0) then
      stop "An argument to dpotrf has an illegal value"
    else if (present(info)) then
      info = local_info
    else if (local_info > 0) then
      stop "Matrix to be factored is not positive definite"
    end if

  end subroutine

  subroutine ComputeEigenSymmetric(a, w, v, info)
    ! Compute the eigenvalues or full eigendecomposition of a symmetric matrix.
    !
    ! For a symmetric matrix `a` of size `(n, n)` computes the length `n` array `w` such
    ! that `w` contains the eigenvalues of `a` in ascending order. Optionally, an array
    ! `v` is additionally computed such that the columns of `v` contains the orthonormal
    ! eigenvectors of `a` in the corresponding order to the eigenvalues in `w`.
    !
    ! For the full eigendecomposition, the computed factors are such that on exit
    !
    !   a = matmul(matmul(v, w_full), transpose(v))`
    !
    ! where `w_full` is a `(n, n)` diagonal matrix with `w` along its diagonal.
    !
    ! Uses the LAPACK subroutine dsyev.
    !
    ! Arguments:
    !
    ! a: Array containing symmetric matrix to be decomposed, of size `(n, n)`.
    ! w: Array which on exit will contain the `n` eigenvalues of `a` in ascending order.
    ! v: Optional. Array which if present on exit will contain a size `(n, n)`
    !   orthogonal matrix with columns corresponding to the eigenvectors of `a` in
    !   equivalent order to the eigenvalues in `w`.
    ! info: Optional. If present, any non-zero output from the `info` diagnostic
    !  argument (corresponding to failure conditions) to the LAPACK routine dsyev is
    !  outputted here to allow appropriate handling in downstream code. Otherwise if
    !  not present, execution is halted and an error message printed on failure
    !  conditions.

    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in) :: w(:)
    real(dp), intent(out), optional :: v(:, :)
    integer(i4b), intent(out), optional :: info

    integer(i4b) :: n, lwork, local_info
    real(dp), allocatable :: work(:), a_copy(:, :)

    n = size(a, 1)
    if (size(a, 2) /= n) then
      stop "Matrix to be decomposed `a` must be square."
    end if
    if (size(w, 1) /= n) then
        stop "`w` must be of size `n` where `n = size(a, 1)`"
    end if
    if (present(v)) then
      if (size(v, 1) /= n .or. size(v, 2) /= n) then
        stop "`v` must be of same size as `a`"
      end if
    end if

    if (present(v)) then
      ! `v` present therefore compute eigenvectors as well as eigenvalues
      ! dsyev computes eigenvectors in place on `a` argument therefore copy `a` to `v`
      v = a
      ! dsyev workspace query to get optimal lwork
      allocate(work(1))
      call dsyev("V", "U", n, v, n, w, work, -1, local_info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
      ! compute eigendecomposition of `a`
      call dsyev("V", "U", n, v, n, w, work, lwork, local_info)
    else
      ! `v` not present therefore compute only eigenvalues
      ! dsyev destroys `a` argument therefore copy `a`
      allocate(a_copy(n, n))
      a_copy = a
      ! dsyev workspace query to get optimal lwork
      allocate(work(1))
      call dsyev("N", "U", n, a_copy, n, w, work, -1, local_info)
      lwork = int(work(1))
      deallocate(work)
      allocate(work(lwork))
      ! compute eigendecomposition of `a`
      call dsyev("N", "U", n, a_copy, n, w, work, lwork, local_info)
    end if

    if (local_info < 0) then
      stop "An argument to dsyev has an illegal value"
    else if (present(info)) then
      info = local_info
    else if (local_info > 0) then
      stop "Iteration in dsyev did not converge"
    end if

  end subroutine

  subroutine InvertTriangular(a, inv_a, a_is_lower, info)
    ! Compute inverse of a triangular matrix.
    !
    ! Uses the LAPACK subroutine dtrtri.
    !
    ! Arguments:
    !
    ! a: Array containing triangular matrix to invert of size `(n, n)`.
    ! inv_a: Array which on output will contain inverse of `a`.
    ! a_is_lower: Whether `a` is lower triangular (True) or upper triangular (False).
    ! info: Optional. If present, any non-zero output from the `info` diagnostic
    !  argument (corresponding to failure conditions) to the LAPACK routine dtrtri is
    !  outputted here to allow appropriate handling in downstream code. Otherwise if
    !  not present, execution is halted and an error message printed on failure
    !  conditions.

    real(dp), intent(in) :: a(:, :)
    logical, intent(in) :: a_is_lower
    real(dp), intent(out) :: inv_a(:, :)
    integer(i4b), intent(out), optional :: info

    integer(i4b) :: n, local_info

    n = size(a, 1)
    if (size(a, 2) /= n) then
      stop "Matrix to be inverted must be square"
    end if

    ! dtrtri computes inverse in place therefore copy relevant triangle of a to inv_a
    call CopyTriangular(a, inv_a, a_is_lower)

    ! call dtrtri with `uplo` set based on `a_is_lower`
    if (a_is_lower) then
      call dtrtri("l", "n", n, inv_a, n, local_info)
    else
      call dtrtri("u", "n", n, inv_a, n, local_info)
    end if

    if (local_info < 0) then
      stop "Argument to dtrtri has invalid value"
    else if (present(info)) then
      info = local_info
    else if (local_info > 0) then
        stop "Triangular matrix to be inverted is singular with zero on diagonal"
    end if

  end subroutine

  subroutine InvertSymmetric(a, inv_a, info)
    ! Compute inverse of a symmetric matrix.
    !
    ! Uses the LAPACK subroutines dsytrf and dsytri.
    !
    ! Arguments:
    !
    ! a: Array containing symmetric matrix to invert of size `(n, n)`.
    ! inv_a: Array which on output will contain inverse of `a` of size `(n, n)`.
    ! info: Optional. If present, any non-zero output from the `info` diagnostic
    !  argument (corresponding to failure conditions) to the LAPACK routines dsytrf and
    !  dsytri is outputted here to allow appropriate handling in downstream code.
    !  Otherwise if not present, execution is halted and an error message printed on
    !  failure conditions.

    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: inv_a(:, :)
    integer(i4b), intent(out), optional :: info

    real(dp), allocatable :: work(:)
    integer(i4b), allocatable :: ipiv(:)
    integer(i4b) :: i, j, n, lwork, local_info

    n = size(a, 1)
    if (size(a, 2) /= n) then
      stop "Matrix to be inverted `a` must be square"
    end if
    if (size(inv_a, 1) /= n .or. size(inv_a, 2) /= n) then
      stop "Array to output matrix inverse in `inv_a` must match size of `a`."
    end if

    ! dsytri and dsytrf operate on `a` argument in-place therefore copy `a` to `inv_a`
    inv_a = a

    ! allocate storate for pivot array
    allocate(ipiv(n))

    ! Workspace query to get optimal lwork
    allocate(work(1))
    call dsytrf("u", n, inv_a, n, ipiv, work, -1, local_info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))

    ! Compute factorization of `a` using Bunch-Kaufman diagonal pivoting method.
    ! We use `uplo="u"` here which means only the upper triangle of `a` is referenced.
    ! Assuming `a` is actually symmetric we can arbitrarily use either the upper or
    ! lower triangle.
    call dsytrf("u", n, inv_a, n, ipiv, work, lwork, local_info)

    if (local_info < 0) then
      stop "Argument to dsytrf has invalid value"
    else if (present(info)) then
      info = local_info
      return
    else if (local_info > 0) then
      stop "Symmetric matrix to be inverted is singular"
    end if

    ! dsytri requires work array of size n
    deallocate(work)
    allocate(work(n))

    ! Compute inverse of `a` from factorization output from dsytrf
    ! On exit only upper triangular elements defined so need to symmetrize
    call dsytri("u", n, inv_a, n, ipiv, work, local_info)

    ! Copy elements of `inv_a` above diagonal to tranposed positions below diagonal
    do i = 2, n
      do j = 1, i - 1
        inv_a(i, j) = inv_a(j, i)
      end do
    end do

    if (local_info < 0) then
      stop "Argument to dsytri has invalid value"
    else if (present(info)) then
      info = local_info
    else if (local_info > 0) then
      stop "Symmetric matrix to be inverted is singular"
    end if

  end subroutine

  subroutine CheckSolveMatricesConformable(a, b, x)
    ! Check dimensions of arrays for linear system `matmul(a, x) = b` are conformable.
    real(dp), intent(in) :: a(:, :), b(:, :), x(:, :)
    if (size(a, 2) /= size(a, 1)) then
      stop "Matrix `a` must be square"
    end if
    if (size(b, 1) /= size(a, 2)) then
      stop "Right-hand-side matrix `b` must be of conformable size to matrix `a`"
    end if
    if (size(x, 1) /= size(b, 1) .or. size(x, 2) /= size(b, 2)) then
      stop "Solution matrix `x` must be of conformable size to matrices `a` and `b`"
    end if
  end subroutine

  subroutine SolveMultiple(a, b, x, info)
    ! Solve a linear system with a matrix right-hand-side.
    !
    ! Solves the linear systems
    !
    !    matmul(a, x) = b
    !
    ! for a general square matrix `a` and matrix `b`.
    !
    ! Uses the LAPACK subroutine dgsev.
    !
    ! Arguments:
    !
    ! a: Array containing triangular matrix of size `(m, m)`.
    ! b: Array containing matrix right-hand-side of linear system, of size `(m, n)`.
    ! x: Array which on exit will contain solution of linear system, of size `(m, n)`.
    ! info: Optional. If present, any non-zero output from the `info` diagnostic
    !  argument (corresponding to failure conditions) to the LAPACK routine dgsev is
    !  outputted here to allow appropriate handling in downstream code. Otherwise if
    !  not present, execution is halted and an error message printed on failure
    !  conditions.

    real(dp), intent(in) :: a(:, :), b(:, :)
    real(dp), intent(out) :: x(:, :)
    integer(i4b), intent(out), optional :: info

    real(dp), allocatable :: lu(:, :)
    integer(i4b), allocatable :: ipiv(:)
    integer(i4b) :: m, n, local_info

    m = size(a, 1)
    n = size(b, 2)

    ! check dimensions of arrays are conformable
    call CheckSolveMatricesConformable(a, b, x)

    ! `a` argument to dgesv overwritten with LU factors therefore assign `a` to `lu`
    allocate(lu(m, m))
    lu = a

    ! `b` argument to dgesv overwritten with solution `x` therefore assign `b` to `x`
    x = b

    ! `ipiv` pivots argument to dgesv should be of size `m`
    allocate(ipiv(m))

    ! Compute solution to linear system by performing pivoted LU decomposition of `a`
    call dgesv(m, n, lu, m, ipiv, x, m, local_info)

    if (local_info < 0) then
      stop "Argument to dgesv has invalid value"
    else if (present(info)) then
      info = local_info
    else if (local_info > 0) then
      stop "Matrix to be inverted in linear system is singular"
    end if

  end subroutine

  subroutine SolveSingle(a, b, x, info)
    ! Solve a linear system with a single vector right-hand-side.
    !
    ! Solves the linear systems
    !
    ! for a general square matrix `a` and vector `b`.
    !
    ! Uses the LAPACK subroutine dgsev.
    !
    ! Arguments:
    !
    ! a: Array containing matrix of size `(n, n)`.
    ! b: Array containing vector right-hand-side of linear system, of size `(n,)`.
    ! x: Array which on exit will contain solution of linear system, of size `(n,)`.
    ! info: Optional. If present, any non-zero output from the `info` diagnostic
    !  argument (corresponding to failure conditions) to the LAPACK routine dgsev is
    !  outputted here to allow appropriate handling in downstream code. Otherwise if
    !  not present, execution is halted and an error message printed on failure
    !  conditions.

    real(dp), intent(in) :: a(:, :)
    real(dp), target, intent(in) :: b(:)
    real(dp), target, intent(out) :: x(:)
    integer(i4b), intent(out), optional :: info

    real (dp), pointer :: b_2d(:, :), x_2d(:, :)
    integer(i4b) :: n

    ! Create rank-2 pointers to b and x to allow passing to SolveMultiple
    n = size(a, 1)
    b_2d(1:n, 1:1) => b
    x_2d(1:n, 1:1) => x
    if (present(info)) then
      call SolveMultiple(a, b_2d, x_2d, info)
    else
      call SolveMultiple(a, b_2d, x_2d)
    end if

  end subroutine

  subroutine SolveTriangularMultiple(a, b, x, a_is_lower, solve_transposed)
    ! Solve a triangular linear system with a matrix right-hand-side.
    !
    ! Solves one of the linear systems
    !
    !    matmul(a, x) = b
    !
    ! or
    !
    !    matmul(transpose(a), x) = b
    !
    ! for a upper or lower triangular matrix `a` and matrix `b`.
    !
    ! Uses the LAPACK subroutine dtrsm.
    !
    ! Arguments:
    !
    ! a: Array containing triangular matrix of size `(m, m)`.
    ! b: Array containing matrix right-hand-side of linear system, of size `(m, n)`.
    ! x: Array which on exit will contain solution of linear system, of size `(m, n)`.
    ! a_is_lower: Whether `a` is lower triangular (True) or upper triangular (False).
    ! solve_transposed: Whether to solve `matmul(a, x) = b` (True) or the transposed
    !   system `matmul(transpose(a), x) = b` (False).

    real(dp), intent(in) :: a(:, :), b(:, :)
    real(dp), intent(out) :: x(:, :)
    logical, intent(in) :: a_is_lower, solve_transposed

    integer(i4b) :: m, n
    m = size(a, 1)
    n = size(b, 2)

    ! check dimensions of arrays are conformable
    call CheckSolveMatricesConformable(a, b, x)

    ! dtrsm computes solution in place in `x` argument therefore copy RHS `b` to `x`
    x = b

    ! call dtrsm with `uplo` set based on `a_is_lower` and `trans` set based on
    ! `solve_transposed`
    if (a_is_lower .and. solve_transposed) then
      call dtrsm("l", "l", "t", "n", m, n, 1.0_dp, a, m, x, m)
    else if (a_is_lower) then
      call dtrsm("l", "l", "n", "n", m, n, 1.0_dp, a, m, x, m)
    else if ((.not. a_is_lower) .and. solve_transposed) then
      call dtrsm("l", "u", "t", "n", m, n, 1.0_dp, a, m, x, m)
    else
      call dtrsm("l", "u", "n", "n", m, n, 1.0_dp, a, m, x, m)
    end if

  end subroutine

  subroutine SolveTriangularSingle(a, b, x, a_is_lower, solve_transposed)
    ! Solve a triangular linear system with a single vector right-hand-side.
    !
    ! Solves one of the linear systems
    !
    !    matmul(a, x) = b
    !
    ! or
    !
    !    matmul(transpose(a), x) = b
    !
    ! for a upper or lower triangular matrix `a` and vector `b`.
    !
    ! Uses the LAPACK subroutine dtrsm.
    !
    ! Arguments:
    !
    ! a: Array containing triangular matrix of size `(n, n)`.
    ! b: Array containing vector right-hand-side of linear system, of size `(n,)`.
    ! x: Array which on exit will contain solution of linear system, of size `(n,)`.
    ! a_is_lower: Whether `a` is lower triangular (True) or upper triangular (False).
    ! solve_transposed: Whether to solve `matmul(a, x) = b` (True) or the transposed
    !   system `matmul(transpose(a), x) = b` (False).

    real(dp), intent(in) :: a(:, :)
    real(dp), target, intent(in) :: b(:)
    real(dp), target, intent(out) :: x(:)
    logical, intent(in) :: a_is_lower, solve_transposed

    real (dp), pointer :: b_2d(:, :), x_2d(:, :)
    integer(i4b) :: n

    ! Create rank-2 pointers to b and x to allow passing to SolveTriangularMultiple
    n = size(a, 1)
    b_2d(1:n, 1:1) => b
    x_2d(1:n, 1:1) => x
    call SolveTriangularMultiple(a, b_2d, x_2d, a_is_lower, solve_transposed)

  end subroutine

  function LogAbsDet(a) result(log_abs_det)
    ! Compute the natural logarithm of the absolute value of the determinant of a matrix
    !
    ! Returns negative infinity if `a` matrix is singular.
    !
    ! Uses LAPACK routine dgetrf.
    !
    ! Arguments:
    !
    ! a: Array containing square matrix to compute log absolute determinant of.

    real(dp), intent(in) :: a(:, :)
    real(dp) :: log_abs_det

    integer(i4b) :: i, n, info
    real(dp), allocatable :: lu(:, :)
    integer(i4b), allocatable :: ipiv(:)

    n = size(a, 1)
    if (size(a, 2) /= n) then
      stop "Matrix `a` to compute log determinant of must be square"
    end if

    allocate(lu(n, n), ipiv(n))

    ! `a` argument overwritten by `lu` factors in dgetrf therefore copy `a` to `lu`
    lu = a

    ! Compute partially pivoted LU factorization of `a`
    ! On exit, `lu` will contain upper-triangular matrix `u` in elements above and on
    ! diagonal and off-diagonal elements of unit diagonal lower-triangular matrix `l`
    ! below the diagonal. The array `ipiv` will contain the pivot indices.
    call dgetrf(n, n, lu, n, ipiv, info)

    if (info < 0) then
      stop "Argument to dgetrf has invalid value"
    else if (info > 0) then
      ! matrix is singular therefore log determinant is negative infinity
      log_abs_det = ieee_value(log_abs_det, ieee_negative_inf)
    else
      ! log absolute determinant is equal to sum logarithms of absolute values of
      ! diagonal elements of `u` (as `l` is unit diagonal and so zero log determinant
      ! and the permutation matrix represented by `ipiv` has determinant +/- 1 so zero
      ! log absolute determinant)
      log_abs_det = 0.0_dp
      do i = 1, n
        log_abs_det = log_abs_det + log(abs(lu(i, i)))
      end do
    end if

  end function

end module