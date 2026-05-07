!===============================================================
!  Test Dense Direct Solver (MKL + LAPACK dgesv)
!  Automatic PASS / FAIL check
!===============================================================
program test_denseSolver
  use LinearSolverBase_mod, only: dp
  use MatrixDense_mod
  use MKLDirectSolver_mod
  implicit none

  type(DenseMatrix)       :: A
  type(MKLDirectSolver)   :: solver
  real(dp), allocatable   :: b(:), x(:), x_exact(:)

  integer :: n, i
  real(dp) :: err, norm_exact, tol

  print *, "=============================================="
  print *, " Test Dense Direct Solver (LAPACK dgesv)"
  print *, "=============================================="

  tol = 1.0d-10

  !----------------------------
  ! Problem size
  !----------------------------
  n = 3
  call A%init(n, .false.)
  allocate(b(n), x(n), x_exact(n))

  !----------------------------
  ! Fill dense matrix A
  !----------------------------
  A%A = reshape( [ &
     3.0_dp,  2.0_dp, -1.0_dp, &
     2.0_dp, -2.0_dp,  4.0_dp, &
    -1.0_dp,  0.5_dp, -1.0_dp  &
  ], shape(A%A), order=[2,1] )

  ! RHS
  b = [ 1.0_dp, -2.0_dp, 0.0_dp ]

  ! Exact solution (from MATLAB / analytic)
  x_exact = [ 1.0_dp, -2.0_dp, -2.0_dp ]

  print *, "Dense matrix A:"
  call A%print()

  !----------------------------
  ! Solve
  !----------------------------
  call solver%attach_matrix(A)
  call solver%analyze()
  call solver%factor()
  call solver%solve(b, x)

  !----------------------------
  ! Compute relative error
  !----------------------------
  err = 0.0_dp
  norm_exact = 0.0_dp
  do i = 1, n
    err        = err        + (x(i) - x_exact(i))**2
    norm_exact = norm_exact +  x_exact(i)**2
  end do

  err        = sqrt(err)
  norm_exact = sqrt(norm_exact)
  err        = err / norm_exact

  !----------------------------
  ! Report
  !----------------------------
  print *, "Computed solution x:"
  do i = 1, n
    print '(A,I1,A,F12.6)', "  x(",i,") = ", x(i)
  end do

  print *, "Relative error ||x-x_exact|| / ||x_exact|| =", err

  if (err < tol) then
    print *, ">>> TEST PASSED (tol =", tol, ")"
  else
    print *, ">>> TEST FAILED (tol =", tol, ")"
    stop 1
  end if

  call solver%free()

end program test_denseSolver
