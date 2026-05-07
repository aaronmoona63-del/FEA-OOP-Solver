program test_csr_pardiso
  use Types
  use MatrixCSR
  use SolverBase
  use PardisoSolverMod
  implicit none

  type(CSRMatrix) :: A
  class(DirectSolverType), allocatable :: solver
  integer :: n
  real(prec), allocatable :: rhs(:), x(:), r(:)

  n = 3
  allocate(rhs(n), x(n), r(n))

  !---------------------------------------------
  ! A = [ 4 1 2
  !       1 3 0
  !       2 0 5 ]
  !
  ! CSR (1-based):
  ! row1: (1,1)=4, (1,2)=1, (1,3)=2
  ! row2: (2,1)=1, (2,2)=3
  ! row3: (3,1)=2, (3,3)=5
  !
  ! row_ptr = [1, 4, 6, 8]
  ! col_ind = [1,2,3,  1,2,  1,3]
  ! values  = [4,1,2,  1,3,  2,5]
  !---------------------------------------------

  A%n = n
  allocate(A%row_ptr(n+1))
  allocate(A%col_ind(7))
  allocate(A%values(7))

  A%row_ptr = (/ 1, 4, 6, 8 /)
  A%col_ind = (/ 1, 2, 3,  1, 2,  1, 3 /)
  A%values  = (/ 4.0_prec, 1.0_prec, 2.0_prec,  &
                 1.0_prec, 3.0_prec,            &
                 2.0_prec, 5.0_prec /)

  rhs = (/ 7.0_prec, 4.0_prec, 7.0_prec /)

  !---------------------------------------------
  ! 创建 PARDISO 求解器 & 求解
  !---------------------------------------------
  allocate(PardisoSolver :: solver)

  call solver%factorize(A)
  call solver%solve(A, rhs, x)

  ! 计算残差 r = A x - rhs
  call A%matvec(x, r)
  r = r - rhs

  print *, "Solution x ="
  print '(3F12.6)', x

  print *, "Residual r = A*x - b ="
  print '(3E16.6)', r

end program test_csr_pardiso
