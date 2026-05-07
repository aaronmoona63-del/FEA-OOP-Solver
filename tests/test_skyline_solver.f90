program test_skyline_solver
  use Types
  use MatrixSkyline
  use SolverBase
  use SkylineSolverMod
  implicit none

  integer :: n, i
  real(prec), allocatable :: rhs(:), x(:)
  class(DirectSolverType), allocatable :: solver   !<-- FIX 1
  type(SkylineMatrix) :: K                         !<-- FIX 2

  !---------------------------------------------
  ! 1. SPD 5×5 测试矩阵
  !---------------------------------------------
  n = 5
  allocate(rhs(n), x(n))
  rhs = 1.0_prec

  ! 初始化 skyline profile
  call K%init_from_profile(n, [1,1,2,3,4])

  ! 设置矩阵（Tridiagonal SPD）
  ! row 1
  K%a(K%prof_ptr(2)-1) = 4.0_prec
  ! row 2
  K%a(K%prof_ptr(2))   = -1.0_prec
  K%a(K%prof_ptr(3)-1) =  4.0_prec
  ! row 3
  K%a(K%prof_ptr(3))   = -1.0_prec
  K%a(K%prof_ptr(4)-1) =  4.0_prec
  ! row 4
  K%a(K%prof_ptr(4))   = -1.0_prec
  K%a(K%prof_ptr(5)-1) =  4.0_prec
  ! row 5
  K%a(K%prof_ptr(5))   = -1.0_prec
  K%a(K%prof_ptr(6)-1) =  4.0_prec

  !---------------------------------------------
  ! 2. 求解
  !---------------------------------------------
  allocate(SkylineSolver :: solver)

  call solver%factorize(K)
  call solver%solve(K, rhs, x)

  !---------------------------------------------
  ! 3. 输出
  !---------------------------------------------
  write(*,*) "Solution ="
  do i = 1, n
     write(*,'(F12.6)') x(i)
  end do

end program test_skyline_solver
