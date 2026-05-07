program demo_skyline_csr
  use Types
  use MatrixSkyline
  use MatrixCSR
  use SolverBase
  use SkylineSolverMod
  implicit none

  type(SkylineMatrix) :: Ask
  type(CSRMatrix)     :: Acsr
  class(SolverType), allocatable :: solver

  integer :: n
  real(prec), allocatable :: x(:), b(:), r(:)

  n = 3

  !------------------------------
  ! 1. 构造 Skyline 矩阵 A
  !    A = [ 4 1 0
  !          1 3 1
  !          0 1 2 ]
  ! profile 设为满下三角: 每行长度 1,2,3 → 总长度 6
  !------------------------------
  Ask%n = n
  Ask%m = n
  allocate(Ask%prof_ptr(n+1))
  allocate(Ask%a(6))

  ! prof_ptr(i) 指向第 i 行在 a 中的起始位置
  ! 行长度: 1,2,3 → 起点: 1,2,4,7
  Ask%prof_ptr = (/ 1, 2, 4, 7 /)

  ! a 的布局：
  ! row1: (1,1) → 4
  ! row2: (2,1),(2,2) → 1,3
  ! row3: (3,1),(3,2),(3,3) → 0,1,2
  Ask%a = (/ 4.0_prec, &
             1.0_prec, 3.0_prec, &
             0.0_prec, 1.0_prec, 2.0_prec /)

  allocate(b(n), x(n), r(n))
  b = (/ 1.0_prec, 2.0_prec, 3.0_prec /)

  !------------------------------
  ! 2. 创建 SkylineSolver，并 factorize + solve
  !------------------------------
  allocate(SkylineSolver :: solver)

  call solver%factorize(Ask)
  call solver%solve(Ask, b, x)

  print *, 'Solution x = '
  print '(3F12.6)', x

  !------------------------------
  ! 3. 构造同一个 A 的 CSR，用 matvec 检查残差 r = A x - b
  !------------------------------
  call build_csr_example(Acsr)
  call Acsr%matvec(x, r)
  r = r - b

  print *, 'Residual r = A x - b = '
  print '(3E16.6)', r

  call solver%destroy()

contains

  ! 构造与 Skyline 相同的 3x3 矩阵 A 的 CSR 表示
  subroutine build_csr_example(A)
    type(CSRMatrix), intent(out) :: A

    A%n = 3
    A%m = 3
    allocate(A%row_ptr(4))
    allocate(A%col_ind(7))
    allocate(A%values(7))

    ! 按行：
    ! row1: (1,1)=4, (1,2)=1
    ! row2: (2,1)=1, (2,2)=3, (2,3)=1
    ! row3: (3,2)=1, (3,3)=2
    A%row_ptr = (/ 1, 3, 6, 8 /)
    A%col_ind = (/ 1,2,  1,2,3,  2,3 /)
    A%values  = (/ 4.0_prec, 1.0_prec, &
                   1.0_prec, 3.0_prec, 1.0_prec, &
                   1.0_prec, 2.0_prec /)
  end subroutine build_csr_example

end program demo_skyline_csr
