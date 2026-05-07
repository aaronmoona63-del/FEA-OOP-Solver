program test_skyline_lu
  use Types
  use MatrixSkyline
  use SkylineLUSolverMod
  implicit none

  type(SkylineMatrix)     :: K
  type(SkylineLUSolver)   :: solver
  integer :: n

  real(prec) :: rhs(5), x(5)

  !----------------------------------------------------------
  ! 构造一个 5×5 对称 skyline 系统
  ! 对示例矩阵：
  ! |K11 K12 K13  0   0 |
  ! |K21 K22 K23 K24  0 |
  ! |K31 K32 K33 K34 K35|
  ! | 0  K42 K43 K44 K45|
  ! | 0   0  K53 K54 K55|
  !
  ! AUPP = [K12, K13, K23, K24, K34, K35, K45]   (长度 7)
  ! JPOIN = [0, 1, 3, 5, 7]                     (长度 5)
  !----------------------------------------------------------

  n = 5

  ! 手工设置 skyline 结构（不调用 init）
  K%nrows = n
  K%ncols = n
  K%is_symmetric = .false.
  K%unsymmetric  = .true.

  !-------------------------------------------
  ! allocate skyline arrays
  !-------------------------------------------
  allocate(K%jpoin(0:n))
  allocate(K%diag(n))
  allocate(K%aupp(7))
  allocate(K%alow(7))

  ! K%jpoin = [0, 1, 3, 5, 7]
  ! K%diag = [4.0, 5.0, 6.0, 7.0, 8.0]
  ! K%aupp =[1.0000,2.0000,1.5000,0.5000,3.0000,1.0000,0.2000] 
  ! K%alow = K%aupp
  ! rhs = [1, 1, 1, 1, 1]

  
  K%aupp =  [0.2000, 0.4000, 0.3000, 0.1000,  0.,  0.,  0.5000, 0.7000]
  K%alow = K%aupp
  K%diag = [1.0, 2.0, 3.0, 4.0, 5.0]
  K%jpoin = [0, 1, 3, 6, 8]
  rhs = [1, 1, 1, 1, 1]

  call K%print()


  !-------------------------------------------
  ! Factorize + Solve
  !-------------------------------------------
  call solver%factorize(K)

  print*, 'diag'
  print '(5f12.5)', K%diag
  print*, 'alow'
  print '(5f12.5)', K%alow
  print*, 'aupp'
  print '(5f12.5)', K%aupp
  print*, 'jpoin'
  print '(5i6)',    K%jpoin


  call solver%solve(K, rhs, x)

  print *, "Computed solution x = "
  print *, x


end program test_skyline_lu
