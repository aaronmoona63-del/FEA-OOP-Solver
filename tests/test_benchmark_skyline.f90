program test_benchmark_skyline
  use Types, only : prec
  use MeshIO_mod, only : load_mesh_from_txt
  use Assembler_OOP_mod, only : assemble_global_stiffness
  use BoundaryCondition_mod, only : apply_benchmark_loads
  use SparseCOO_mod, only : SparseCOO
  use SparseSkyline_mod, only : SparseSkyline
  use SparseConvert_mod, only : COO_to_Skyline
  use LinearSolverBase_mod, only : LinearSolver, LinearSolverOptions
  use LinearSolverFactory_mod, only : create_linear_solver
  implicit none

  type(SparseCOO) :: Kcoo
  type(SparseSkyline), target :: Ksky
  class(LinearSolver), allocatable :: solver
  type(LinearSolverOptions) :: opts
  
  real(prec), allocatable :: rhs(:), x(:), coords(:,:)
  integer, allocatable :: conn(:,:)
  integer :: n_nodes, n_elements, n_unknowns
  real(prec) :: props(2)
  real(8) :: t_start, t_end, t_assembly, t_convert, t_solve
  
  integer :: i
  real(prec) :: big_number

  print *, "--- [OOP FEM Engine: Traditional Skyline Solver Test] ---"
  
  call load_mesh_from_txt('/work/tests/mesh_scale.txt', n_nodes, n_elements, coords, conn, props)
  n_unknowns = n_nodes * 3
  
  ! --- 1. 装配阶段 ---
  call cpu_time(t_start)
  call Kcoo%init(n_unknowns)
  call assemble_global_stiffness(Kcoo, n_elements, n_nodes, coords, conn, props)
  
  allocate(rhs(n_unknowns), x(n_unknowns))
  call apply_benchmark_loads(n_nodes, coords, rhs)
  
  ! 施加位移约束 (使用大数惩罚法，直接作用于 COO 矩阵)
  big_number = 1.0e15_prec
  do i = 1, n_nodes
      if (abs(coords(3, i) - 0.0_prec) < 1.0e-5_prec) then
          call Kcoo%add_entry(3*i-2, 3*i-2, big_number)
          call Kcoo%add_entry(3*i-1, 3*i-1, big_number)
          call Kcoo%add_entry(3*i,   3*i,   big_number)
      end if
  end do
  call cpu_time(t_end)
  t_assembly = t_end - t_start

  ! --- 2. 转换阶段 (COO -> Skyline) ---
  print *, "=> 正在转换矩阵格式为 Skyline (这可能很占内存)..."
  call cpu_time(t_start)
  call Kcoo%preprocess() 
  ! 👉 注意：必须传入 .true.，告诉程序只提取下三角(对称)
  Ksky = COO_to_Skyline(Kcoo, .true.)
  call cpu_time(t_end)
  t_convert = t_end - t_start

  ! --- 3. Skyline 求解阶段 ---
  print *, "=> 正在进行 Skyline 直接分解与求解..."
  call cpu_time(t_start)
  opts%solver_family = "skyline"
  call create_linear_solver(Ksky, opts, solver)
  call solver%analyze(); call solver%factor(); call solver%solve(rhs, x)
  call cpu_time(t_end)
  t_solve = t_end - t_start
  ! 检验解的正确性 (打印顶部中心点的一个位移)
  print *, " 顶部受拉点位移测试值: ", x(n_unknowns)

  ! 🏆 打印压测报告
  print *, "========================================="
  print *, "     传统 Skyline 直接法性能报告         "
  print *, "========================================="
  print *, " 单元数量 : ", n_elements
  print *, " 总自由度 : ", n_unknowns
  print *, " 剖面规模 : ", size(Ksky%al)  ! 👉 打印 Skyline 下三角数组的大小
  print *, " 装配耗时 : ", t_assembly, " 秒"
  print *, " 转换耗时 : ", t_convert,  " 秒"
  print *, " 求解耗时 : ", t_solve,    " 秒"
  print *, "========================================="

  call solver%free()
end program test_benchmark_skyline