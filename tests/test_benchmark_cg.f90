program test_benchmark_cg
  use Types, only : prec
  use MeshIO_mod, only : load_mesh_from_txt
  use Assembler_OOP_mod, only : assemble_global_stiffness
  use BoundaryCondition_mod, only : apply_benchmark_loads, apply_benchmark_bcs_csr 
  use SparseCOO_mod, only : SparseCOO
  use SparseCSR_mod, only : SparseCSR
  use SparseConvert_mod, only : COO_to_CSR
  use LinearSolverBase_mod, only : LinearSolver, LinearSolverOptions
  use LinearSolverFactory_mod, only : create_linear_solver
  implicit none

  type(SparseCOO) :: Kcoo
  type(SparseCSR), target :: Kcsr
  class(LinearSolver), allocatable :: solver
  type(LinearSolverOptions) :: opts
  
  real(prec), allocatable :: rhs(:), x(:), coords(:,:)
  integer, allocatable :: conn(:,:)
  integer :: n_nodes, n_elements, n_unknowns
  real(prec) :: props(2)

  real(8) :: t_start, t_end, t_assembly, t_convert, t_solve

  print *, "--- [OOP FEM Engine: Standard CG (No Precon) Test] ---"
  
  call load_mesh_from_txt('/work/tests/mesh_scale.txt', n_nodes, n_elements, coords, conn, props)
  n_unknowns = n_nodes * 3
  
  ! --- 1. 装配阶段 ---
  call cpu_time(t_start)
  call Kcoo%init(n_unknowns)
  call assemble_global_stiffness(Kcoo, n_elements, n_nodes, coords, conn, props)
  call cpu_time(t_end)
  t_assembly = t_end - t_start

  ! --- 2. 转换阶段 ---
  call cpu_time(t_start)
  call Kcoo%preprocess() 
  Kcsr = COO_to_CSR(Kcoo)
  call cpu_time(t_end)
  t_convert = t_end - t_start

  ! --- 3. 边界与载荷施加 (CSR 置零置一法) ---
  allocate(rhs(n_unknowns), x(n_unknowns))
  call apply_benchmark_loads(n_nodes, coords, rhs)
  call apply_benchmark_bcs_csr(n_nodes, coords, Kcsr, rhs)

  ! --- 4. 标准 CG 求解阶段 ---
  call cpu_time(t_start)
  
  opts%solver_family = "cg"      ! <--- 唯一指定：无预处理的标准 CG
  ! opts%max_iter = 50000        ! 如果它很难收敛，可以放宽最大迭代次数
  
  call create_linear_solver(Kcsr, opts, solver)
  call solver%analyze(); call solver%factor(); call solver%solve(rhs, x)
  
  call cpu_time(t_end)
  t_solve = t_end - t_start

  ! 检验解的正确性
  print *, " 顶部受拉点位移测试值: ", x(n_unknowns)

  ! 🏆 打印压测报告
  print *, "========================================="
  print *, "     标准 CG 法 (无预处理) 性能报告      "
  print *, "========================================="
  print *, " 单元数量 : ", n_elements
  print *, " 总自由度 : ", n_unknowns
  print *, " 非零元素 : ", Kcsr%nnz   
  print *, " 装配耗时 : ", t_assembly, " 秒"
  print *, " 转换耗时 : ", t_convert,  " 秒"
  print *, " 求解耗时 : ", t_solve,    " 秒"
  print *, "========================================="

  call solver%free()
end program test_benchmark_cg