program test_benchmark
  use Types, only : prec
  use MeshIO_mod, only : load_mesh_from_txt
  use Assembler_OOP_mod, only : assemble_global_stiffness
  ! 👉 这里换成新的边界调用模块
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

  print *, "--- [OOP FEM Engine: DIRECT Solver Performance Test] ---"
  
  call load_mesh_from_txt('/work/tests/mesh_scale.txt', n_nodes, n_elements, coords, conn, props)
  n_unknowns = n_nodes * 3
  
  call cpu_time(t_start)
  call Kcoo%init(n_unknowns)
  call assemble_global_stiffness(Kcoo, n_elements, n_nodes, coords, conn, props)
  call cpu_time(t_end)
  t_assembly = t_end - t_start

  call cpu_time(t_start)
  call Kcoo%preprocess() 
  Kcsr = COO_to_CSR(Kcoo)
  call cpu_time(t_end)
  t_convert = t_end - t_start

  ! 👉 这里改成和迭代法一模一样的 CSR 边界处理
  allocate(rhs(n_unknowns), x(n_unknowns))
  call apply_benchmark_loads(n_nodes, coords, rhs)
  call apply_benchmark_bcs_csr(n_nodes, coords, Kcsr, rhs)

  call cpu_time(t_start)
  opts%solver_family = "direct"
  call create_linear_solver(Kcsr, opts, solver)
  call solver%analyze(); call solver%factor(); call solver%solve(rhs, x)
  call cpu_time(t_end)
  t_solve = t_end - t_start
  ! 检验解的正确性 (打印顶部中心点的一个位移)
  print *, " 顶部受拉点位移测试值: ", x(n_unknowns)

  print *, "========================================="
  print *, "       直接法 (PARDISO) 性能报告         "
  print *, "========================================="
  print *, " 单元数量 : ", n_elements
  print *, " 总自由度 : ", n_unknowns          ! <--- 新增这行
  print *, " 非零元素 : ", Kcsr%nnz            ! <--- 新增这行
  print *, " 装配耗时 : ", t_assembly, " 秒"
  print *, " 转换耗时 : ", t_convert,  " 秒"
  print *, " 求解耗时 : ", t_solve,    " 秒"
  print *, "========================================="
  call solver%free()
end program test_benchmark