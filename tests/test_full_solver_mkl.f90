program test_full_solver_mkl
  use Types, only : prec
  use SparseCOO_mod, only : SparseCOO
  use SparseCSR_mod, only : SparseCSR
  use SparseConvert_mod, only : COO_to_CSR
  use AssemblerCOO_mod, only : assemble_all_elements_coo
  use LinearSolverBase_mod, only : LinearSolver, LinearSolverOptions
  use LinearSolverFactory_mod, only : create_linear_solver
  use VTK_Export_mod, only : export_to_vtk
  implicit none

  type(SparseCOO) :: Kcoo
  type(SparseCSR), target :: Kcsr
  class(LinearSolver), allocatable :: solver
  type(LinearSolverOptions) :: opts
  real(prec), allocatable :: rhs(:), x(:)
  logical :: fail
  integer :: n_elements, n_unknowns, i
  real(prec) :: big_number
  real(prec), allocatable :: coords(:,:), conn_vtk(:,:)
  integer, allocatable :: conn(:,:)

  print *, '--- [MKL + VTK Pipeline Test (悬臂梁受力)] ---'
  n_elements = 2
  n_unknowns = 12
  
  call assemble_all_elements_coo(n_elements, n_unknowns, Kcoo, rhs, fail)
  if (fail) stop 'Assembly failed!'

  ! 边界条件：罚函数法固定左端 (Node 1 & Node 4)
  big_number = 1.0e12_prec
  call Kcoo%add_entry(1, 1, big_number)
  call Kcoo%add_entry(2, 2, big_number)
  call Kcoo%add_entry(7, 7, big_number)
  call Kcoo%add_entry(8, 8, big_number)
  
  ! 施加载荷：Node 3 向下受力 (Fy = -1.0) --> 恢复小变形假设
  rhs(6) = rhs(6) - 1.0_prec

  Kcsr = COO_to_CSR(Kcoo)
  allocate(x(n_unknowns)); x = 0.0_prec
  opts%solver_family = "direct"
  call create_linear_solver(Kcsr, opts, solver)
  
  ! -----------------------------------------------
  ! MKL 直接求解器标准三步曲（分析 -> 分解 -> 求解）
  ! -----------------------------------------------
  call solver%analyze()
  call solver%factor()
  call solver%solve(rhs, x)

  print *, 'Tip Displacement (Node 3, Uy):', x(6)

  ! VTK 导出
  allocate(coords(2, 6), conn(4, 2))
  coords(:,1)=[0.,0.]; coords(:,2)=[1.,0.]; coords(:,3)=[2.,0.]
  coords(:,4)=[0.,1.]; coords(:,5)=[1.,1.]; coords(:,6)=[2.,1.]
  conn(:,1)=[1,2,5,4]; conn(:,2)=[2,3,6,5]
  call export_to_vtk('beam_results.vtk', 6, 2, coords, conn, x)
  
  ! 验证结果：受力变小了，判定标准也对应缩小 (1e-6)
  if (x(6) < -1.0e-6_prec) then
      print *, 'test_full_solver_mkl PASSED'
  else
      print *, 'test_full_solver_mkl FAILED: displacement too small'
      stop 1
  end if
  
  call solver%free()
end program test_full_solver_mkl