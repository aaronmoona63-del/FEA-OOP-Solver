program test_tension
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

  print *, '--- [Test 1: Pure Tension (单轴拉伸)] ---'
  n_elements = 2
  n_unknowns = 12
  
  call assemble_all_elements_coo(n_elements, n_unknowns, Kcoo, rhs, fail)
  if (fail) stop 'Assembly failed!'

  ! 边界条件：固定左端 Node 1 (Uy, Ux) 和 Node 4 (Ux)
  big_number = 1.0e12_prec
  call Kcoo%add_entry(1, 1, big_number)
  call Kcoo%add_entry(2, 2, big_number)
  call Kcoo%add_entry(7, 7, big_number)
  
  ! 施加载荷：右侧节点向右拉 (Fx = +1.0)
  rhs(5)  = rhs(5) + 1.0_prec
  rhs(11) = rhs(11) + 1.0_prec

  Kcsr = COO_to_CSR(Kcoo)
  allocate(x(n_unknowns)); x = 0.0_prec
  opts%solver_family = "direct"
  call create_linear_solver(Kcsr, opts, solver)
  
  ! -----------------------------------------------
  ! MKL 直接求解器标准三步曲
  ! -----------------------------------------------
  call solver%analyze()
  call solver%factor()
  call solver%solve(rhs, x)

  print *, 'Tip Displacement Ux (Node 3):', x(5)
  print *, 'Tip Displacement Ux (Node 6):', x(11)
  
  ! 验证材料受拉伸时上下变形均匀
  if (abs(x(5) - x(11)) < 1e-5_prec) then
      print *, 'Tension Test PASSED (Uniform deformation observed!)'
  end if

  ! VTK 导出
  allocate(coords(2, 6), conn(4, 2))
  coords(:,1)=[0.,0.]; coords(:,2)=[1.,0.]; coords(:,3)=[2.,0.]
  coords(:,4)=[0.,1.]; coords(:,5)=[1.,1.]; coords(:,6)=[2.,1.]
  conn(:,1)=[1,2,5,4]; conn(:,2)=[2,3,6,5]
  call export_to_vtk('result_tension.vtk', 6, 2, coords, conn, x)
  
  call solver%free()
end program test_tension