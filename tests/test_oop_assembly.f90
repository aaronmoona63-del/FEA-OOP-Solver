program test_oop_assembly
  use Types, only : prec
  use MeshIO_mod, only : load_mesh_from_txt
  use Assembler_OOP_mod, only : assemble_global_stiffness
  use ResultIO_mod, only : export_vtk_hex8
  ! 👉 注意这里：我们引入了全新且绝对安全的独立后处理模块
  use PostProcessor_mod, only : calc_hex8_von_mises
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
  integer :: i, nid, n_nodes, n_elements, n_unknowns
  real(prec) :: props(2), big_number

  ! 用于应力计算的临时变量
  real(prec), allocatable :: von_mises_list(:)
  real(prec) :: coords_elem(3, 8), u_elem(24), vm
  integer :: el, j

  print *, "--- [OOP FEM Engine: Full Pipeline Test] ---"

  call load_mesh_from_txt('/work/tests/mesh_3d.txt', n_nodes, n_elements, coords, conn, props)
  n_unknowns = n_nodes * 3
  props(1) = 210000.0_prec  
  props(2) = 0.3_prec       

  call Kcoo%init(n_unknowns)
  call assemble_global_stiffness(Kcoo, n_elements, n_nodes, coords, conn, props)

  allocate(rhs(n_unknowns), x(n_unknowns))
  rhs = 0.0_prec; x = 0.0_prec
  rhs(3*11 - 2) = 25000000.0_prec
  rhs(3*9  - 2) = 25000000.0_prec
  rhs(3*10 - 2) = 25000000.0_prec
  rhs(3*12 - 2) = 25000000.0_prec

  big_number = 1.0e15_prec
  do i = 1, 4
      if(i==1) nid = 3; if(i==2) nid = 4
      if(i==3) nid = 6; if(i==4) nid = 7
      call Kcoo%add_entry(3*nid-2, 3*nid-2, big_number)
      call Kcoo%add_entry(3*nid-1, 3*nid-1, big_number)
      call Kcoo%add_entry(3*nid,   3*nid,   big_number)
  end do

  call Kcoo%preprocess() 
  Kcsr = COO_to_CSR(Kcoo)
  opts%solver_family = "direct"
  call create_linear_solver(Kcsr, opts, solver)
  call solver%analyze(); call solver%factor(); call solver%solve(rhs, x)

  print *, "--- Global Solve Completed ---"

  ! =========================================================
  ! 👉 安全应力恢复 (不再依赖 Element_Hex8，直接调用独立计算器)
  ! =========================================================
  allocate(von_mises_list(n_elements))
  do el = 1, n_elements
      do j = 1, 8
          coords_elem(:, j) = coords(:, conn(j, el))
          u_elem(3*j-2)     = x(3*conn(j, el) - 2)
          u_elem(3*j-1)     = x(3*conn(j, el) - 1)
          u_elem(3*j)       = x(3*conn(j, el))
      end do
      ! 直接调用独立函数算应力，完全解耦！
      call calc_hex8_von_mises(coords_elem, props, u_elem, vm)
      von_mises_list(el) = vm
  end do

  print *, "================================"
  print *, "   OOP Result (Node 9)          "
  print *, "================================"
  print *, "Ux: ", x(3*9 - 2)
  print *, "Uy: ", x(3*9 - 1)
  print *, "Uz: ", x(3*9)

  call export_vtk_hex8('/work/tests/oop_result_3d.vtk', n_nodes, n_elements, coords, conn, x, von_mises_list)

  call solver%free()
end program test_oop_assembly