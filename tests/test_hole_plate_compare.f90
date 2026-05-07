program test_hole_plate_compare
  use Types, only : prec
  use MeshIO_mod, only : load_mesh_from_txt
  use Assembler_OOP_mod, only : assemble_global_stiffness
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
  integer :: i, s, hole_node_id, corner_node_id, count_right
  real(prec) :: props(2), min_dist_hole, min_dist_corner, dist
  real(prec) :: hole_apex_ux, corner_ux, l2_norm
  
  ! 用于计算右端面面积和力的变量
  real(prec) :: y_max, y_min, z_max, z_min, face_area, total_force, nodal_force

  real(8) :: t_start, t_end
  character(len=20) :: solver_names(3)

  print *, "========================================================="
  print *, " 🚀 带孔板拉伸：力控制(Neumann)物理极度对标比武 🚀"
  print *, "========================================================="

  ! --- 1. 加载网格 ---
  print *, ">>> [1/3] 正在读取带孔板网格数据..."
  call load_mesh_from_txt('/work/tests/mesh_daikonglashen.txt', n_nodes, n_elements, coords, conn, props)
  n_unknowns = n_nodes * 3
  
  ! --- 寻找两个核心监控点 ---
  hole_node_id = 1
  corner_node_id = 1
  min_dist_hole = 1.0e10_prec
  min_dist_corner = 1.0e10_prec
  
  do i = 1, n_nodes
      ! 测点 A：孔口上方顶点
      dist = sqrt((coords(1, i) - 0.0_prec)**2 + (coords(2, i) - 20.0_prec)**2 + (coords(3, i) - 5.0_prec)**2)
      if (dist < min_dist_hole) then
          min_dist_hole = dist
          hole_node_id = i
      end if
      ! 测点 B：右端面最上角
      dist = sqrt((coords(1, i) - 100.0_prec)**2 + (coords(2, i) - 50.0_prec)**2 + (coords(3, i) - 5.0_prec)**2)
      if (dist < min_dist_corner) then
          min_dist_corner = dist
          corner_node_id = i
      end if
  end do

  ! --- 2. 纯净装配与 Neumann 边界施加 ---
  print *, ">>> [2/3] 正在装配纯刚度矩阵..."
  call Kcoo%init(n_unknowns)
  call assemble_global_stiffness(Kcoo, n_elements, n_nodes, coords, conn, props)
  call Kcoo%preprocess() 
  
  allocate(rhs(n_unknowns), x(n_unknowns))
  rhs = 0.0_prec

  ! =====================================================================
  ! 🌟 核心修改：精准施加 100 MPa 拉伸面力 (Neumann 边界)
  ! =====================================================================
  print *, "      正在自动识别端面并施加 100 MPa 分布拉力..."
  y_max = -1.0e10_prec; y_min = 1.0e10_prec
  z_max = -1.0e10_prec; z_min = 1.0e10_prec
  count_right = 0
  
  ! 第一遍循环：寻找右端面尺寸以计算受力面积
  do i = 1, n_nodes
      if (abs(coords(1, i) - 100.0_prec) < 1.0e-3_prec) then
          count_right = count_right + 1
          if (coords(2, i) > y_max) y_max = coords(2, i)
          if (coords(2, i) < y_min) y_min = coords(2, i)
          if (coords(3, i) > z_max) z_max = coords(3, i)
          if (coords(3, i) < z_min) z_min = coords(3, i)
      end if
  end do
  
  face_area = (y_max - y_min) * (z_max - z_min)
  total_force = 100.0_prec * face_area ! 总拉力 = 应力 * 面积
  nodal_force = total_force / real(count_right, prec) ! 平均分配到面上各节点
  
  ! 第二遍循环：将拉力加载到 RHS，并对左侧进行位移固支
  do i = 1, n_nodes
      if (abs(coords(1, i) - 100.0_prec) < 1.0e-3_prec) then
          ! 右端：施加向右 (X正向) 的节点力
          rhs((i-1)*3 + 1) = rhs((i-1)*3 + 1) + nodal_force
      else if (abs(coords(1, i) - (-100.0_prec)) < 1.0e-3_prec) then
          ! 左端：对称移项法锁死 XYZ 自由度
          call apply_dirichlet_coo_sym(Kcoo, rhs, i, 1, 0.0_prec)
          call apply_dirichlet_coo_sym(Kcoo, rhs, i, 2, 0.0_prec)
          call apply_dirichlet_coo_sym(Kcoo, rhs, i, 3, 0.0_prec)
      end if
  end do
  ! =====================================================================

  print *, "      正在分发良态矩阵至 CSR 格式..."
  Kcsr = COO_to_CSR(Kcoo)

  ! --- 3. 求解器大比武循环 ---
  print *, ">>> [3/3] 启动现代求解器竞技..."
  print *, "---------------------------------------------------------------------------------"
  print *, " 求解器名称       | 耗时(s) | 孔顶X位移(mm) | 右上角X位移(mm) |  全场 L2 范数"
  print *, "---------------------------------------------------------------------------------"
  
  solver_names = ["Pardiso (Direct)", "CG (No Precond) ", "PCG (Jacobi)    "]
  
  do s = 1, 3
      x = 0.0_prec

      select case(s)
      case(1)
          opts%solver_family = "direct"
          call create_linear_solver(Kcsr, opts, solver)
      case(2)
          opts%solver_family = "cg" 
          opts%tol = 1.0e-7_prec
          opts%max_iter = 100000 
          call create_linear_solver(Kcsr, opts, solver)
      case(3)
          opts%solver_family = "iter" 
          opts%tol = 1.0e-7_prec
          opts%max_iter = 100000 
          call create_linear_solver(Kcsr, opts, solver)
      end select

      call cpu_time(t_start)
      call solver%analyze()
      call solver%factor()
      call solver%solve(rhs, x)
      call cpu_time(t_end)

      hole_apex_ux = x((hole_node_id-1)*3 + 1)
      corner_ux    = x((corner_node_id-1)*3 + 1)
      l2_norm      = sqrt(sum(x**2))

      write(*, '(A18, A, F7.4, A, F13.8, A, F15.8, A, F14.8)') &
          solver_names(s), " | ", t_end - t_start, " | ", hole_apex_ux, " | ", corner_ux, " | ", l2_norm

      call solver%free()
  end do
  print *, "---------------------------------------------------------------------------------"

contains
  subroutine apply_dirichlet_coo_sym(matK, f, node_id, dof, bnd_val)
      type(SparseCOO), intent(inout) :: matK
      real(prec), intent(inout) :: f(:)
      integer, intent(in) :: node_id, dof
      real(prec), intent(in) :: bnd_val
      integer :: eq_idx, k, r, c
      
      eq_idx = (node_id - 1) * 3 + dof
      
      do k = 1, matK%nnz
          r = matK%row(k)
          c = matK%col(k)
          
          if (c == eq_idx .and. r /= eq_idx) then
              f(r) = f(r) - matK%val(k) * bnd_val
              matK%val(k) = 0.0_prec
          else if (r == eq_idx .and. c /= eq_idx) then
              matK%val(k) = 0.0_prec
          else if (r == eq_idx .and. c == eq_idx) then
              matK%val(k) = 1.0_prec
          end if
      end do
      f(eq_idx) = bnd_val
  end subroutine apply_dirichlet_coo_sym

end program test_hole_plate_compare