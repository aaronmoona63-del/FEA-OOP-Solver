program main_hole_plate_analysis
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
  integer :: i, f, n1, n2, n3, n4
  integer :: f_nodes(4, 6)
  real(prec) :: props(2)
  real(prec) :: y_len, z_len, face_area, nodal_force
  real(8) :: t_start, t_end
  real(8) :: time_pardiso, time_cg, time_pcg

  print *, "========================================================="
  print *, " 🟥 [实验组] 带孔畸变板全自动多求解器大比武 🟥"
  print *, "========================================================="

  ! 读取带孔板网格
  call load_mesh_from_txt('/work/tests/mesh_daikonglashen.txt', n_nodes, n_elements, coords, conn, props)
  n_unknowns = n_nodes * 3
  if (props(1) < 1.0_prec) props(1) = 210000.0_prec
  if (props(2) < 0.01_prec) props(2) = 0.3_prec

  print *, ">>> [1/3] 正在装配带孔板全局刚度矩阵..."
  call Kcoo%init(n_unknowns)
  call assemble_global_stiffness(Kcoo, n_elements, n_nodes, coords, conn, props)
  call Kcoo%preprocess() 
  
  allocate(rhs(n_unknowns), x(n_unknowns))
  rhs = 0.0_prec

  do i = 1, n_elements
      f_nodes(:,1) = [1,4,3,2]; f_nodes(:,2) = [5,6,7,8]
      f_nodes(:,3) = [1,2,6,5]; f_nodes(:,4) = [3,4,8,7]
      f_nodes(:,5) = [2,3,7,6]; f_nodes(:,6) = [1,5,8,4]
      do f = 1, 6
          n1 = conn(f_nodes(1,f), i); n2 = conn(f_nodes(2,f), i)
          n3 = conn(f_nodes(3,f), i); n4 = conn(f_nodes(4,f), i)
          if (abs(coords(1,n1)-100.0_prec) < 1e-3_prec .and. abs(coords(1,n2)-100.0_prec) < 1e-3_prec .and. &
              abs(coords(1,n3)-100.0_prec) < 1e-3_prec .and. abs(coords(1,n4)-100.0_prec) < 1e-3_prec) then
              y_len = max(coords(2,n1), coords(2,n2), coords(2,n3), coords(2,n4)) - min(coords(2,n1), coords(2,n2), coords(2,n3), coords(2,n4))
              z_len = max(coords(3,n1), coords(3,n2), coords(3,n3), coords(3,n4)) - min(coords(3,n1), coords(3,n2), coords(3,n3), coords(3,n4))
              face_area = y_len * z_len
              nodal_force = (100.0_prec * face_area) / 4.0_prec
              rhs((n1-1)*3+1) = rhs((n1-1)*3+1) + nodal_force
              rhs((n2-1)*3+1) = rhs((n2-1)*3+1) + nodal_force
              rhs((n3-1)*3+1) = rhs((n3-1)*3+1) + nodal_force
              rhs((n4-1)*3+1) = rhs((n4-1)*3+1) + nodal_force
          end if
      end do
  end do

  do i = 1, n_nodes
      if (abs(coords(1, i) - (-100.0_prec)) < 1.0e-3_prec) then
          call apply_dirichlet_coo_sym(Kcoo, rhs, i, 1, 0.0_prec)
          call apply_dirichlet_coo_sym(Kcoo, rhs, i, 2, 0.0_prec)
          call apply_dirichlet_coo_sym(Kcoo, rhs, i, 3, 0.0_prec)
      end if
  end do

  Kcsr = COO_to_CSR(Kcoo)
  opts%tol = 1.0e-7_prec
  opts%max_iter = 100000 

  print *, " "
  print *, ">>> [2/3] 🚀 开始求解器大比武..."

  ! ====== 战将 1：Intel MKL Pardiso ======
  x = 0.0_prec
  opts%solver_family = "direct"
  call create_linear_solver(Kcsr, opts, solver)
  call cpu_time(t_start); call solver%analyze(); call solver%factor(); call solver%solve(rhs, x); call cpu_time(t_end)
  time_pardiso = t_end - t_start
  print *, "    [1] Pardiso (Direct) 耗时: ", time_pardiso, " 秒"
  call solver%free(); deallocate(solver)

  ! ====== 战将 2：标准 CG ======
  x = 0.0_prec  
  opts%solver_family = "cg"
  call create_linear_solver(Kcsr, opts, solver)
  call cpu_time(t_start); call solver%analyze(); call solver%factor(); call solver%solve(rhs, x); call cpu_time(t_end)
  time_cg = t_end - t_start
  print *, "    [2] Standard CG      耗时: ", time_cg, " 秒"
  call solver%free(); deallocate(solver)

  ! ====== 战将 3：Jacobi PCG ======
  x = 0.0_prec  
  opts%solver_family = "iter" 
  call create_linear_solver(Kcsr, opts, solver)
  call cpu_time(t_start); call solver%analyze(); call solver%factor(); call solver%solve(rhs, x); call cpu_time(t_end)
  time_pcg = t_end - t_start
  print *, "    [3] Jacobi PCG       耗时: ", time_pcg, " 秒"
  call solver%free(); deallocate(solver)

  print *, "========================================================="
  print *, " 🎉 带孔板测试结束！请记录屏幕上的时间和迭代次数！"
  print *, "========================================================="

contains
  subroutine apply_dirichlet_coo_sym(matK, f, node_id, dof, bnd_val)
      type(SparseCOO), intent(inout) :: matK
      real(prec), intent(inout) :: f(:)
      integer, intent(in) :: node_id, dof
      real(prec), intent(in) :: bnd_val
      integer :: eq_idx, k, r, c
      eq_idx = (node_id - 1) * 3 + dof
      do k = 1, matK%nnz
          r = matK%row(k); c = matK%col(k)
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
end program main_hole_plate_analysis