program main_hole_400k
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
  real(prec) :: props(2), y_len, z_len, face_area, nodal_force
  real(8) :: t_start, t_end
  real(prec) :: ux_val
  real(prec) :: x_min, x_max

  ! ================= 用于记录性能和位移的文件名 =================
  character(len=50) :: txt_filename, disp_filename
  character(len=20) :: current_mesh_name
  real(8)           :: memory_coo_mb, memory_csr_mb
  real(8)           :: solve_duration
  ! ==============================================================

  print *, "========================================================="
  print *, " 🟧 [次大规模] ~40万自由度带孔板 标度律黄金点核验 🟧"
  print *, "========================================================="

  ! 读取 40 万网格文件
  call load_mesh_from_txt('/work/tests/mesh_hole_400k.txt', n_nodes, n_elements, coords, conn, props)
  n_unknowns = n_nodes * 3
  if (props(1) < 1.0_prec) props(1) = 210000.0_prec
  if (props(2) < 0.01_prec) props(2) = 0.3_prec

  x_min = minval(coords(1, :))
  x_max = maxval(coords(1, :))
  print *, "🔍 自动探测几何边界成功: X_min = ", x_min, " | X_max = ", x_max

  print *, ">>> [1/2] 正在装配全局刚度矩阵..."
  call Kcoo%init(n_unknowns)
  call assemble_global_stiffness(Kcoo, n_elements, n_nodes, coords, conn, props)
  call Kcoo%preprocess() 
  
  allocate(rhs(n_unknowns), x(n_unknowns))
  rhs = 0.0_prec

  ! 施加拉伸载荷 (右端面 100 MPa)
  do i = 1, n_elements
      f_nodes(:,1) = [1,4,3,2]; f_nodes(:,2) = [5,6,7,8]
      f_nodes(:,3) = [1,2,6,5]; f_nodes(:,4) = [3,4,8,7]
      f_nodes(:,5) = [2,3,7,6]; f_nodes(:,6) = [1,5,8,4]
      do f = 1, 6
          n1 = conn(f_nodes(1,f), i); n2 = conn(f_nodes(2,f), i)
          n3 = conn(f_nodes(3,f), i); n4 = conn(f_nodes(4,f), i)
          if (abs(coords(1,n1)-x_max) < 1e-3_prec .and. abs(coords(1,n2)-x_max) < 1e-3_prec .and. &
              abs(coords(1,n3)-x_max) < 1e-3_prec .and. abs(coords(1,n4)-x_max) < 1e-3_prec) then
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

  ! 施加固支边界 (左端面全固定)
  do i = 1, n_nodes
      if (abs(coords(1, i) - x_min) < 1.0e-3_prec) then
          call apply_dirichlet_coo_sym(Kcoo, rhs, i, 1, 0.0_prec)
          call apply_dirichlet_coo_sym(Kcoo, rhs, i, 2, 0.0_prec)
          call apply_dirichlet_coo_sym(Kcoo, rhs, i, 3, 0.0_prec)
      end if
  end do

  Kcsr = COO_to_CSR(Kcoo)
  opts%tol = 1.0e-6_prec
  opts%max_iter = 100000 
  
  ! 设定日志与位移输出的文件名
  txt_filename = "Chapter4_Benchmark_Log.txt"
  disp_filename = "displacement_hole_400k.txt"  ! <--- 这里是新加的位移结果文件名
  current_mesh_name = "Hole_400k"
  
  memory_coo_mb = real(Kcoo%nnz, 8) * 16.0_8 / 1048576.0_8
  memory_csr_mb = real(Kcsr%nnz * 12 + (n_unknowns + 1) * 4, 8) / 1048576.0_8

  print *, ">>> [2/2] 🚀 启动大比武..."

  ! 我们以 Pardiso 的结果作为标准输出
  print *, "    [1] Intel MKL Pardiso (直接法)..."
  x = 0.0_prec
  opts%solver_family = "direct"
  call create_linear_solver(Kcsr, opts, solver)
  call cpu_time(t_start); call solver%analyze(); call solver%factor(); call solver%solve(rhs, x); call cpu_time(t_end)
  ux_val = get_right_face_ux(n_nodes, coords, x, x_max)
  solve_duration = t_end - t_start
  print *, "    ✅ 耗时: ", solve_duration, " 秒 | Ux: ", ux_val, " mm"
  
  ! --- 记录性能日志 ---
  call write_benchmark_log(txt_filename, current_mesh_name, "Pardiso", n_unknowns, Kcsr%nnz, solve_duration, memory_coo_mb, memory_csr_mb)
  
  ! --- [核心新增] 输出所有节点的位移到 TXT 文件 ---
  print *, "    💾 正在导出全节点位移数据至 TXT 文件..."
  call write_nodal_displacements(disp_filename, n_nodes, coords, x)
  print *, "    ✅ 导出完成！文件已保存为: ", trim(disp_filename)
  
  call solver%free(); deallocate(solver)

  ! === 后面的 CG 和 PCG 省略位移导出，只记录求解时间 ===
  ! 2: CG
  print *, "    [2] Standard CG (无预处理)..."
  x = 0.0_prec  
  opts%solver_family = "cg"
  call create_linear_solver(Kcsr, opts, solver)
  call cpu_time(t_start); call solver%analyze(); call solver%factor(); call solver%solve(rhs, x); call cpu_time(t_end)
  solve_duration = t_end - t_start
  call write_benchmark_log(txt_filename, current_mesh_name, "Std_CG", n_unknowns, Kcsr%nnz, solve_duration, memory_coo_mb, memory_csr_mb)
  call solver%free(); deallocate(solver)

  ! 3: PCG
  print *, "    [3] Jacobi PCG (对角预处理)..."
  x = 0.0_prec  
  opts%solver_family = "iter" 
  call create_linear_solver(Kcsr, opts, solver)
  call cpu_time(t_start); call solver%analyze(); call solver%factor(); call solver%solve(rhs, x); call cpu_time(t_end)
  solve_duration = t_end - t_start
  call write_benchmark_log(txt_filename, current_mesh_name, "Jacobi_PCG", n_unknowns, Kcsr%nnz, solve_duration, memory_coo_mb, memory_csr_mb)
  call solver%free(); deallocate(solver)

contains

  ! ==================================================================
  ! 新增：将全网格节点坐标与位移输出到 TXT 文件的子程序
  ! ==================================================================
  subroutine write_nodal_displacements(filename, nn, c, u)
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(in)          :: nn
      real(prec), intent(in)       :: c(:,:), u(:)
      integer :: i, f_unit

      f_unit = 110
      open(unit=f_unit, file=filename, status='replace', action='write')

      ! 写入完美的对齐表头
      write(f_unit, '(A8,6A16)') "Node_ID", "X_Coord", "Y_Coord", "Z_Coord", "U_x", "U_y", "U_z"

      ! 遍历所有节点写入坐标和 XYZ 三个方向的位移 (科学计数法保留6位小数)
      do i = 1, nn
          write(f_unit, '(I8,6ES16.6)') i, c(1,i), c(2,i), c(3,i), &
                                        u((i-1)*3+1), u((i-1)*3+2), u((i-1)*3+3)
      end do

      close(f_unit)
  end subroutine write_nodal_displacements


  ! (原有的记录性能日志的子程序)
  subroutine write_benchmark_log(filename, mesh_name, solver_name, dofs, nnz, time_solve, mem_coo, mem_csr)
      implicit none
      character(len=*), intent(in) :: filename, mesh_name, solver_name
      integer, intent(in)          :: dofs, nnz
      real(8), intent(in)          :: time_solve, mem_coo, mem_csr
      integer :: file_unit = 105
      logical :: file_exists

      inquire(file=filename, exist=file_exists)
      open(unit=file_unit, file=filename, status='unknown', position='append', action='write')

      if (.not. file_exists) then
          write(file_unit, '(A)') "=========================================================================================================="
          write(file_unit, '(A)') "                                       第四章 大规模有限元基准测试性能日志                                "
          write(file_unit, '(A)') "=========================================================================================================="
          write(file_unit, '(A20,A15,A15,A15,A15,A15,A15)') "Mesh_Name", "Solver", "DOFs", "NNZ", "Solve_Time(s)", "COO_Mem(MB)", "CSR_Mem(MB)"
          write(file_unit, '(A)') "----------------------------------------------------------------------------------------------------------"
      end if

      write(file_unit, '(A20,A15,I15,I15,F15.6,F15.2,F15.2)') &
          trim(mesh_name), trim(solver_name), dofs, nnz, time_solve, mem_coo, mem_csr

      close(file_unit)
  end subroutine write_benchmark_log


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

  function get_right_face_ux(nn, c, u, xmax) result(avg_ux)
      integer, intent(in) :: nn
      real(prec), intent(in) :: c(:,:), u(:)
      real(prec), intent(in) :: xmax
      real(prec) :: avg_ux
      integer :: i, cnt
      avg_ux = 0.0_prec
      cnt = 0
      do i = 1, nn
          if (abs(c(1, i) - xmax) < 1.0e-3_prec) then
              avg_ux = avg_ux + u((i-1)*3 + 1)
              cnt = cnt + 1
          end if
      end do
      if (cnt > 0) avg_ux = avg_ux / real(cnt, prec)
  end function get_right_face_ux

end program main_hole_400k