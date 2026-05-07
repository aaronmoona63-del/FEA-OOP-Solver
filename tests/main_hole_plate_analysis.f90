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

  print *, "========================================================="
  print *, " ⚙️ OOFEM 计算引擎：纯净拉伸位移场全流程分析 ⚙️"
  print *, "========================================================="

  ! --- 1. 读取网格与材料属性 ---
  call load_mesh_from_txt('/work/tests/mesh_daikonglashen.txt', n_nodes, n_elements, coords, conn, props)
  n_unknowns = n_nodes * 3
  if (props(1) < 1.0_prec) props(1) = 210000.0_prec
  if (props(2) < 0.01_prec) props(2) = 0.3_prec

  ! --- 2. 全局刚度矩阵装配 ---
  print *, ">>> [1/4] 正在装配全局刚度矩阵..."
  call Kcoo%init(n_unknowns)
  call assemble_global_stiffness(Kcoo, n_elements, n_nodes, coords, conn, props)
  call Kcoo%preprocess() 
  
  ! --- 3. 施加纯正的 Neumann 分布面载荷 ---
  print *, ">>> [2/4] 正在依据单元面积积分，精确施加 100 MPa 拉伸力..."
  allocate(rhs(n_unknowns), x(n_unknowns))
  rhs = 0.0_prec

  do i = 1, n_elements
      f_nodes(:,1) = [1,4,3,2]; f_nodes(:,2) = [5,6,7,8]
      f_nodes(:,3) = [1,2,6,5]; f_nodes(:,4) = [3,4,8,7]
      f_nodes(:,5) = [2,3,7,6]; f_nodes(:,6) = [1,5,8,4]
      
      do f = 1, 6
          n1 = conn(f_nodes(1,f), i); n2 = conn(f_nodes(2,f), i)
          n3 = conn(f_nodes(3,f), i); n4 = conn(f_nodes(4,f), i)
          
          if (abs(coords(1,n1)-100.0_prec) < 1e-3_prec .and. &
              abs(coords(1,n2)-100.0_prec) < 1e-3_prec .and. &
              abs(coords(1,n3)-100.0_prec) < 1e-3_prec .and. &
              abs(coords(1,n4)-100.0_prec) < 1e-3_prec) then
              
              y_len = max(coords(2,n1), coords(2,n2), coords(2,n3), coords(2,n4)) - &
                      min(coords(2,n1), coords(2,n2), coords(2,n3), coords(2,n4))
              z_len = max(coords(3,n1), coords(3,n2), coords(3,n3), coords(3,n4)) - &
                      min(coords(3,n1), coords(3,n2), coords(3,n3), coords(3,n4))
              
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

  ! --- 4. 启动自研 PCG 求解器求解 ---
  print *, ">>> [3/4] 启动自研 PCG 求解器，求解未知位移场..."
  x = 0.0_prec
  opts%solver_family = "iter" 
  opts%tol = 1.0e-7_prec
  opts%max_iter = 100000 
  call create_linear_solver(Kcsr, opts, solver)
  
  call cpu_time(t_start); call solver%analyze(); call solver%factor(); call solver%solve(rhs, x); call cpu_time(t_end)
  print *, "      ✅ 位移场求解完成！耗时: ", t_end - t_start, " 秒"

  ! --- 5. 纯净位移导出 ---
  print *, ">>> [4/4] 正在导出纯净位移场至 VTK 文件..."
  call export_vtk_displacement("hole_plate_displacement.vtk", n_nodes, n_elements, coords, conn, x)
  
  print *, "========================================================="
  print *, " 🎉 完美竣工！结果已保存至：build/hole_plate_displacement.vtk 🎉"
  print *, "========================================================="
  ! =====================================================================
  ! 📊 专属定制：提取中轴线数据用于 Python 绘图
  ! =====================================================================
  print *, ">>> [DATA] 正在提取路径曲线数据至 txt 文件..."
  
  ! 1. 提取顶边 X 方向位移 (Y=50, Z=5, 沿 X 轴)
  open(unit=88, file='path_disp_ux.txt', status='replace')
  do i = 1, n_nodes
      if (abs(coords(2, i) - 50.0_prec) < 1.0e-3_prec .and. abs(coords(3, i) - 5.0_prec) < 1.0e-3_prec) then
          write(88, '(2F15.6)') coords(1, i), x((i-1)*3+1)  ! 写入: X坐标, Ux位移
      end if
  end do
  close(88)

  print *, "      ✅ 中轴线位移数据已保存至：build/path_disp_ux.txt"
  ! =====================================================================

  call solver%free()

contains
  ! 对称移项边界处理
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

  ! 纯位移 VTK 导出 (无应力计算)
  subroutine export_vtk_displacement(filename, nn, ne, coords, conn, u)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: nn, ne
      real(prec), intent(in) :: coords(:,:), u(:)
      integer, intent(in) :: conn(:,:)
      integer :: i, vtk_unit = 50

      open(unit=vtk_unit, file=filename, status='replace', action='write')
      write(vtk_unit, '(A)') "# vtk DataFile Version 3.0"
      write(vtk_unit, '(A)') "OOP FEM Displacement Results"
      write(vtk_unit, '(A)') "ASCII"
      write(vtk_unit, '(A)') "DATASET UNSTRUCTURED_GRID"
      
      write(vtk_unit, '(A, I0, A)') "POINTS ", nn, " float"
      do i = 1, nn
          write(vtk_unit, '(3E15.6)') coords(1,i), coords(2,i), coords(3,i)
      end do
      
      write(vtk_unit, '(A, I0, X, I0)') "CELLS ", ne, ne * 9
      do i = 1, ne
          write(vtk_unit, '(I2, 8I8)') 8, conn(1:8, i) - 1
      end do
      
      write(vtk_unit, '(A, I0)') "CELL_TYPES ", ne
      do i = 1, ne
          write(vtk_unit, '(I3)') 12 
      end do
      
      write(vtk_unit, '(A, I0)') "POINT_DATA ", nn
      write(vtk_unit, '(A)') "VECTORS Displacement float"
      do i = 1, nn
          write(vtk_unit, '(3E15.6)') u((i-1)*3+1), u((i-1)*3+2), u((i-1)*3+3)
      end do
      
      close(vtk_unit)
  end subroutine export_vtk_displacement
end program main_hole_plate_analysis