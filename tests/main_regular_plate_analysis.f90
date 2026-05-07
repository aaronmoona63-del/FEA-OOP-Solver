program main_regular_plate_analysis
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
  
  ! 专用于记录三个求解器的时间
  real(8) :: time_pardiso, time_cg, time_pcg

  print *, "========================================================="
  print *, " 🟩 [对照组] 规则板全自动多求解器大比武 (One-Click Benchmark) 🟩"
  print *, "========================================================="

  call load_mesh_from_txt('/work/tests/mesh_regular_plate.txt', n_nodes, n_elements, coords, conn, props)
  n_unknowns = n_nodes * 3
  if (props(1) < 1.0_prec) props(1) = 210000.0_prec
  if (props(2) < 0.01_prec) props(2) = 0.3_prec

  print *, ">>> [1/4] 正在装配规则网格全局刚度矩阵 (只需装配一次!)..."
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
  print *, ">>> [2/4] 🚀 开始求解器大比武..."

  ! ====== 战将 1：Intel MKL Pardiso (直接法) ======
  x = 0.0_prec  ! 务必清零初始位移
  opts%solver_family = "direct"
  call create_linear_solver(Kcsr, opts, solver)
  call cpu_time(t_start); call solver%analyze(); call solver%factor(); call solver%solve(rhs, x); call cpu_time(t_end)
  time_pardiso = t_end - t_start
  print *, "    [1] Pardiso (Direct) 耗时: ", time_pardiso, " 秒"
  call solver%free()
  deallocate(solver) ! 释放内存，为下一个求解器腾出空间

  ! ====== 战将 2：标准 CG (无预处理) ======
  x = 0.0_prec  
  opts%solver_family = "cg"
  call create_linear_solver(Kcsr, opts, solver)
  call cpu_time(t_start); call solver%analyze(); call solver%factor(); call solver%solve(rhs, x); call cpu_time(t_end)
  time_cg = t_end - t_start
  print *, "    [2] Standard CG      耗时: ", time_cg, " 秒"
  call solver%free()
  deallocate(solver)

  ! ====== 战将 3：Jacobi PCG (自研预处理) ======
  x = 0.0_prec  
  opts%solver_family = "iter" 
  call create_linear_solver(Kcsr, opts, solver)
  call cpu_time(t_start); call solver%analyze(); call solver%factor(); call solver%solve(rhs, x); call cpu_time(t_end)
  time_pcg = t_end - t_start
  print *, "    [3] Jacobi PCG       耗时: ", time_pcg, " 秒"
  
  ! 导出最后一次跑完的位移结果
  call export_vtk_displacement("regular_plate_displacement.vtk", n_nodes, n_elements, coords, conn, x)
  call solver%free()
  deallocate(solver)

  ! --- 3. 自动将耗时写入文本，供 Python 读取 ---
  print *, ">>> [3/4] 正在将比武数据写入 benchmark_time.txt..."
  open(unit=99, file='benchmark_time.txt', status='replace')
  write(99, '(F10.6)') time_pardiso
  write(99, '(F10.6)') time_cg
  write(99, '(F10.6)') time_pcg
  close(99)

  print *, "========================================================="
  print *, " 🎉 自动化测试竣工！数据已保存，随时可运行 Python 绘图！"
  print *, "========================================================="

contains
  ! (保持原有的 apply_dirichlet_coo_sym 和 export_vtk_displacement 不变)
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
end program main_regular_plate_analysis