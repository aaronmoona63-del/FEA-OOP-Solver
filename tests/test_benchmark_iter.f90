program test_benchmark_iter
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
  integer :: i
  real(prec) :: props(2)

  ! === 新增：用于存储二维格式位移的数组 ===
  real(8), allocatable :: nodalDisps(:,:)
  ! ========================================

  real(8) :: t_start, t_end, t_assembly, t_convert, t_solve

  print *, "--- [OOP FEM Engine: ITERATIVE Solver Performance Test] ---"
  
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

  ! --- 3. 边界与载荷施加阶段 ---
  allocate(rhs(n_unknowns), x(n_unknowns))
  call apply_benchmark_loads(n_nodes, coords, rhs)
  call apply_benchmark_bcs_csr(n_nodes, coords, Kcsr, rhs)

  ! --- 4. 迭代法求解阶段 ---
  call cpu_time(t_start)
  
  opts%solver_family = "iter"
  
  call create_linear_solver(Kcsr, opts, solver)
  call solver%analyze(); call solver%factor(); call solver%solve(rhs, x)
  
  call cpu_time(t_end)
  t_solve = t_end - t_start

  ! 检验解的正确性
  print *, " 顶部受拉点位移测试值: ", x(n_unknowns)

  ! 打印压测报告
  print *, "========================================="
  print *, "       迭代法 (Iterative) 性能报告       "
  print *, "========================================="
  print *, " 单元数量 : ", n_elements
  print *, " 总自由度 : ", n_unknowns
  print *, " 非零元素 : ", Kcsr%nnz
  print *, " 装配耗时 : ", t_assembly, " 秒"
  print *, " 转换耗时 : ", t_convert,  " 秒"
  print *, " 求解耗时 : ", t_solve,    " 秒"
  print *, "========================================="

  call solver%free()

  ! =========================================================
  ! 🚀 新增：将全场位移导出为 VTK 文件，用于 ParaView 渲染！
  ! =========================================================
  print *, "=> 正在提取全场位移并生成 VTK 文件..."
  
  ! 分配二维位移数组 (节点数, 3个方向)
  allocate(nodalDisps(n_nodes, 3))
  
  ! 将求解器算出的一维向量 x 重组为二维格式
  do i = 1, n_nodes
      nodalDisps(i, 1) = x(3*i - 2) ! Ux
      nodalDisps(i, 2) = x(3*i - 1) ! Uy
      nodalDisps(i, 3) = x(3*i)     ! Uz
  end do
  
  ! 调用内部子程序生成 VTK
  call ExportToVTK('displacement_hole_135k.vtk', n_nodes, n_elements, coords, conn, nodalDisps)
  
  ! 释放临时数组
  deallocate(nodalDisps)

  ! =========================================================
  ! 原有代码：提取中心轴 Z 方向位移
  ! =========================================================
  print *, "=> 正在提取中心轴节点位移数据到 my_solver_uz.txt ..."
  open(unit=99, file='my_solver_uz.txt', status='replace')
  do i = 1, n_nodes
      if (abs(coords(1, i) - 5.0_prec) < 1.0e-4_prec .and. &
          abs(coords(2, i) - 5.0_prec) < 1.0e-4_prec) then
          write(99, *) coords(3, i), x(3*i)
      end if
  end do
  close(99)
  print *, "=> 提取完成！"

contains

  ! ======================================================================
  ! 内部子程序：将 Hex8 网格和位移结果输出为 ParaView 可读的 VTK 文件
  ! 注意：适配了你的 coords(3, N) 和 conn(8, N) 内存排列方式
  ! ======================================================================
  subroutine ExportToVTK(filename, numNodes, numElems, nodes, elements, U)
      character(len=*), intent(in) :: filename
      integer, intent(in)          :: numNodes, numElems
      real(8), intent(in)          :: nodes(3, numNodes)    ! (X,Y,Z), N
      integer, intent(in)          :: elements(8, numElems) ! (8个节点), E
      real(8), intent(in)          :: U(numNodes, 3)        ! N, (Ux,Uy,Uz)
      
      integer :: j, fileUnit
      
      fileUnit = 88
      open(unit=fileUnit, file=filename, status='replace', action='write')
      
      ! 1. Header
      write(fileUnit, '(A)') '# vtk DataFile Version 3.0'
      write(fileUnit, '(A)') 'FEM 3D Hex8 Displacement Results'
      write(fileUnit, '(A)') 'ASCII'
      write(fileUnit, '(A)') 'DATASET UNSTRUCTURED_GRID'
      
      ! 2. 节点坐标 POINTS
      write(fileUnit, '(A, I0, A)') 'POINTS ', numNodes, ' float'
      do j = 1, numNodes
          write(fileUnit, '(3(E15.6, 1X))') nodes(1, j), nodes(2, j), nodes(3, j)
      end do
      
      ! 3. 单元连接 CELLS
      write(fileUnit, '(A, I0, 1X, I0)') 'CELLS ', numElems, numElems * 9
      do j = 1, numElems
          ! VTK 节点索引从 0 开始，Fortran 从 1 开始，必须减 1
          write(fileUnit, '(I2, 8(1X, I0))') 8, elements(1, j)-1, elements(2, j)-1, &
                                                elements(3, j)-1, elements(4, j)-1, &
                                                elements(5, j)-1, elements(6, j)-1, &
                                                elements(7, j)-1, elements(8, j)-1
      end do
      
      ! 4. 单元类型 CELL_TYPES
      write(fileUnit, '(A, I0)') 'CELL_TYPES ', numElems
      do j = 1, numElems
          write(fileUnit, '(I2)') 12  ! 12 代表 VTK_HEXAHEDRON
      end do
      
      ! 5. 位移场 POINT_DATA
      write(fileUnit, '(A, I0)') 'POINT_DATA ', numNodes
      write(fileUnit, '(A)') 'VECTORS Displacement float'
      do j = 1, numNodes
          write(fileUnit, '(3(E15.6, 1X))') U(j, 1), U(j, 2), U(j, 3)
      end do
      
      close(fileUnit)
      print *, "✅ [成功] 完整三维 VTK 文件已保存为: ", trim(filename)

  end subroutine ExportToVTK

end program test_benchmark_iter