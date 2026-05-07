program test_two_elements_3d
  use Types, only : prec
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

  real(prec), allocatable :: rhs(:), x(:)
  real(prec) :: nodes(3, 12), E, nu, c1, detJ, weight
  real(prec) :: D(6, 6), Ke(24, 24), B(6, 24), dNdxi(3, 8), dNdx(3, 8)
  real(prec) :: Jac(3, 3), invJac(3, 3), temp_DB(6, 24)
  real(prec) :: K_global(36, 36)
  real(prec) :: xi, eta, zeta, big_number
  real(prec) :: gp_loc(2)
  
  ! 应力后处理
  real(prec) :: de(24), temp_stress(6), sigma(6)
  
  integer :: el_cnn(8, 2), fixNodes(4)
  integer :: n_nodes, n_elements, n_unknowns
  integer :: i, j, el, r, c, rr, cc, gp1, gp2, gp3, nid
  
  print *, '--- [MATLAB Benchmark: 3D Hex8 Elements Tension] ---'
  
  n_nodes = 12
  n_elements = 2
  n_unknowns = n_nodes * 3
  
  ! 节点坐标 (与 MATLAB 完全一致)
  nodes(1,:) = [1000.0_prec, 1000.0_prec,    0.0_prec,    0.0_prec, 1000.0_prec,    0.0_prec,    0.0_prec, 1000.0_prec, 2000.0_prec, 2000.0_prec, 2000.0_prec, 2000.0_prec]
  nodes(2,:) = [   0.0_prec,    0.0_prec,    0.0_prec,    0.0_prec, 1000.0_prec, 2000.0_prec, 2000.0_prec, 1000.0_prec, 1000.0_prec, 1000.0_prec,    0.0_prec,    0.0_prec]
  nodes(3,:) = [   0.0_prec, 1000.0_prec, 1000.0_prec,    0.0_prec, 1000.0_prec, 1000.0_prec,    0.0_prec,    0.0_prec, 1000.0_prec,    0.0_prec, 1000.0_prec,    0.0_prec]
  
  ! 单元连结 (8节点)
  el_cnn(:,1) = [8, 5, 6, 7, 1, 2, 3, 4]
  el_cnn(:,2) = [2, 11, 12, 1, 5, 9, 10, 8]
  
  ! 材料参数
  E = 210000.0_prec
  nu = 0.3_prec
  
  ! 3D 弹性本构矩阵 D (6x6)
  c1 = E / ((1.0_prec + nu) * (1.0_prec - 2.0_prec * nu))
  D = 0.0_prec
  D(1,1) = c1 * (1.0_prec - nu); D(1,2) = c1 * nu;               D(1,3) = c1 * nu
  D(2,1) = c1 * nu;              D(2,2) = c1 * (1.0_prec - nu);  D(2,3) = c1 * nu
  D(3,1) = c1 * nu;              D(3,2) = c1 * nu;               D(3,3) = c1 * (1.0_prec - nu)
  D(4,4) = c1 * (1.0_prec - 2.0_prec * nu) / 2.0_prec  ! G 剪切模量
  D(5,5) = D(4,4); D(6,6) = D(4,4)

  K_global = 0.0_prec
  allocate(rhs(n_unknowns), x(n_unknowns))
  rhs = 0.0_prec
  x = 0.0_prec
  
  gp_loc = [-1.0_prec/sqrt(3.0_prec), 1.0_prec/sqrt(3.0_prec)]
  
  ! =========================================================
  ! 1. 组装 3D 全局刚度矩阵
  ! =========================================================
  do el = 1, n_elements
      Ke = 0.0_prec
      do gp1 = 1, 2
          do gp2 = 1, 2
              do gp3 = 1, 2
                  xi = gp_loc(gp1); eta = gp_loc(gp2); zeta = gp_loc(gp3)
                  weight = 1.0_prec * 1.0_prec * 1.0_prec
                  
                  ! 形函数偏导数 (3 x 8)
                  dNdxi(1,:) = [ -(1.0_prec-eta)*(1.0_prec-zeta)/8.0_prec,  (1.0_prec-eta)*(1.0_prec-zeta)/8.0_prec,  (1.0_prec+eta)*(1.0_prec-zeta)/8.0_prec, -(1.0_prec+eta)*(1.0_prec-zeta)/8.0_prec, &
                                 -(1.0_prec-eta)*(1.0_prec+zeta)/8.0_prec,  (1.0_prec-eta)*(1.0_prec+zeta)/8.0_prec,  (1.0_prec+eta)*(1.0_prec+zeta)/8.0_prec, -(1.0_prec+eta)*(1.0_prec+zeta)/8.0_prec ]
                  dNdxi(2,:) = [ -(1.0_prec-xi)*(1.0_prec-zeta)/8.0_prec,  -(1.0_prec+xi)*(1.0_prec-zeta)/8.0_prec,   (1.0_prec+xi)*(1.0_prec-zeta)/8.0_prec,   (1.0_prec-xi)*(1.0_prec-zeta)/8.0_prec,  &
                                 -(1.0_prec-xi)*(1.0_prec+zeta)/8.0_prec,  -(1.0_prec+xi)*(1.0_prec+zeta)/8.0_prec,   (1.0_prec+xi)*(1.0_prec+zeta)/8.0_prec,   (1.0_prec-xi)*(1.0_prec+zeta)/8.0_prec ]
                  dNdxi(3,:) = [ -(1.0_prec-xi)*(1.0_prec-eta)/8.0_prec,   -(1.0_prec+xi)*(1.0_prec-eta)/8.0_prec,   -(1.0_prec+xi)*(1.0_prec+eta)/8.0_prec,   -(1.0_prec-xi)*(1.0_prec+eta)/8.0_prec,   &
                                  (1.0_prec-xi)*(1.0_prec-eta)/8.0_prec,    (1.0_prec+xi)*(1.0_prec-eta)/8.0_prec,    (1.0_prec+xi)*(1.0_prec+eta)/8.0_prec,    (1.0_prec-xi)*(1.0_prec+eta)/8.0_prec ]
                  
                  ! 3x3 雅可比矩阵
                  Jac = 0.0_prec
                  do i = 1, 8
                      do rr = 1, 3
                          do cc = 1, 3
                              Jac(rr,cc) = Jac(rr,cc) + dNdxi(rr,i) * nodes(cc, el_cnn(i,el))
                          end do
                      end do
                  end do
                  
                  ! 行列式与求逆 (3x3)
                  detJ = Jac(1,1)*(Jac(2,2)*Jac(3,3) - Jac(2,3)*Jac(3,2)) &
                       - Jac(1,2)*(Jac(2,1)*Jac(3,3) - Jac(2,3)*Jac(3,1)) &
                       + Jac(1,3)*(Jac(2,1)*Jac(3,2) - Jac(2,2)*Jac(3,1))
                       
                  invJac(1,1) =  (Jac(2,2)*Jac(3,3) - Jac(2,3)*Jac(3,2)) / detJ
                  invJac(1,2) = -(Jac(1,2)*Jac(3,3) - Jac(1,3)*Jac(3,2)) / detJ
                  invJac(1,3) =  (Jac(1,2)*Jac(2,3) - Jac(1,3)*Jac(2,2)) / detJ
                  invJac(2,1) = -(Jac(2,1)*Jac(3,3) - Jac(2,3)*Jac(3,1)) / detJ
                  invJac(2,2) =  (Jac(1,1)*Jac(3,3) - Jac(1,3)*Jac(3,1)) / detJ
                  invJac(2,3) = -(Jac(1,1)*Jac(2,3) - Jac(1,3)*Jac(2,1)) / detJ
                  invJac(3,1) =  (Jac(2,1)*Jac(3,2) - Jac(2,2)*Jac(3,1)) / detJ
                  invJac(3,2) = -(Jac(1,1)*Jac(3,2) - Jac(1,2)*Jac(3,1)) / detJ
                  invJac(3,3) =  (Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1)) / detJ
                  
                  dNdx = matmul(invJac, dNdxi)
                  
                  ! 组装 6x24 的 B 矩阵
                  B = 0.0_prec
                  do i = 1, 8
                      B(1, 3*i-2) = dNdx(1,i); B(1, 3*i-1) = 0.0_prec;  B(1, 3*i)   = 0.0_prec
                      B(2, 3*i-2) = 0.0_prec;  B(2, 3*i-1) = dNdx(2,i); B(2, 3*i)   = 0.0_prec
                      B(3, 3*i-2) = 0.0_prec;  B(3, 3*i-1) = 0.0_prec;  B(3, 3*i)   = dNdx(3,i)
                      
                      B(4, 3*i-2) = dNdx(2,i); B(4, 3*i-1) = dNdx(1,i); B(4, 3*i)   = 0.0_prec
                      B(5, 3*i-2) = 0.0_prec;  B(5, 3*i-1) = dNdx(3,i); B(5, 3*i)   = dNdx(2,i)
                      B(6, 3*i-2) = dNdx(3,i); B(6, 3*i-1) = 0.0_prec;  B(6, 3*i)   = dNdx(1,i)
                  end do
                  
                  temp_DB = matmul(D, B)
                  Ke = Ke + matmul(transpose(B), temp_DB) * detJ * weight
              end do
          end do
      end do
      
      ! 组装到稠密全局矩阵
      do i = 1, 8
          do j = 1, 8
              r = 3*(el_cnn(i,el)-1)
              c = 3*(el_cnn(j,el)-1)
              do rr = 1, 3
                  do cc = 1, 3
                      K_global(r+rr, c+cc) = K_global(r+rr, c+cc) + Ke(3*i-3+rr, 3*j-3+cc)
                  end do
              end do
          end do
      end do
  end do

  ! =========================================================
  ! 2. 边界条件施加 (3D)
  ! =========================================================
  ! Neumann 边界: 在端面 (节点 11, 9, 10, 12) 施加 X 方向面力
  ! 总面积 10^6 mm^2，受力 100 MPa，总力 10^8 N，每个节点分担 2.5e7 N
  rhs(3*11 - 2) = 25000000.0_prec
  rhs(3*9  - 2) = 25000000.0_prec
  rhs(3*10 - 2) = 25000000.0_prec
  rhs(3*12 - 2) = 25000000.0_prec

  ! Dirichlet 边界: 固定节点 3, 4, 6, 7
  fixNodes = [3, 4, 6, 7]
  big_number = 1.0e15_prec
  do i = 1, 4
      nid = fixNodes(i)
      K_global(3*nid-2, 3*nid-2) = big_number
      K_global(3*nid-1, 3*nid-1) = big_number
      K_global(3*nid,   3*nid)   = big_number
      rhs(3*nid-2) = 0.0_prec
      rhs(3*nid-1) = 0.0_prec
      rhs(3*nid)   = 0.0_prec
  end do

  ! =========================================================
  ! 3. 求解线性方程组
  ! =========================================================
  call Kcoo%init(n_unknowns)
  do i = 1, n_unknowns
      do j = 1, n_unknowns    
          if (abs(K_global(i,j)) > 1.0e-12_prec) then
              call Kcoo%add_entry(i, j, K_global(i,j))
          end if
      end do
  end do

  Kcsr = COO_to_CSR(Kcoo)
  opts%solver_family = "direct"
  call create_linear_solver(Kcsr, opts, solver)
  call solver%analyze(); call solver%factor(); call solver%solve(rhs, x)

  print *, '================================'
  print *, '   3D Displacement (Node 9)     '
  print *, '================================'
  print *, 'Ux: ', x(3*9 - 2)
  print *, 'Uy: ', x(3*9 - 1)
  print *, 'Uz: ', x(3*9)
  
  ! 导出定制的 3D VTK
  call export_vtk_3d('output_3d_benchmark.vtk', n_nodes, n_elements, nodes, el_cnn, x)
  print *, '3D VTK file generated successfully.'
  
  call solver%free()

contains

  ! =========================================================
  ! 专属 3D VTK 导出子程序 (修正 Intel Fortran 严格语法限制)
  ! =========================================================
  subroutine export_vtk_3d(filename, nnode, nelem, coords, conn, u)
      character(len=*), intent(in) :: filename
      integer, intent(in)          :: nnode, nelem
      ! 修正点：使用显式的传入变量 nnode 和 nelem 来定义数组大小
      real(prec), intent(in)       :: coords(3, nnode)
      integer, intent(in)          :: conn(8, nelem)
      real(prec), intent(in)       :: u(nnode * 3)
      
      integer :: iunit, id

      open(newunit=iunit, file=filename, status='replace')
      write(iunit, '(A)') '# vtk DataFile Version 3.0'
      write(iunit, '(A)') 'FEM 3D Hex8 Results'
      write(iunit, '(A)') 'ASCII'
      write(iunit, '(A)') 'DATASET UNSTRUCTURED_GRID'

      write(iunit, '(A, I8, A)') 'POINTS ', nnode, ' float'
      do id = 1, nnode
          write(iunit, '(3(E14.6, 1X))') coords(1, id), coords(2, id), coords(3, id)
      end do

      write(iunit, '(A, I8, I8)') 'CELLS ', nelem, nelem * 9
      do id = 1, nelem
          write(iunit, '(I2, 8(1X, I8))') 8, (conn(1,id)-1), (conn(2,id)-1), (conn(3,id)-1), &
                                             (conn(4,id)-1), (conn(5,id)-1), (conn(6,id)-1), &
                                             (conn(7,id)-1), (conn(8,id)-1)
      end do

      write(iunit, '(A, I8)') 'CELL_TYPES ', nelem
      do id = 1, nelem
          write(iunit, '(I2)') 12  ! 12 代表 VTK_HEXAHEDRON
      end do

      write(iunit, '(A, I8)') 'POINT_DATA ', nnode
      write(iunit, '(A)') 'SCALARS Node_ID int 1'
      write(iunit, '(A)') 'LOOKUP_TABLE default'
      do id = 1, nnode
          write(iunit, '(I8)') id
      end do

      write(iunit, '(A)') 'VECTORS Displacement float'
      do id = 1, nnode
          write(iunit, '(3(E14.6, 1X))') u(3*id-2), u(3*id-1), u(3*id)
      end do
      
      write(iunit, '(A)') 'SCALARS u_magnitude float 1'
      write(iunit, '(A)') 'LOOKUP_TABLE default'
      do id = 1, nnode
          write(iunit, '(E14.6)') sqrt(u(3*id-2)**2 + u(3*id-1)**2 + u(3*id)**2)
      end do
      close(iunit)
  end subroutine export_vtk_3d

end program test_two_elements_3d