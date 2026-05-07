program test_four_elements_matlab
  use Types, only : prec
  use SparseCOO_mod, only : SparseCOO
  use SparseCSR_mod, only : SparseCSR
  use SparseConvert_mod, only : COO_to_CSR
  use LinearSolverBase_mod, only : LinearSolver, LinearSolverOptions
  use LinearSolverFactory_mod, only : create_linear_solver
  use VTK_Export_mod, only : export_to_vtk
  implicit none

  type(SparseCOO) :: Kcoo
  type(SparseCSR), target :: Kcsr
  class(LinearSolver), allocatable :: solver
  type(LinearSolverOptions) :: opts

  real(prec), allocatable :: rhs(:), x(:)
  real(prec) :: nodes(2, 9), E, nu, thickness, detJ, weight
  real(prec) :: D(3, 3), Ke(8, 8), B(3, 8), dNdxi(2, 4), dNdx(2, 4)
  real(prec) :: Jac(2, 2), invJac(2, 2)
  real(prec) :: K_global(18, 18)
  real(prec) :: xi, eta, big_number
  real(prec) :: gp_loc(2)
  real(prec) :: mag7, mag8, mag9  
  
  ! 应力后处理变量
  real(prec) :: de(8)        
  real(prec) :: temp_stress(3)
  real(prec) :: sigma(3)     
  
  integer :: el_cnn(4, 4), conn(4, 4)
  integer :: n_nodes, n_elements, n_unknowns, i, j, el, r, c, gp1, gp2, nid
  
  print *, '--- [MATLAB Benchmark: Displacements & Stress] ---'
  
  n_nodes = 9
  n_elements = 4
  n_unknowns = n_nodes * 2
  
  nodes(1,:) = [0.0_prec, 0.0_prec, 0.0_prec, 1000.0_prec, 1000.0_prec, 1000.0_prec, 2000.0_prec, 2000.0_prec, 2000.0_prec]
  nodes(2,:) = [2000.0_prec, 1000.0_prec, 0.0_prec, 2000.0_prec, 1000.0_prec, 0.0_prec, 2000.0_prec, 1000.0_prec, 0.0_prec]
  
  el_cnn(:,1) = [1, 2, 5, 4]
  el_cnn(:,2) = [2, 3, 6, 5]
  el_cnn(:,3) = [5, 6, 9, 8]
  el_cnn(:,4) = [4, 5, 8, 7]
  
  ! 修正为真实的物理参数 (钢材弹性模量)
  E = 210000.0_prec
  nu = 0.3_prec
  thickness = 1.0_prec
  
  D = 0.0_prec
  D(1,1) = 1.0_prec
  D(1,2) = nu
  D(2,1) = nu
  D(2,2) = 1.0_prec
  D(3,3) = (1.0_prec - nu) / 2.0_prec
  D = D * (E / (1.0_prec - nu**2))

  K_global = 0.0_prec
  allocate(rhs(n_unknowns), x(n_unknowns))
  rhs = 0.0_prec
  x = 0.0_prec
  
  gp_loc = [-1.0_prec/sqrt(3.0_prec), 1.0_prec/sqrt(3.0_prec)]
  
  ! 1. 组装全局刚度矩阵
  do el = 1, n_elements
      Ke = 0.0_prec
      do gp1 = 1, 2
          do gp2 = 1, 2
              xi = gp_loc(gp1)
              eta = gp_loc(gp2)
              weight = 1.0_prec * 1.0_prec
              
              dNdxi(1,:) = [-(1.0_prec-eta)/4.0_prec,  (1.0_prec-eta)/4.0_prec, (1.0_prec+eta)/4.0_prec, -(1.0_prec+eta)/4.0_prec]
              dNdxi(2,:) = [-(1.0_prec-xi)/4.0_prec,  -(1.0_prec+xi)/4.0_prec,  (1.0_prec+xi)/4.0_prec,   (1.0_prec-xi)/4.0_prec]
              
              Jac = 0.0_prec
              do i = 1, 4
                  Jac(1,1) = Jac(1,1) + dNdxi(1,i) * nodes(1, el_cnn(i,el))
                  Jac(1,2) = Jac(1,2) + dNdxi(1,i) * nodes(2, el_cnn(i,el))
                  Jac(2,1) = Jac(2,1) + dNdxi(2,i) * nodes(1, el_cnn(i,el))
                  Jac(2,2) = Jac(2,2) + dNdxi(2,i) * nodes(2, el_cnn(i,el))
              end do
              
              detJ = Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1)
              invJac(1,1) =  Jac(2,2)/detJ
              invJac(1,2) = -Jac(1,2)/detJ
              invJac(2,1) = -Jac(2,1)/detJ
              invJac(2,2) =  Jac(1,1)/detJ
              dNdx = matmul(invJac, dNdxi)
              
              B = 0.0_prec
              do i = 1, 4
                  B(1, 2*i-1) = dNdx(1,i)
                  B(2, 2*i)   = dNdx(2,i)
                  B(3, 2*i-1) = dNdx(2,i)
                  B(3, 2*i)   = dNdx(1,i)
              end do
              Ke = Ke + matmul(transpose(B), matmul(D, B)) * detJ * weight * thickness
          end do
      end do
      
      do i = 1, 4
          do j = 1, 4
              r = 2*(el_cnn(i,el)-1)
              c = 2*(el_cnn(j,el)-1)
              K_global(r+1, c+1) = K_global(r+1, c+1) + Ke(2*i-1, 2*j-1)
              K_global(r+1, c+2) = K_global(r+1, c+2) + Ke(2*i-1, 2*j)
              K_global(r+2, c+1) = K_global(r+2, c+1) + Ke(2*i,   2*j-1)
              K_global(r+2, c+2) = K_global(r+2, c+2) + Ke(2*i,   2*j)
          end do
      end do
  end do

  ! 2. 施加载荷与边界条件
  rhs(13) =  50000.0_prec
  rhs(15) = 100000.0_prec
  rhs(17) =  50000.0_prec
  big_number = 1.0e15_prec
  
  do i = 1, 6
      K_global(i, i) = big_number
      rhs(i) = 0.0_prec
  end do

  ! 3. 求解线性方程组
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
  
  call solver%analyze()
  call solver%factor()
  call solver%solve(rhs, x)

  mag7 = sqrt(x(13)**2 + x(14)**2)
  mag8 = sqrt(x(15)**2 + x(16)**2)
  mag9 = sqrt(x(17)**2 + x(18)**2)

  print *, '================================'
  print *, '   Displacement Magnitude (U)   '
  print *, '================================'
  print *, 'Node 7 U:    ', mag7
  print *, 'Node 8 U:    ', mag8
  print *, 'Node 9 U:    ', mag9

  ! =========================================================
  ! 新增模块：应力恢复 (Stress Recovery) 稳健版
  ! =========================================================
  print *, '================================'
  print *, '   Element Stresses (Center)    '
  print *, '================================'
  do el = 1, n_elements
      de = 0.0_prec
      do i = 1, 4
          nid = el_cnn(i, el)
          de(2*i - 1) = x(2*nid - 1)  ! u
          de(2*i)     = x(2*nid)      ! v
      end do
      
      xi = 0.0_prec
      eta = 0.0_prec
      
      dNdxi(1,:) = [-(1.0_prec-eta)/4.0_prec,  (1.0_prec-eta)/4.0_prec, (1.0_prec+eta)/4.0_prec, -(1.0_prec+eta)/4.0_prec]
      dNdxi(2,:) = [-(1.0_prec-xi)/4.0_prec,  -(1.0_prec+xi)/4.0_prec,  (1.0_prec+xi)/4.0_prec,   (1.0_prec-xi)/4.0_prec]
      
      Jac = 0.0_prec
      do i = 1, 4
          Jac(1,1) = Jac(1,1) + dNdxi(1,i) * nodes(1, el_cnn(i,el))
          Jac(1,2) = Jac(1,2) + dNdxi(1,i) * nodes(2, el_cnn(i,el))
          Jac(2,1) = Jac(2,1) + dNdxi(2,i) * nodes(1, el_cnn(i,el))
          Jac(2,2) = Jac(2,2) + dNdxi(2,i) * nodes(2, el_cnn(i,el))
      end do
      
      detJ = Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1)
      invJac(1,1) =  Jac(2,2)/detJ
      invJac(1,2) = -Jac(1,2)/detJ
      invJac(2,1) = -Jac(2,1)/detJ
      invJac(2,2) =  Jac(1,1)/detJ
      dNdx = matmul(invJac, dNdxi)
      
      B = 0.0_prec
      do i = 1, 4
          B(1, 2*i-1) = dNdx(1,i)
          B(2, 2*i)   = dNdx(2,i)
          B(3, 2*i-1) = dNdx(2,i)
          B(3, 2*i)   = dNdx(1,i)
      end do
      
      ! 拆分矩阵乘法，防止部分编译器报错
      temp_stress = matmul(B, de)
      sigma = matmul(D, temp_stress)
      
      ! 使用最通用的打印方式
      print *, 'Element ', el, ' | Sig_x: ', sigma(1), ' | Sig_y: ', sigma(2), ' | Tau_xy: ', sigma(3)
  end do
  
  ! 导出 VTK
  conn(:,1) = el_cnn(:,1)
  conn(:,2) = el_cnn(:,2)
  conn(:,3) = el_cnn(:,3)
  conn(:,4) = el_cnn(:,4)
  call export_to_vtk('matlab_benchmark_results.vtk', n_nodes, n_elements, nodes, conn, x)
  print *, 'VTK file [matlab_benchmark_results.vtk] generated successfully.'
  
  call solver%free()
end program test_four_elements_matlab