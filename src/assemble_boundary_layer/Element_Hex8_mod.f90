module Element_Hex8_mod
  use Types, only : prec
  use ElementBase_mod, only : ElementBase
  implicit none
  private
  public :: Element_Hex8

  ! 继承基类
  type, extends(ElementBase) :: Element_Hex8
  contains
      ! 实现合同里规定的 get_stiffness_matrix 方法
      procedure :: get_stiffness_matrix => hex8_get_stiffness_matrix
  end type Element_Hex8

contains

  subroutine hex8_get_stiffness_matrix(this, coords, props, Ke)
      class(Element_Hex8), intent(in) :: this
      real(prec), intent(in)          :: coords(:,:)  ! Hex8 应该是 (3, 8) 的大小
      real(prec), intent(in)          :: props(:)     ! 我们约定 props(1)是E，props(2)是nu
      real(prec), intent(inout)       :: Ke(:,:)      ! 返回 (24, 24) 的矩阵

      ! ====== 这里面全是你之前写过的核心物理逻辑，原封不动搬过来 ======
      real(prec) :: E, nu, c1, detJ, weight
      real(prec) :: D(6, 6), B(6, 24), dNdxi(3, 8), dNdx(3, 8)
      real(prec) :: Jac(3, 3), invJac(3, 3), temp_DB(6, 24)
      real(prec) :: xi, eta, zeta
      real(prec) :: gp_loc(2)
      integer :: i, rr, cc, gp1, gp2, gp3

      ! 1. 解析材料属性
      E  = props(1)
      nu = props(2)
      Ke = 0.0_prec

      ! 2. 计算 3D 弹性本构矩阵 D (6x6)
      c1 = E / ((1.0_prec + nu) * (1.0_prec - 2.0_prec * nu))
      D = 0.0_prec
      D(1,1) = c1 * (1.0_prec - nu); D(1,2) = c1 * nu;               D(1,3) = c1 * nu
      D(2,1) = c1 * nu;              D(2,2) = c1 * (1.0_prec - nu);  D(2,3) = c1 * nu
      D(3,1) = c1 * nu;              D(3,2) = c1 * nu;               D(3,3) = c1 * (1.0_prec - nu)
      D(4,4) = c1 * (1.0_prec - 2.0_prec * nu) / 2.0_prec
      D(5,5) = D(4,4); D(6,6) = D(4,4)

      gp_loc = [-1.0_prec/sqrt(3.0_prec), 1.0_prec/sqrt(3.0_prec)]

      ! 3. 2x2x2 高斯积分大循环
      do gp1 = 1, 2
          do gp2 = 1, 2
              do gp3 = 1, 2
                  xi = gp_loc(gp1); eta = gp_loc(gp2); zeta = gp_loc(gp3)
                  weight = 1.0_prec * 1.0_prec * 1.0_prec
                  
                  dNdxi(1,:) = [ -(1.0_prec-eta)*(1.0_prec-zeta)/8.0_prec,  (1.0_prec-eta)*(1.0_prec-zeta)/8.0_prec,  (1.0_prec+eta)*(1.0_prec-zeta)/8.0_prec, -(1.0_prec+eta)*(1.0_prec-zeta)/8.0_prec, &
                                 -(1.0_prec-eta)*(1.0_prec+zeta)/8.0_prec,  (1.0_prec-eta)*(1.0_prec+zeta)/8.0_prec,  (1.0_prec+eta)*(1.0_prec+zeta)/8.0_prec, -(1.0_prec+eta)*(1.0_prec+zeta)/8.0_prec ]
                  dNdxi(2,:) = [ -(1.0_prec-xi)*(1.0_prec-zeta)/8.0_prec,  -(1.0_prec+xi)*(1.0_prec-zeta)/8.0_prec,   (1.0_prec+xi)*(1.0_prec-zeta)/8.0_prec,   (1.0_prec-xi)*(1.0_prec-zeta)/8.0_prec,  &
                                 -(1.0_prec-xi)*(1.0_prec+zeta)/8.0_prec,  -(1.0_prec+xi)*(1.0_prec+zeta)/8.0_prec,   (1.0_prec+xi)*(1.0_prec+zeta)/8.0_prec,   (1.0_prec-xi)*(1.0_prec+zeta)/8.0_prec ]
                  dNdxi(3,:) = [ -(1.0_prec-xi)*(1.0_prec-eta)/8.0_prec,   -(1.0_prec+xi)*(1.0_prec-eta)/8.0_prec,   -(1.0_prec+xi)*(1.0_prec+eta)/8.0_prec,   -(1.0_prec-xi)*(1.0_prec+eta)/8.0_prec,   &
                                  (1.0_prec-xi)*(1.0_prec-eta)/8.0_prec,    (1.0_prec+xi)*(1.0_prec-eta)/8.0_prec,    (1.0_prec+xi)*(1.0_prec+eta)/8.0_prec,    (1.0_prec-xi)*(1.0_prec+eta)/8.0_prec ]
                  
                  Jac = 0.0_prec
                  do i = 1, 8
                      do rr = 1, 3
                          do cc = 1, 3
                              Jac(rr,cc) = Jac(rr,cc) + dNdxi(rr,i) * coords(cc, i)  ! 注意这里改为了直接取 coords 数组
                          end do
                      end do
                  end do
                  
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
      
  end subroutine hex8_get_stiffness_matrix

end module Element_Hex8_mod