module PostProcessor_mod
  use Types, only : prec
  implicit none
  private
  public :: calc_hex8_von_mises

contains

  ! 极其安全的独立方法：不修改任何底层架构，直接传入纯数组进行后处理
  subroutine calc_hex8_von_mises(coords, props, u_elem, von_mises)
      real(prec), intent(in)  :: coords(3, 8)
      real(prec), intent(in)  :: props(2)
      real(prec), intent(in)  :: u_elem(24)
      real(prec), intent(out) :: von_mises

      real(prec) :: E, nu, c1, detJ
      real(prec) :: D(6, 6), B(6, 24), dNdxi(3, 8), dNdx(3, 8)
      real(prec) :: Jac(3, 3), invJac(3, 3)
      real(prec) :: stress(6)
      real(prec) :: xi, eta, zeta
      integer :: i, rr, cc

      E = props(1); nu = props(2)
      c1 = E / ((1.0_prec + nu) * (1.0_prec - 2.0_prec * nu))
      D = 0.0_prec
      D(1,1) = c1*(1.0_prec-nu); D(1,2) = c1*nu;             D(1,3) = c1*nu
      D(2,1) = c1*nu;            D(2,2) = c1*(1.0_prec-nu);  D(2,3) = c1*nu
      D(3,1) = c1*nu;            D(3,2) = c1*nu;             D(3,3) = c1*(1.0_prec-nu)
      D(4,4) = c1*(1.0_prec-2.0_prec*nu)/2.0_prec
      D(5,5) = D(4,4); D(6,6) = D(4,4)

      ! 取几何中心点计算应力 (xi=0, eta=0, zeta=0)
      xi = 0.0_prec; eta = 0.0_prec; zeta = 0.0_prec
      
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
                  Jac(rr,cc) = Jac(rr,cc) + dNdxi(rr,i) * coords(cc, i)
              end do
          end do
      end do
      
      detJ = Jac(1,1)*(Jac(2,2)*Jac(3,3) - Jac(2,3)*Jac(3,2)) - Jac(1,2)*(Jac(2,1)*Jac(3,3) - Jac(2,3)*Jac(3,1)) + Jac(1,3)*(Jac(2,1)*Jac(3,2) - Jac(2,2)*Jac(3,1))
           
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
      
      ! 应力矩阵相乘
      stress = matmul(D, matmul(B, u_elem))
      
      ! 冯·米塞斯等效应力
      von_mises = sqrt( 0.5_prec * ( (stress(1)-stress(2))**2 + (stress(2)-stress(3))**2 + (stress(3)-stress(1))**2 &
                               + 6.0_prec * (stress(4)**2 + stress(5)**2 + stress(6)**2) ) )
  end subroutine calc_hex8_von_mises

end module PostProcessor_mod