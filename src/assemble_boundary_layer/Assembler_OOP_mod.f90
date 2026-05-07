module Assembler_OOP_mod
  use Types, only : prec
  use SparseCOO_mod, only : SparseCOO
  use Element_Hex8_mod, only : Element_Hex8
  implicit none
  private
  public :: assemble_global_stiffness

contains

  ! 纯粹的数据加工厂：接收 6 个参数，不依赖任何全局变量
  subroutine assemble_global_stiffness(Kcoo, n_elements, n_nodes, coords, conn, props)
      type(SparseCOO), intent(inout) :: Kcoo
      integer, intent(in)            :: n_elements, n_nodes
      real(prec), intent(in)         :: coords(3, n_nodes)
      integer, intent(in)            :: conn(8, n_elements)
      real(prec), intent(in)         :: props(2)

      integer :: el, i, j
      real(prec) :: coords_elem(3, 8), Ke(24, 24)
      type(Element_Hex8) :: hex8_elem
      integer :: gdofs(24)

      print *, "--- Starting Pure OOP Global Assembly ---"
      do el = 1, n_elements
          ! 提取局部坐标和全局映射
          do i = 1, 8
              coords_elem(:, i) = coords(:, conn(i, el))
              gdofs(3*i-2) = 3 * conn(i, el) - 2
              gdofs(3*i-1) = 3 * conn(i, el) - 1
              gdofs(3*i)   = 3 * conn(i, el)
          end do
          
          ! 呼叫 Hex8 计算物理属性
          call hex8_elem%get_stiffness_matrix(coords_elem, props, Ke)
          
          ! 组装到全局
          do i = 1, 24
              do j = 1, 24
                  if (abs(Ke(i, j)) > 1.0e-12_prec) then
                      call Kcoo%add_entry(gdofs(i), gdofs(j), Ke(i, j))
                  end if
              end do
          end do
      end do
      print *, "--- Assembly Completed Successfully ---"
  end subroutine assemble_global_stiffness
end module Assembler_OOP_mod