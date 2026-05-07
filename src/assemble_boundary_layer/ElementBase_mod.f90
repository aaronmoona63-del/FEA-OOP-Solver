module ElementBase_mod
  use Types, only : prec
  implicit none
  private
  public :: ElementBase

  ! 定义抽象基类
  type, abstract :: ElementBase
  contains
      ! 延迟绑定（纯虚函数），强制所有的子类都必须实现这个计算 Ke 的接口
      procedure(calc_Ke_interface), deferred :: get_stiffness_matrix
  end type ElementBase

  ! 接口的“合同”长这样
  abstract interface
      subroutine calc_Ke_interface(this, coords, props, Ke)
          import :: ElementBase, prec
          class(ElementBase), intent(in) :: this
          ! coords: 当前单元的节点坐标，维度是 (ndims, nnode)
          real(prec), intent(in)         :: coords(:,:) 
          ! props: 当前单元的材料和几何属性，比如 [E, nu, thickness]
          real(prec), intent(in)         :: props(:)    
          ! Ke: 返回组装好的局部刚度矩阵
          real(prec), intent(inout)      :: Ke(:,:)     
      end subroutine calc_Ke_interface
  end interface

end module ElementBase_mod