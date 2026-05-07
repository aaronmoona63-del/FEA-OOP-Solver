module BoundaryCondition_mod
  use Types, only : prec
  use SparseCSR_mod, only : SparseCSR
  implicit none
  private
  public :: apply_benchmark_loads, apply_benchmark_bcs_csr

contains

  ! =========================================================
  ! 1. 独立施加载荷到 RHS 向量 (全网格自适应版)
  ! =========================================================
  subroutine apply_benchmark_loads(n_nodes, coords, rhs)
      integer, intent(in)       :: n_nodes
      real(prec), intent(in)    :: coords(3, n_nodes)
      real(prec), intent(inout) :: rhs(n_nodes * 3)
      integer :: i
      real(prec) :: w_x, w_y, total_weight

      rhs = 0.0_prec
      total_weight = 0.0_prec

      ! 第一遍：扫描顶面，计算几何总权重 (基于形函数面积分配原理)
      do i = 1, n_nodes
          if (abs(coords(3, i) - 50.0_prec) < 1.0e-5_prec) then
              w_x = 1.0_prec; w_y = 1.0_prec
              ! 处于 X 边界的节点，管辖面积减半
              if (abs(coords(1,i) - 0.0_prec) < 1.0e-5_prec .or. abs(coords(1,i) - 10.0_prec) < 1.0e-5_prec) w_x = 0.5_prec
              ! 处于 Y 边界的节点，管辖面积减半
              if (abs(coords(2,i) - 0.0_prec) < 1.0e-5_prec .or. abs(coords(2,i) - 10.0_prec) < 1.0e-5_prec) w_y = 0.5_prec
              total_weight = total_weight + (w_x * w_y)
          end if
      end do

      ! 第二遍：将 10,000 N 的总拉力，严格按权重分配给所有顶面节点
      do i = 1, n_nodes
          if (abs(coords(3, i) - 50.0_prec) < 1.0e-5_prec) then
              w_x = 1.0_prec; w_y = 1.0_prec
              if (abs(coords(1,i) - 0.0_prec) < 1.0e-5_prec .or. abs(coords(1,i) - 10.0_prec) < 1.0e-5_prec) w_x = 0.5_prec
              if (abs(coords(2,i) - 0.0_prec) < 1.0e-5_prec .or. abs(coords(2,i) - 10.0_prec) < 1.0e-5_prec) w_y = 0.5_prec
              
              ! 无论网格多密，总拉力永远是 10000 N
              rhs(3*i) = 10000.0_prec * (w_x * w_y) / total_weight
          end if
      end do
  end subroutine apply_benchmark_loads
  ! =========================================================
  ! 2. 在 CSR 矩阵上施加完美对称的位移约束 (置零置一法)
  ! =========================================================
  subroutine apply_benchmark_bcs_csr(n_nodes, coords, Kcsr, rhs)
      integer, intent(in)            :: n_nodes
      real(prec), intent(in)         :: coords(3, n_nodes)
      type(SparseCSR), intent(inout) :: Kcsr
      real(prec), intent(inout)      :: rhs(n_nodes * 3)

      logical, allocatable :: is_fixed(:)
      integer :: i, j, k, n_unknowns

      n_unknowns = n_nodes * 3
      allocate(is_fixed(n_unknowns))
      is_fixed = .false.

      ! 1. 扫描找出所有被固定的自由度 (Z=0底面)
      do i = 1, n_nodes
          if (abs(coords(3, i) - 0.0_prec) < 1.0e-5_prec) then
              is_fixed(3*i - 2) = .true.
              is_fixed(3*i - 1) = .true.
              is_fixed(3*i)     = .true.
          end if
      end do

      ! 2. 强行将固定自由度对应的载荷清零
      do i = 1, n_unknowns
          if (is_fixed(i)) then
              rhs(i) = 0.0_prec
          end if
      end do

      ! 3. 遍历 CSR 修改矩阵元素 (行、列全置 0，对角线置 1)
      ! 这一步堪称算法艺术，O(NNZ) 复杂度完美搞定对称置零
      do i = 1, n_unknowns
          do k = Kcsr%row_ptr(i), Kcsr%row_ptr(i+1) - 1
              j = Kcsr%col_ind(k)
              ! 如果当前行被固定，或者当前列被固定
              if (is_fixed(i) .or. is_fixed(j)) then
                  if (i == j) then
                      Kcsr%val(k) = 1.0_prec  ! 对角线置1
                  else
                      Kcsr%val(k) = 0.0_prec  ! 非对角线彻底清零
                  end if
              end if
          end do
      end do

      deallocate(is_fixed)
      print *, "=> 独立模块: CSR 对称约束(置零置一法) 已成功施加！"
  end subroutine apply_benchmark_bcs_csr

end module BoundaryCondition_mod