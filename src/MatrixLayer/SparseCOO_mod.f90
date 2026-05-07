module SparseCOO_mod
  use SparseMatrixBase_mod
  use MatrixUtilities_mod, only: merge_sort_2i_r8, unique_merge_2i_r8
  implicit none
  private
  public :: SparseCOO

  integer, parameter :: dp = kind(1.0d0)

  !=========================================================
  !   COO 稀疏格式（动态构建 + 预处理）
  !=========================================================
  type, extends(SparseMatrixBase) :: SparseCOO
  integer :: nnz = 0
  integer :: capacity = 0
  integer, allocatable :: row(:)
  integer, allocatable :: col(:)
  real(dp), allocatable :: val(:)
contains
  procedure :: init        => coo_init
  procedure :: to_dense    => coo_to_dense
  procedure :: spmv        => coo_spmv
  procedure :: get_nnz     => coo_get_nnz

  ! COO 特有接口
  procedure :: add_entry   => coo_add_entry
  procedure :: preprocess  => coo_preprocess

  ! 输出供调试：
  procedure :: print   => print_coo
  end type SparseCOO

contains
!=========================================================
!  COO 初始化
!=========================================================
  subroutine coo_init(self, n)
    class(SparseCOO), intent(inout) :: self
    integer, intent(in) :: n

    self%n        = n
    self%nnz      = 0
    self%capacity = 128

    if (allocated(self%row)) deallocate(self%row)
    if (allocated(self%col)) deallocate(self%col)
    if (allocated(self%val)) deallocate(self%val)

    allocate(self%row(self%capacity))
    allocate(self%col(self%capacity))
    allocate(self%val(self%capacity))
  end subroutine coo_init


!=========================================================
!  自动扩容
!=========================================================
  subroutine coo_grow(self)
    class(SparseCOO), intent(inout) :: self

    integer :: newcap
    integer, allocatable :: r(:), c(:)
    real(dp), allocatable :: v(:)

    newcap = max(128, merge(self%capacity*2, 128, self%capacity>0))

    ! 只拷贝已使用的 nnz 部分
    allocate(r(self%nnz), c(self%nnz), v(self%nnz))
    if (self%nnz > 0) then
      r = self%row(1:self%nnz)
      c = self%col(1:self%nnz)
      v = self%val(1:self%nnz)
    end if

    deallocate(self%row, self%col, self%val)
    allocate(self%row(newcap))
    allocate(self%col(newcap))
    allocate(self%val(newcap))

    if (self%nnz > 0) then
      self%row(1:self%nnz) = r
      self%col(1:self%nnz) = c
      self%val(1:self%nnz) = v
    end if

    self%capacity = newcap
  end subroutine coo_grow


!=========================================================
!  添加非零元素 A(i,j) += value
!=========================================================
  subroutine coo_add_entry(self, i, j, value)
    class(SparseCOO), intent(inout) :: self
    integer, intent(in) :: i, j
    real(dp), intent(in) :: value

    if (self%nnz == self%capacity) call coo_grow(self)

    self%nnz = self%nnz + 1
    self%row(self%nnz) = i
    self%col(self%nnz) = j
    self%val(self%nnz) = value
  end subroutine coo_add_entry


!=========================================================
!  COO → Dense
!=========================================================
  function coo_to_dense(self) result(A)
    class(SparseCOO), intent(in) :: self
    real(dp), allocatable :: A(:,:)
    integer :: k

    allocate(A(self%n, self%n))
    A = 0.0_dp

    do k = 1, self%nnz
      A(self%row(k), self%col(k)) = A(self%row(k), self%col(k)) + self%val(k)
    end do
  end function coo_to_dense


!=========================================================
!  COO SpMV: y = A*x
!=========================================================
  subroutine coo_spmv(self, x, y)
    class(SparseCOO), intent(in) :: self
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    integer :: k

    y = 0.0_dp
    do k = 1, self%nnz
      y(self%row(k)) = y(self%row(k)) + self%val(k) * x(self%col(k))
    end do
  end subroutine coo_spmv


!=========================================================
!  获取 nnz
!=========================================================
  function coo_get_nnz(self) result(nz)
    class(SparseCOO), intent(in) :: self
    integer :: nz
    nz = self%nnz
  end function coo_get_nnz


!=========================================================
!  COO 预处理：排序 + 去重
!  排序规则：按 (row, col) 升序
!  去重规则：相同 (row,col) 的 val 累加
!=========================================================
  subroutine coo_preprocess(self)
    class(SparseCOO), intent(inout) :: self

    if (self%nnz <= 1) return

    ! 只对前 nnz 个有效元素排序
    call merge_sort_2i_r8(self%row, self%col, self%val, 1, self%nnz)

    ! 合并重复 (row,col)
    call unique_merge_2i_r8(self%row, self%col, self%val, self%nnz)

    ! 注意：capacity 不变，nnz 可能变小，前 nnz 部分有效
  end subroutine coo_preprocess

! 输出 供调试：
  subroutine print_coo(self)
    use iso_fortran_env, only: wp => real64
    implicit none
    class(SparseCOO), intent(in) :: self
    integer :: k

    print *, "================ COO (row, col, val) ================"
    print *, "nnz =", self%nnz

    if (self%nnz == 0) then
      print *, "<empty>"
      return
    end if

    print '(A4,3A12)', " idx", "row", "col", "val"
    do k = 1, self%nnz
      print '(I4,2I12,F12.4)', k, self%row(k), self%col(k), self%val(k)
    end do

    print *, "====================================================="
  end subroutine print_coo
end module SparseCOO_mod
