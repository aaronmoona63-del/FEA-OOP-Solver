module SparseConvert_mod
  use SparseMatrixBase_mod
  use SparseCOO_mod
  use SparseCSR_mod
  use SparseSkyline_mod
  use MatrixUtilities_mod, only: merge_sort_i_r8
  implicit none
  private

  public :: COO_to_CSR
  public :: COO_to_Skyline
  public :: CSR_to_COO
  public :: Skyline_to_COO

!
!   CSR 存储格式：
!       row_ptr(i)            —— 第 i 行的起始位置（1-based）
!       col_ind(k)            —— 第 k 个非零元素的列号
!       val(k)                —— 第 k 个非零元素的数值
!
!   COO 存储格式：
!       row(p), col(p), val(p)
!       每个 p 对应一个非零元素的 (行, 列, 值)

contains

!====================================================================
!  COO → CSR
!====================================================================
  function COO_to_CSR(mtxcoo) result(mtxcsr)
    class(SparseCOO), intent(in) :: mtxcoo ! COO 按行且行内按列有序, 没有重复的 (row, col)
    type(SparseCSR) :: mtxcsr

    ! local variables
    integer :: n, nnz
    integer :: i, k, r, p
    integer, allocatable :: row_count(:), next_ptr(:)

    n   = mtxcoo%n
    nnz = mtxcoo%nnz

    call mtxcsr%init(n)
    mtxcsr%nnz = nnz

    allocate(row_count(n))
    row_count = 0

    ! Step 1: row_count
    do k = 1, nnz
      row_count( mtxcoo%row(k) ) = row_count( mtxcoo%row(k) ) + 1
    end do

    ! Step 2: prefix sum -> row_ptr
    allocate(mtxcsr%row_ptr(n+1))
    mtxcsr%row_ptr(1) = 1

    do i = 1, n
      mtxcsr%row_ptr(i+1) = mtxcsr%row_ptr(i) + row_count(i)
    end do

    allocate(mtxcsr%col_ind(nnz))
    allocate(mtxcsr%val(nnz))

    allocate(next_ptr(n))
    next_ptr = mtxcsr%row_ptr(1:n)

    ! Step 3: fill CSR
    do k = 1, nnz
      r = mtxcoo%row(k)
      p = next_ptr(r)

      mtxcsr%col_ind(p) = mtxcoo%col(k)
      mtxcsr%val(p)     = mtxcoo%val(k)

      next_ptr(r) = p + 1
    end do
  end function COO_to_CSR

!====================================================================
!  最终版：COO → Skyline（不使用指针）
!
!  输入：
!      mtxcoo    = COO 格式矩阵
!      symmetric = 是否按实对称方式存储（只存 AL）
!      force_sort = 若 COO 未排序是否自动排序
!
!  输出：
!      mtxsky    = Skyline 格式矩阵
!
!====================================================================
!   COO → Skyline
!====================================================================
  function COO_to_Skyline(mtxcoo, symmetric) result(mtxsky)
    class(SparseCOO), intent(in) :: mtxcoo
    logical, intent(in) :: symmetric       ! false = unsym, true = symmetric storage

    type(SparseSkyline) :: mtxsky
    integer :: n, nnz
    integer :: i, j, k
    integer, allocatable :: first_row(:), first_col(:)
    integer :: is, height
    integer :: lastend, thisend, pos
    integer :: jh, is_j

    real(8) :: v

    n   = mtxcoo%n
    nnz = mtxcoo%nnz

    call mtxsky%init(n)

    !------------------------------------------------
    ! Step 1：构建 skyline profile：first_row(j) & first_col(j)
    !------------------------------------------------
    allocate(first_row(n))
    allocate(first_col(n))

    first_row = [(j, j=1,n)]   ! 列 j 向上扫描的最小行号（默认对角线）
    first_col = [(j, j=1,n)]   ! 行 j 向左扫描的最小列号（默认对角线）

    do k = 1, nnz
      i = mtxcoo%row(k)
      j = mtxcoo%col(k)

      if (i < j) then
        ! A(i,j) 属于上三角 → 列扫描用来更新 first_row(j)
        if (i < first_row(j)) first_row(j) = i
      else if (i > j) then
        ! A(i,j) 属于下三角 → 行扫描用来更新 first_col(i)
        if (j < first_col(i)) first_col(i) = j
      end if
    end do

    !------------------------------------------------
    ! Step 2：综合行列得到准确profile,
    !                 min(first_row(j), first_col(j))
    !        并构建 jp（列指针）
    !------------------------------------------------
    allocate(mtxsky%jp(n))
    mtxsky%jp = 0

    do j = 1, n
      is = min(first_row(j), first_col(j))
      height = j - is
      if (j == 1) then
        mtxsky%jp(j) = height
      else
        mtxsky%jp(j) = mtxsky%jp(j-1) + height
      end if
    end do

    mtxsky%nnz = mtxsky%jp(n) + n    ! add diagonal count

    ! allocate storage
    allocate(mtxsky%ad(n))
    allocate(mtxsky%al(mtxsky%jp(n)))
    allocate(mtxsky%au(mtxsky%jp(n)))

    mtxsky%ad = 0.0d0
    mtxsky%al = 0.0d0
    mtxsky%au = 0.0d0

    !------------------------------------------------
    ! Step 3：扫描 COO 填 AD + AL + AU（按 Skyline 结构）
    !------------------------------------------------
    do k = 1, nnz
      i = mtxcoo%row(k)
      j = mtxcoo%col(k)
      v = mtxcoo%val(k)

      !===============================
      ! 1. 对角线：只进这里一次
      !===============================
      if (i == j) then
        mtxsky%ad(i) = mtxsky%ad(i) + v
        cycle         ! 关键：直接跳到下一个 k，后面分支不再走
      end if
      !===================================================
      ! 3. 此时已经保证 i ≠ j（对角线前面已经 cycle 掉了）
      !    非对称情形：上三角 → AU，下三角 → AL
      !    对称情况： 下三角COO直接舍弃
      !===================================================
      if (i < j) then
        !------------------- 上三角，列 j 填 AU -------------------
        if (j == 1) then
          lastend = 0
        else
          lastend = mtxsky%jp(j-1)
        end if
        thisend = mtxsky%jp(j)
        jh     = thisend - lastend
        is_j   = j - jh

        if (i < is_j ) cycle  ! 不在 skyline
        pos = lastend + (i - is_j + 1)

        mtxsky%au(pos) = mtxsky%au(pos) + v

      else if (i > j) then
        if (symmetric) cycle ! 后面舍弃
        !------------------- 下三角，列 i 填 AL -------------------
        if (i == 1) then
          lastend = 0
        else
          lastend = mtxsky%jp(i-1)
        end if
        thisend = mtxsky%jp(i)
        jh     = thisend - lastend
        is_j   = i - jh

        if (j < is_j ) cycle ! 不在 skyline
        pos = lastend + (j - is_j + 1)

        mtxsky%al(pos) = mtxsky%al(pos) + v
      end if

    end do

    if (symmetric) then
      mtxsky%al = mtxsky%au
    end if

  end function COO_to_Skyline

!====================================================================
!   Skyline → COO（支持 symmetric / unsymmetric）
!====================================================================
  function Skyline_to_COO(mtxsky, symmetric) result(mtxcoo)
    class(SparseSkyline), intent(in) :: mtxsky
    logical, intent(in) :: symmetric      ! .true. : only AL is meaningful

    type(SparseCOO) :: mtxcoo
    integer :: n, j, iidx, base, jh, is, row

    n = mtxsky%n
    call mtxcoo%init(n)

    !===========================================================
    ! 1. 对角线
    !===========================================================
    do j = 1, n
      if (mtxsky%ad(j) /= 0.0d0) call mtxcoo%add_entry(j, j, mtxsky%ad(j))
    end do

    !===========================================================
    ! 2. 处理非对角（AL / AU）
    !
    ! Skyline 列 j 结构：
    !   base = jp(j-1)
    !   jh   = jp(j) - base
    !   is   = j - jh  (最顶端行号)
    !
    ! AL(base + k) 对应    (j,   is + k - 1)
    ! AU(base + k) 对应    (is + k - 1, j)
    !===========================================================
    do j = 1, n

      if (j == 1) then
        base = 0
      else
        base = mtxsky%jp(j-1)
      end if

      jh = mtxsky%jp(j) - base     ! 此列 Skyline 高度

      if (jh <= 0) cycle

      is = j - jh

      do iidx = 1, jh
        row = is + (iidx - 1)
        !==========================
        ! 上三角 (row, j)，来自 AU
        !==========================
        if (mtxsky%au(base + iidx) /= 0.0d0) then
          call mtxcoo%add_entry(row, j, mtxsky%au(base + iidx))
        end if
        !==========================
        ! 下三角 (j, row)，来自 AL
        ! symmetric = .true. 时，AL 不存在，也不输出
        !==========================
        if (.not. symmetric) then
          if (mtxsky%al(base + iidx) /= 0.0d0) then
            call mtxcoo%add_entry(j, row, mtxsky%al(base + iidx))
          end if
        end if
      end do
    end do
  end function Skyline_to_COO


!====================================================================
!   CSR → COO
!
!   转换原理：
!       对于 CSR 的每一行 i，
!       非零范围是 k = row_ptr(i) : row_ptr(i+1)-1
!       将这些元素展开即可。
!====================================================================
  function CSR_to_COO(mtxcsr) result(mtxcoo)
    class(SparseCSR), intent(in) :: mtxcsr
    type(SparseCOO) :: mtxcoo

    integer :: n, nnz
    integer :: i, k, p

    !------------------------------------------------------------
    ! 基本信息
    !------------------------------------------------------------
    n   = mtxcsr%n
    nnz = mtxcsr%nnz

    call mtxcoo%init(n)

    ! 分配 COO 三数组
    if (allocated(mtxcoo%row)) deallocate(mtxcoo%row)
    if (allocated(mtxcoo%col)) deallocate(mtxcoo%col)
    if (allocated(mtxcoo%val)) deallocate(mtxcoo%val)
    allocate(mtxcoo%row(nnz))
    allocate(mtxcoo%col(nnz))
    allocate(mtxcoo%val(nnz))

    ! p 指向 COO 写入位置（1-based）
    p = 1

    !------------------------------------------------------------
    ! 核心步骤：展开 CSR 每行的非零
    ! row_ptr(i) ~ row_ptr(i+1)-1 是第 i 行的所有非零
    !------------------------------------------------------------
    do i = 1, n
      do k = mtxcsr%row_ptr(i), mtxcsr%row_ptr(i+1) - 1

        ! COO 的行号（显式）
        mtxcoo%row(p) = i

        ! COO 的列号从 CSR 直接取
        mtxcoo%col(p) = mtxcsr%col_ind(k)

        ! COO 的数值
        mtxcoo%val(p) = mtxcsr%val(k)

        p = p + 1
      end do
    end do

    ! 最终非零数
    mtxcoo%nnz = nnz

  end function CSR_to_COO

end module SparseConvert_mod
