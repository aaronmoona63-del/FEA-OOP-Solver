module SparseSkyline_mod
    use SparseMatrixBase_mod
    implicit none
    private
    public :: SparseSkyline

    !============================================================
    ! Skyline 稀疏矩阵格式（只负责存储和基础操作）
    !
    ! 结构:
    !   ad(j)                 对角 j
    !   al(k)                 下三角（行向）
    !   au(k)                 上三角（列向）
    !   jp(j)                 第 j 列的累积下三角存储偏移
    !
    ! 注意：
    !   本模块不负责构造 Skyline（如 build_from_coo）
    !   这些应当交给 SparseConvert_mod 实现。
    !============================================================
    type, extends(SparseMatrixBase) :: SparseSkyline
        logical :: is_symmetric = .true.    ! Skyline 通常用于对称矩阵
        integer :: nnz = 0                  ! 结构非零（对角 + 下三角）
        integer, allocatable :: jp(:)       ! 列指针（长度 n）
        real(8), allocatable :: al(:)       ! 下三角存储（长度 jp(n)）
        real(8), allocatable :: au(:)       ! 上三角存储（同 al）
        real(8), allocatable :: ad(:)       ! 对角线（长度 n）
    contains
        procedure :: init     => sky_init
        procedure :: to_dense => sky_to_dense
        procedure :: spmv     => sky_spmv
        procedure :: get_nnz  => sky_get_nnz

        ! Skyline 不支持 add_entry，也不支持动态修改结构

        ! 输出供调试：
        procedure :: print   => print_skyline
    end type SparseSkyline


contains
!============================================================
! 初始化（只设置 n，并清空存储）
! 实际分配 jp/al/au/ad 由构造模块完成（如从 COO 转换）
!============================================================
subroutine sky_init(self, n)
    class(SparseSkyline), intent(inout) :: self
    integer, intent(in) :: n

    self%n = n
    self%nnz = 0
    self%is_symmetric = .true.

    if (allocated(self%jp)) deallocate(self%jp)
    if (allocated(self%al)) deallocate(self%al)
    if (allocated(self%au)) deallocate(self%au)
    if (allocated(self%ad)) deallocate(self%ad)
end subroutine sky_init


!============================================================
! get_nnz：返回结构非零数（不含显式对称）
!============================================================
function sky_get_nnz(self) result(nz)
    class(SparseSkyline), intent(in) :: self
    integer :: nz
    nz = self%nnz
end function sky_get_nnz


!============================================================
! to_dense：Skyline → Dense
! 用于调试、验证格式转换
!============================================================
function sky_to_dense(self) result(A)
    class(SparseSkyline), intent(in) :: self
    real(8), allocatable :: A(:,:)

    integer :: n, j, base, jh, is, iidx, row

    n = self%n
    allocate(A(n,n))
    A = 0.0d0

    ! 对角线
    do j = 1, n
        A(j,j) = self%ad(j)
    end do

    ! 上下三角（使用 skyline 存储）
    do j = 1, n
        if (j == 1) then
            base = 0
        else
            base = self%jp(j-1)
        end if

        jh = self%jp(j) - base
        if (jh > 0) then
            is = j - jh
            do iidx = 1, jh
                row = is + (iidx - 1)

                ! 下三角元素 A(j,row)
                A(j, row) = A(j, row) + self%al(base + iidx)

                ! 上三角元素 A(row,j)
                A(row, j) = A(row, j) + self%au(base + iidx)
            end do
        end if
    end do
end function sky_to_dense


!============================================================
! SpMV: y = A * x
! 基于 skyline 存储直接计算
!============================================================
subroutine sky_spmv(self, x, y)
    class(SparseSkyline), intent(in) :: self
    real(8), intent(in)  :: x(:)
    real(8), intent(out) :: y(:)

    integer :: n, j, base, jh, is, iidx, row

    n = self%n
    y = 0.0d0

    do j = 1, n
        ! 对角项
        y(j) = y(j) + self%ad(j) * x(j)

        ! skyline 下三角与上三角
        if (j == 1) then
            base = 0
        else
            base = self%jp(j-1)
        end if

        jh = self%jp(j) - base
        if (jh > 0) then
            is = j - jh

            do iidx = 1, jh
                row = is + (iidx - 1)

                ! 下三角 A(j,row)
                y(j)   = y(j)   + self%al(base + iidx) * x(row)
                ! 上三角 A(row,j)
                y(row) = y(row) + self%au(base + iidx) * x(j)
            end do
        end if
    end do
end subroutine sky_spmv


subroutine print_skyline(self)
    class(SparseSkyline), intent(in) :: self
    
    print *, "================ Skyline Matrix ================"
    print *, "n   =", self%n
    print *, "nnz =", self%nnz

    !---------------------------------------------
    ! jp 指针
    !---------------------------------------------
    print *, "jp (column skyline pointers):"
    write(*,'( *(I10,1X) )') self%jp

    !---------------------------------------------
    ! ad（对角线）
    !---------------------------------------------
    print *, "ad (diagonal):"
    write(*,'( *(F12.4,1X) )') self%ad

    !---------------------------------------------
    ! al（下三角）
    !---------------------------------------------
    print *, "al (lower skyline block):"
    if (allocated(self%al)) then
        write(*,'( *(F12.4,1X) )') self%al
    else
        print *, "<not allocated>"
    end if

    !---------------------------------------------
    ! au（上三角）
    !---------------------------------------------
    print *, "au (upper skyline block):"
    if (allocated(self%au)) then
        write(*,'( *(F12.4,1X) )') self%au
    else
        print *, "<not allocated>"
    end if
    print *, "================================================="
end subroutine print_skyline

end module SparseSkyline_mod
