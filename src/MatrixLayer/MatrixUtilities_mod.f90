module MatrixUtilities_mod
!!
!!  通用工具模块：矩阵结构算法的基础工具
!!  -----------------------------------------
!!  提供：
!!     * merge_sort_2i_r8    —— 三数组(i,j,val)按(i,j)排序
!!     * unique_merge_2i_r8  —— 合并重复(i,j)，val累加
!!     * merge_sort_i_r8     —— 一数组(i,val)按i排序（CSR行内用）
!!
!!  完全独立，不依赖 SparseCOO / CSR / Skyline
!!
implicit none
private

public :: merge_sort_2i_r8
public :: unique_merge_2i_r8
public :: merge_sort_i_r8

contains
!===============================================================
!  归并排序（通用版）：按 (i,j) 升序排序
!===============================================================
recursive subroutine merge_sort_2i_r8(i1, i2, r8, L, R)
    integer, intent(inout) :: i1(:), i2(:)
    real(8), intent(inout) :: r8(:)
    integer, intent(in)    :: L, R

    integer :: mid
    if (R - L + 1 <= 16) then
        call insertion_sort_2i_r8(i1, i2, r8, L, R)
        return
    end if

    mid = (L + R)/2
    call merge_sort_2i_r8(i1, i2, r8, L, mid)
    call merge_sort_2i_r8(i1, i2, r8, mid+1, R)
    call merge_blocks_2i_r8(i1, i2, r8, L, mid, R)
end subroutine merge_sort_2i_r8


subroutine insertion_sort_2i_r8(i1, i2, r8, L, R)
    integer, intent(inout) :: i1(:), i2(:)
    real(8), intent(inout) :: r8(:)
    integer, intent(in)    :: L, R

    integer :: k, j, ti1, ti2
    real(8) :: tv

    do k = L+1, R
        ti1 = i1(k)
        ti2 = i2(k)
        tv  = r8(k)
        j = k - 1

        do while (j >= L .and. ( i1(j) > ti1 .or. (i1(j)==ti1 .and. i2(j) > ti2) ))
            i1(j+1) = i1(j)
            i2(j+1) = i2(j)
            r8(j+1) = r8(j)
            j = j - 1
        end do

        i1(j+1) = ti1
        i2(j+1) = ti2
        r8(j+1) = tv
    end do
end subroutine insertion_sort_2i_r8


subroutine merge_blocks_2i_r8(i1, i2, r8, L, mid, R)
    integer, intent(inout) :: i1(:), i2(:)
    real(8), intent(inout) :: r8(:)
    integer, intent(in)    :: L, mid, R

    integer :: n, k, p, q
    integer, allocatable :: t1(:), t2(:)
    real(8), allocatable :: tv(:)

    n = R - L + 1
    allocate(t1(n), t2(n), tv(n))

    p = L
    q = mid + 1
    k = 1

    do while (p <= mid .and. q <= R)
        if ( i1(p) < i1(q) .or. (i1(p)==i1(q) .and. i2(p) <= i2(q)) ) then
            t1(k)=i1(p); t2(k)=i2(p); tv(k)=r8(p)
            p = p + 1
        else
            t1(k)=i1(q); t2(k)=i2(q); tv(k)=r8(q)
            q = q + 1
        end if
        k = k + 1
    end do

    do while (p <= mid)
        t1(k)=i1(p); t2(k)=i2(p); tv(k)=r8(p)
        k = k + 1; p = p + 1
    end do
    do while (q <= R)
        t1(k)=i1(q); t2(k)=i2(q); tv(k)=r8(q)
        k = k + 1; q = q + 1
    end do

    i1(L:R) = t1
    i2(L:R) = t2
    r8(L:R) = tv
end subroutine merge_blocks_2i_r8


!===============================================================
!  unique_merge_2i_r8：合并重复 (i1,i2)，r8 相加
!===============================================================
subroutine unique_merge_2i_r8(i1, i2, r8, nnz)
    integer, intent(inout) :: i1(:), i2(:)
    real(8), intent(inout) :: r8(:)
    integer, intent(inout) :: nnz

    integer :: k, p

    if (nnz <= 1) return

    p = 1
    do k = 2, nnz
        if (i1(k)==i1(p) .and. i2(k)==i2(p)) then
            r8(p) = r8(p) + r8(k)
        else
            p = p + 1
            i1(p) = i1(k)
            i2(p) = i2(k)
            r8(p) = r8(k)
        end if
    end do

    nnz = p
end subroutine unique_merge_2i_r8


!===============================================================
!  merge_sort_i_r8：用于 CSR 行内 (col,val) 排序
!===============================================================
recursive subroutine merge_sort_i_r8(i1, r8, L, R)
    integer, intent(inout) :: i1(:)
    real(8), intent(inout) :: r8(:)
    integer, intent(in)    :: L, R

    integer :: mid
    if (R-L+1 <= 16) then
        call insertion_sort_i_r8(i1, r8, L, R)
        return
    end if

    mid = (L+R)/2
    call merge_sort_i_r8(i1, r8, L, mid)
    call merge_sort_i_r8(i1, r8, mid+1, R)
    call merge_blocks_i_r8(i1, r8, L, mid, R)
end subroutine merge_sort_i_r8


subroutine insertion_sort_i_r8(i1, r8, L, R)
    integer, intent(inout) :: i1(:)
    real(8), intent(inout) :: r8(:)
    integer, intent(in)    :: L, R

    integer :: k, j, ti
    real(8) :: tv

    do k = L+1, R
        ti = i1(k)
        tv = r8(k)
        j  = k - 1
        do while (j >= L .and. i1(j) > ti)
            i1(j+1) = i1(j)
            r8(j+1) = r8(j)
            j = j - 1
        end do
        i1(j+1) = ti
        r8(j+1) = tv
    end do
end subroutine insertion_sort_i_r8


subroutine merge_blocks_i_r8(i1, r8, L, mid, R)
    integer, intent(inout) :: i1(:)
    real(8), intent(inout) :: r8(:)
    integer, intent(in)    :: L, mid, R

    integer :: n, p, q, k
    integer, allocatable :: t1(:)
    real(8), allocatable :: tv(:)

    n = R - L + 1
    allocate(t1(n), tv(n))

    p = L
    q = mid + 1
    k = 1

    do while (p <= mid .and. q <= R)
        if (i1(p) <= i1(q)) then
            t1(k) = i1(p)
            tv(k) = r8(p)
            p = p + 1
        else
            t1(k) = i1(q)
            tv(k) = r8(q)
            q = q + 1
        end if
        k = k + 1
    end do

    do while (p <= mid)
        t1(k) = i1(p)
        tv(k) = r8(p)
        k = k + 1
        p = p + 1
    end do

    do while (q <= R)
        t1(k) = i1(q)
        tv(k) = r8(q)
        k = k + 1
        q = q + 1
    end do

    i1(L:R) = t1
    r8(L:R) = tv
end subroutine merge_blocks_i_r8

end module MatrixUtilities_mod
