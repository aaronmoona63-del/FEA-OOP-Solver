program test_sort_random
    use Sort_mod
    implicit none

    integer, allocatable :: row(:), col(:)
    real(8), allocatable :: val(:)

    integer :: nnz, i
    integer, allocatable :: key_before(:), key_after(:)
    real(8) :: rtmp

    nnz = 1000
    allocate(row(nnz), col(nnz), val(nnz))
    allocate(key_before(nnz), key_after(nnz))

    !-----------------------------------------
    ! 生成随机 COO 数据（含重复）
    !-----------------------------------------
    do i = 1, nnz
        call random_number(rtmp) ! 0.0 ≤ x < 1.0
        row(i) = int(rtmp * 10.0d0) + 1   ! 行号 1~10
        call random_number(rtmp)
        col(i) = int(rtmp * 10.0d0) + 1   ! 列号 1~10
        call random_number(rtmp)
        val(i) = rtmp                    ! 值在 0~1
    end do

    !-----------------------------------------
    ! 生成排序前的 key = row*10000 + col
    !-----------------------------------------
    do i = 1, nnz
        key_before(i) = row(i) * 10000 + col(i)
    end do

    !-----------------------------------------
    ! 调用排序 + 去重
    !-----------------------------------------
    call COO_sort_and_unique(row, col, val, nnz)

    !-----------------------------------------
    ! 生成排序后的 key
    !-----------------------------------------
    do i = 1, nnz
        key_after(i) = row(i) * 10000 + col(i)
    end do

    !-----------------------------------------
    ! 测试 1：检查是否非降序
    !-----------------------------------------
    do i = 2, nnz
        if (key_after(i) < key_after(i-1)) then
            print *, "Sort failed: order incorrect at index ", i
            stop 10
        end if
    end do

    !-----------------------------------------
    ! 测试 2：重复项是否被合并
    !         方法：检查没有连续相同键
    !-----------------------------------------
    do i = 2, nnz
        if (row(i)==row(i-1) .and. col(i)==col(i-1)) then
            print *, "Unique failed: duplicate found at ", i
            stop 20
        end if
    end do

    !-----------------------------------------
    ! 如果运行到这里，则测试成功
    !-----------------------------------------
    print *, "test_sort_random PASSED with nnz =", nnz
end program test_sort_random