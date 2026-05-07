program test_coo_csr_skyline
    use SparseMatrixBase_mod
    use SparseCOO_mod
    use SparseCSR_mod
    use SparseSkyline_mod
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    integer :: n, i, j
    real(dp), allocatable :: K(:,:)         ! 参考稠密矩阵
    type(SparseCOO)        :: C             ! 原始 COO
    type(SparseCSR)        :: S             ! 原始 CSR
    type(SparseSkyline)    :: Sky           ! 原始 Skyline

    class(SparseMatrixBase), pointer :: B   ! 多态指针

    ! 各种格式恢复成的 dense
    real(dp), allocatable :: A_ref(:,:), A_csr(:,:), A_csr_coo(:,:)
    real(dp), allocatable :: A_sky(:,:), A_sky_coo(:,:), A_sky_csr(:,:)

    ! SpMV 测试
    real(dp), allocatable :: x(:), y_ref(:), y_coo(:), y_csr(:), y_sky(:)
    real(dp) :: err

    print *, "==============================================="
    print *, "  COO ↔ CSR ↔ Skyline 全流程一致性测试"
    print *, "==============================================="

    !-------------------------------------------
    ! 1. 构造参考稠密矩阵 K（与前面例子相同）
    !-------------------------------------------
    n = 5
    allocate(K(n,n))

    K = reshape([ &
        1.0d0, 0.2d0, 0.4d0, 0.1d0, 0.0d0, &
        0.2d0, 2.0d0, 0.3d0, 0.0d0, 0.0d0, &
        0.4d0, 0.3d0, 3.0d0, 0.0d0, 0.5d0, &
        0.1d0, 0.0d0, 0.0d0, 4.0d0, 0.7d0, &
        0.0d0, 0.0d0, 0.5d0, 0.7d0, 5.0d0 ], shape=[n,n])

    print *, ">>> 原始稠密矩阵 K:"
    call print_matrix(K)

    !-------------------------------------------
    ! 2. Dense → COO
    !-------------------------------------------
    print *, " "
    print *, ">>> 构造 SparseCOO（Dense → COO）"

    call C%init(n)
    do i = 1, n
        do j = 1, n
            if (K(i,j) /= 0.0d0) then
                call C%add_entry(i, j, K(i,j))
            end if
        end do
    end do
    print *, "COO nnz =", C%get_nnz()

    A_ref = C%to_dense()
    print *, "COO → Dense (A_ref):"
    call print_matrix(A_ref)

    !-------------------------------------------
    ! 3. COO → CSR，并检查
    !-------------------------------------------
    print *, " "
    print *, ">>> COO → CSR（显式 build_from_coo）"

    call S%init(n)
    call S%build_from_coo(C)

    A_csr = S%to_dense()
    err = maxval(abs(A_csr - A_ref))
    print *, "CSR → Dense 与参考差值 max|A_csr - A_ref| =", err

    !-------------------------------------------
    ! 4. CSR → COO，并检查
    !-------------------------------------------
    print *, " "
    print *, ">>> CSR → COO（调用 S%to_coo，多态返回）"

    B => S%to_coo()
    select type (C2 => B)
    type is (SparseCOO)
        A_csr_coo = C2%to_dense()
        err = maxval(abs(A_csr_coo - A_ref))
        print *, "CSR→COO→Dense 与参考差值 max|A_csr_coo - A_ref| =", err
    class default
        print *, "ERROR: S%to_coo() 没有返回 SparseCOO 类型！"
        stop
    end select

    !-------------------------------------------
    ! 5. COO → Skyline，并检查
    !-------------------------------------------
    print *, " "
    print *, ">>> COO → Skyline（调用 build_from_coo）"

    call Sky%init(n)
    call Sky%build_from_coo(C, .true.)    ! 当前假设对称

    A_sky = Sky%to_dense()
    err = maxval(abs(A_sky - A_ref))
    print *, "Skyline → Dense 与参考差值 max|A_sky - A_ref| =", err

    !-------------------------------------------
    ! 6. Skyline → COO，并检查
    !-------------------------------------------
    print *, " "
    print *, ">>> Skyline → COO（Sky%to_coo，多态返回）"

    B => Sky%to_coo()
    select type (C3 => B)
    type is (SparseCOO)
        A_sky_coo = C3%to_dense()
        err = maxval(abs(A_sky_coo - A_ref))
        print *, "Skyline→COO→Dense 差值 max|A_sky_coo - A_ref| =", err
    class default
        print *, "ERROR: Sky%to_coo() 没有返回 SparseCOO 类型！"
        stop
    end select

    !-------------------------------------------
    ! 7. Skyline → CSR，并检查
    !-------------------------------------------
    print *, " "
    print *, ">>> Skyline → CSR（Sky%to_csr，多态返回）"

    B => Sky%to_csr()
    select type (S2 => B)
    type is (SparseCSR)
        A_sky_csr = S2%to_dense()
        err = maxval(abs(A_sky_csr - A_ref))
        print *, "Skyline→CSR→Dense 差值 max|A_sky_csr - A_ref| =", err
    class default
        print *, "WARNING: Sky%to_csr() 暂未实现有效返回（ptr 为空？）"
    end select

    !-------------------------------------------
    ! 8. SpMV 一致性测试： y = A * x
    !-------------------------------------------
    print *, " "
    print *, ">>> SpMV 一致性测试：y = A * x"

    allocate(x(n), y_ref(n), y_coo(n), y_csr(n), y_sky(n))
    x = (/ 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0 /)

    ! 参考结果：dense matmul
    y_ref = matmul(A_ref, x)

    call C%spmv(x, y_coo)
    call S%spmv(x, y_csr)
    call Sky%spmv(x, y_sky)

    print *, "y_ref (Dense) :", y_ref
    print *, "y_coo (COO)   :", y_coo
    print *, "y_csr (CSR)   :", y_csr
    print *, "y_sky (Sky)   :", y_sky

    print *, "max|y_coo - y_ref| =", maxval(abs(y_coo - y_ref))
    print *, "max|y_csr - y_ref| =", maxval(abs(y_csr - y_ref))
    print *, "max|y_sky - y_ref| =", maxval(abs(y_sky - y_ref))

    print *, "==============================================="
    print *, "  COO ↔ CSR ↔ Skyline 全流程一致性测试结束"
    print *, "==============================================="

contains

    subroutine print_matrix(A)
        real(dp), intent(in) :: A(:,:)
        integer :: i, j
        do i = 1, size(A,1)
            write(*,'(100F10.3)') (A(i,j), j=1,size(A,2))
        end do
    end subroutine print_matrix

end program test_coo_csr_skyline
