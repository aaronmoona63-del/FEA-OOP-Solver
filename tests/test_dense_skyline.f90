program test_dense2skyline_full
    use MatrixDense_mod
    use SparseSkyline_mod
    implicit none

    integer :: n
    type(DenseMatrix) :: Ksym
    type(SparseSkyline) :: Sky
    double precision, allocatable :: al(:), au(:), ad(:)
    integer, allocatable :: jp(:)
    double precision, allocatable :: Arec(:,:), x(:), y(:)

    print *, "==========================================="
    print *, "   Test Dense → Skyline & SparseSkyline"
    print *, "==========================================="

    !-----------------------------------------------
    ! 1) 构建 dense 矩阵（对称）
    !-----------------------------------------------
    n = 5
    call Ksym%init(n, .true.)

    Ksym%A = reshape([ &
        1.0d0, 0.2d0, 0.4d0, 0.1d0, 0.0d0, &
        0.2d0, 2.0d0, 0.3d0, 0.0d0, 0.0d0, &
        0.4d0, 0.3d0, 3.0d0, 0.0d0, 0.5d0, &
        0.1d0, 0.0d0, 0.0d0, 4.0d0, 0.7d0, &
        0.0d0, 0.0d0, 0.5d0, 0.7d0, 5.0d0 ], shape=[5,5])

    print *, ">>> Original Dense matrix:"
    call Ksym%print()

    !-----------------------------------------------
    ! 2) Dense → Skyline（通过 DenseMatrix_mod）
    !-----------------------------------------------
    print *, " "
    print *, ">>> Converting Dense → Skyline (MatrixDense_mod)..."

    call Ksym%to_skyline(al, au, ad, jp)

    print *, "jp = ", jp
    print *, "ad = ", ad
    print *, "al = ", al
    print *, "au = ", au

    !-----------------------------------------------
    ! 3) 构造 SparseSkyline 对象（测试 SparseSkyline_mod）
    !-----------------------------------------------
    print *, " "
    print *, ">>> Build SparseSkyline from Dense-produced skyline arrays"

    call Sky%init(n)

    ! 直接写入 skyline 结构
    Sky%jp = jp
    Sky%ad = ad
    Sky%al = al
    Sky%au = au
    Sky%is_symmetric = .true.
    Sky%nnz = size(al) + n

    !-----------------------------------------------
    ! 4) SparseSkyline → Dense 再还原（验证正确性）
    !-----------------------------------------------
    print *, " "
    print *, ">>> Convert Skyline → Dense (SparseSkyline_mod)"

    Arec = Sky%to_dense()

    print *, "Recovered Dense matrix from Skyline:"
    call print_matrix(Arec)

    !-----------------------------------------------
    ! 5) SpMV 测试： y = A*x
    !-----------------------------------------------
    allocate(x(n), y(n))
    x = (/ 1d0, 2d0, 3d0, 4d0, 5d0 /)

    call Sky%spmv(x, y)

    print *, " "
    print *, ">>> Test SpMV  y = A*x"
    print *, "x =", x
    print *, "y =", y

    print *, "==========================================="
    print *, "        Dense → Skyline tests done."
    print *, "==========================================="

contains

    !-----------------------------------------------
    ! 辅助打印矩阵
    !-----------------------------------------------
    subroutine print_matrix(A)
        double precision, intent(in) :: A(:,:)
        integer :: i, j
        do i = 1, size(A,1)
            write(*,'(100F10.3)') (A(i,j), j=1,size(A,2))
        end do
    end subroutine print_matrix

end program test_dense2skyline_full
