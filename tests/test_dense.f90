program test_dense
    use MatrixDense_mod
    implicit none

    type(DenseMatrix) :: Ksym, Kunsym
    integer :: n
    double precision, allocatable :: al(:), au(:), ad(:)
    integer, allocatable :: jp(:)

    print *, "==========================================="
    print *, "  Test DenseMatrix → Skyline Conversion"
    print *, "==========================================="

    ! -----------------------------------------------------
    ! 1) 测试 **对称矩阵**（和你的 MATLAB 示例一致）
    ! -----------------------------------------------------
    n = 5
    call Ksym%init(n, .true.)   ! 对称矩阵

    Ksym%A = reshape([ &
        1.0d0, 0.2d0, 0.4d0, 0.1d0, 0.0d0, &
        0.2d0, 2.0d0, 0.3d0, 0.0d0, 0.0d0, &
        0.4d0, 0.3d0, 3.0d0, 0.0d0, 0.5d0, &
        0.1d0, 0.0d0, 0.0d0, 4.0d0, 0.7d0, &
        0.0d0, 0.0d0, 0.5d0, 0.7d0, 5.0d0 ], shape=[n,n])

    print *, ">>> Symmetric Dense Matrix Ksym:"
    call Ksym%print()

    ! --- 转换到 skyline ---
    call Ksym%to_skyline(al, au, ad, jp)

    print *, " "
    print *, ">>> Skyline jp:"
    print *, jp
    print *, ">>> ad:"
    print *, ad
    print *, ">>> al:"
    print *, al
    print *, ">>> au:"
    print *, au


    ! -----------------------------------------------------
    ! 2) 测试 **非对称矩阵**
    !   这里也使用你的 MATLAB 示例的非对称矩阵
    ! -----------------------------------------------------
    call Kunsym%init(n, .false.)   ! 非对称矩阵

    Kunsym%A = reshape([ &
        1.0d0, 0.2d0, 0.4d0, 0.1d0, 0.0d0, &
        0.3d0, 2.0d0, 0.3d0, 0.0d0, 0.0d0, &
        0.4d0, 0.3d0, 3.0d0, 0.0d0, 0.5d0, &
        0.1d0, 3.0d0, 0.0d0, 4.0d0, 0.7d0, &
        0.0d0, 3.0d0, 0.5d0, 0.7d0, 5.0d0 ], shape=[5,5])

    print *, " "
    print *, ">>> Unsymmetric Dense Matrix Kunsym:"
    call Kunsym%print()

    ! --- 转换到 skyline ---
    call Kunsym%to_skyline(al, au, ad, jp)

    print *, " "
    print *, ">>> Skyline jp:"
    print *, jp
    print *, ">>> ad:"
    print *, ad
    print *, ">>> al:"
    print *, al
    print *, ">>> au:"
    print *, au

    print *, "==========================================="
    print *, "  Dense → Skyline tests completed."
    print *, "==========================================="

end program test_dense