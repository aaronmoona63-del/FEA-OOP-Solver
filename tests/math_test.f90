! ==============================================================================
!  项目名称: Fortran 核心求解器底层数学校验
!  功能描述: 使用 Pardiso 求解 5x5 异质病态矩阵，对标 MATLAB 结果
! ==============================================================================
program math_test
    implicit none
    external :: pardisoinit, pardiso

    integer(8) :: pt(64) 
    integer :: iparm(64)
    integer :: maxfct, mnum, mtype, phase, n, nrhs, error, msglvl
    integer :: i
    integer, allocatable :: perm(:) 

    ! 动态数组：稀疏矩阵 (CSR 格式) 及 向量
    integer, allocatable :: ia(:), ja(:)
    real(8), allocatable :: a(:), b(:), x(:)

    ! --- 1. 问题初始化 ---
    n = 5
    nrhs = 1
    maxfct = 1
    mnum = 1
    
    ! 分配内存 (ia 长度 n+1, ja 和 a 长度为 13)
    allocate(ia(n+1), ja(13), a(13), perm(n), b(n), x(n))

    ! --- 2. 灌入 5x5 异质病态矩阵 (CSR 格式) ---
    ia = [1, 3, 6, 9, 12, 14] 
    
    ja = [1, 2, &                 ! 第1行
          1, 2, 3, &              ! 第2行
          2, 3, 4, &              ! 第3行
          3, 4, 5, &              ! 第4行
          4, 5]                   ! 第5行
          
    a  = [ 2000.0d0, -1000.0d0, &                     ! 第1行
          -1000.0d0,  1001.0d0,    -1.0d0, &          ! 第2行
             -1.0d0,     2.0d0,    -1.0d0, &          ! 第3行
             -1.0d0,  1001.0d0, -1000.0d0, &          ! 第4行
          -1000.0d0,  1000.0d0]                       ! 第5行

    ! --- 3. 设置右端项 b ---
    b = [100.0d0, 0.0d0, 50.0d0, 0.0d0, 200.0d0]
    x = 0.0d0
    perm = 0

    ! --- 4. 配置 Pardiso (mtype=11 非对称/全元素实数矩阵) ---
    mtype = 11  
    pt = 0 
    iparm = 0 
    call pardisoinit(pt, mtype, iparm)

    iparm(1)  = 1    
    iparm(3)  = 1    
    iparm(28) = 0    
    iparm(35) = 0    

    ! --- 5. 执行求解 ---
    print *, "=== Fortran Core Solver Verification ==="
    print *, "Matrix Dimension: ", n, ", Non-zeros: ", 13
    print *, "Calling Pardiso Direct Solver..."
    print *, "----------------------------------------"
    
    phase = 13       
    msglvl = 0       
    call pardiso(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
                 perm, nrhs, iparm, msglvl, b, x, error)

    if (error /= 0) then
        print *, "Error: Pardiso failed with error code ", error
        stop
    end if

    ! --- 6. 格式化输出结果 ---
    print *, "Solution vector x:"
    do i = 1, n
        write(*, '(A,I1,A,F15.10)') "  x(", i, ") = ", x(i)
    end do
    print *, "========================================"

    ! --- 7. 内存清理 ---
    phase = -1 
    call pardiso(pt, maxfct, mnum, mtype, phase, n, a, ia, ja, &
                 perm, nrhs, iparm, msglvl, b, x, error)
    
    deallocate(ia, ja, a, perm, b, x)

end program math_test