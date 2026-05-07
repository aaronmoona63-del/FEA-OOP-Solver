program test_mkl_direct_solver
    use LinearSolverBase_mod
    use SparseCSR_mod
    use MKLDirectSolver_mod
    implicit none

    !integer, parameter :: dp = kind(1.0d0)

    !========================================
    ! 测试 A：Intel 官方非对称系统（mtype=11）
    !========================================
    call test_unsymmetric()

    !========================================
    ! 测试 B：Intel 官方对称不定系统（mtype=-2）
    !========================================
    call test_symmetric()

contains

! ===============================================================
!  Test A：5×5 实非对称系统（Intel PARDISO 示例）
! ===============================================================
subroutine test_unsymmetric()
    use MatrixDense_mod
    type(SparseCSR)         :: A
    type(DenseMatrix) :: A_dense
    type(MKLDirectSolver)   :: solver
    real(dp), allocatable   :: b(:), x(:)

    integer :: i

    print *, "=============================================="
    print *, " Test A: Real Unsymmetric System (mtype=11)"
    print *, "=============================================="

    ! Matrix size
    A%n = 5
    A%nnz = 13

    print *, "A%n (before attach) =", A%n

    allocate(A%row_ptr(6))
    allocate(A%col_ind(13))
    allocate(A%val(13))

    ! Official Intel example data (unsymmetric)
    A%row_ptr = [1,4,6,9,12,14]

    A%col_ind = [ &
        1,2,4, &       ! row 1
        1,2, &         ! row 2
        3,4,5, &       ! row 3
        1,3,4, &       ! row 4
        2,5 ]          ! row 5

    A%val = [ &
         1d0, -1d0, -3d0, &
        -2d0,  5d0, &
         4d0,  6d0,  4d0, &
        -4d0,  2d0,  7d0, &
         8d0, -5d0 ]

    allocate(b(5), x(5))
    b = 1d0

    A_dense%A = A%to_dense()
    call A_dense%print()
    !==================================
    ! Setup solver
    !==================================
    call solver%attach_matrix(A)
    solver%mtype = 11       ! real unsymmetric

    call solver%analyze()
    call solver%factor()
    call solver%solve(b, x)

    print *, "Solution x:"
    do i = 1, 5
        print '(A,I1,A,1X,F12.6)', "  x(",i,") = ", x(i)
    end do

end subroutine test_unsymmetric


!===============================================================
!  Test B：8×8 对称不定系统（Intel PARDISO 示例）
!===============================================================
subroutine test_symmetric()
    use MatrixDense_mod
    type(SparseCSR)         :: A
    type(DenseMatrix):: A_dense
    type(MKLDirectSolver)   :: solver
    real(dp), allocatable   :: b(:), x(:)

    integer :: i

    print *, "=============================================="
    print *, " Test B: Real Symmetric Indefinite (mtype=-2)"
    print *, "=============================================="

    A%n = 8
    A%nnz = 18

    allocate(A%row_ptr(9))
    allocate(A%col_ind(18))
    allocate(A%val(18))

    ! Intel official symmetric CSR
    A%row_ptr = [1,5,8,10,12,15,17,18,19]

    A%col_ind = [ &
        1,3,6,7, &
        2,3,5, &
        3,8, &
        4,7, &
        5,6,7, &
        6,8, &
        7, &
        8 ]

    A%val = [ &
         7d0, 1d0, 2d0, 7d0, &
        -4d0, 8d0, 2d0, &
         1d0, 5d0, &
         7d0, 9d0, &
         5d0, 1d0, 5d0, &
        -1d0, 5d0, &
         11d0, &
         5d0 ]

    allocate(b(8), x(8))
    b = 1d0

    A_dense%A = A%to_dense()
    call A_dense%print()

    !==================================
    ! Setup solver
    !==================================
    call solver%attach_matrix(A)
    solver%mtype = 11     ! real symmetric indefinite ! -2   

    call solver%analyze()
    call solver%factor()
    call solver%solve(b, x)

    print *, "Solution x:"
    do i = 1, 8
        print '(A,I1,A,1X,F12.6)', "  x(",i,") = ", x(i)
    end do

end subroutine test_symmetric


end program test_mkl_direct_solver
