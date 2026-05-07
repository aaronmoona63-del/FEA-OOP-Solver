program test_sort_basic
    use Sort_mod
    implicit none

    integer, allocatable :: row(:), col(:)
    real(8), allocatable :: val(:)
    integer :: nnz

    nnz = 8
    allocate(row(nnz), col(nnz), val(nnz))
    row=[3, 1, 2, 3, 1, 2, 3, 1]
    col=[2, 3, 1, 2, 3, 1, 2, 3]
    val=[1, 1, 1, 2, 2, 2, 3, 3]

    call COO_sort_and_unique(row, col, val, nnz)

    ! === EXPECTED RESULTS ===
    ! (1,3) → val=6
    ! (2,1) → val=3
    ! (3,2) → val=6

    if (nnz /= 3) stop 1
    if (row(1) /= 1 .or. col(1) /= 3 .or. abs(val(1)-6.0d0) > 1e-12) stop 2
    if (row(2) /= 2 .or. col(2) /= 1 .or. abs(val(2)-3.0d0) > 1e-12) stop 3
    if (row(3) /= 3 .or. col(3) /= 2 .or. abs(val(3)-6.0d0) > 1e-12) stop 4

end program test_sort_basic
