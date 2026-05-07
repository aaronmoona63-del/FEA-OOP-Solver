module SparseCSR_mod
    use SparseMatrixBase_mod
    implicit none
    private
    public :: SparseCSR

    integer, parameter :: dp = kind(1.0d0)

    type, extends(SparseMatrixBase) :: SparseCSR
        integer :: nnz = 0
        integer, allocatable :: row_ptr(:)
        integer, allocatable :: col_ind(:)
        real(dp), allocatable :: val(:)
    contains
        procedure :: init       => csr_init
        procedure :: to_dense   => csr_to_dense
        procedure :: spmv       => csr_spmv
        procedure :: get_nnz    => csr_get_nnz

        ! 输出供调试：
        procedure :: print   => print_csr
    end type SparseCSR
contains

subroutine csr_init(self, n)
    class(SparseCSR), intent(inout) :: self
    integer, intent(in) :: n

    self%n = n
    self%nnz = 0

    if (allocated(self%row_ptr)) deallocate(self%row_ptr)
    if (allocated(self%col_ind)) deallocate(self%col_ind)
    if (allocated(self%val))     deallocate(self%val)
end subroutine csr_init

!---------------------------------------------------------
function csr_get_nnz(self) result(nz)
    class(SparseCSR), intent(in) :: self
    integer :: nz
    nz = self%nnz
end function csr_get_nnz

!---------------------------------------------------------
function csr_to_dense(self) result(A)
    class(SparseCSR), intent(in) :: self
    real(dp), allocatable :: A(:,:)
    integer :: i, k, s, e

    allocate(A(self%n, self%n))
    A = 0.0_dp

    do i = 1, self%n
        s = self%row_ptr(i)
        e = self%row_ptr(i+1) - 1
        do k = s, e
            A(i, self%col_ind(k)) = A(i, self%col_ind(k)) + self%val(k)
        end do
    end do
end function csr_to_dense

!---------------------------------------------------------
subroutine csr_spmv(self, x, y)
    class(SparseCSR), intent(in) :: self
    real(dp), intent(in)  :: x(:)
    real(dp), intent(out) :: y(:)
    integer :: i, k, s, e

    y = 0.0_dp
    do i = 1, self%n
        s = self%row_ptr(i)
        e = self%row_ptr(i+1) - 1
        do k = s, e
            y(i) = y(i) + self%val(k) * x(self%col_ind(k))
        end do
    end do
end subroutine csr_spmv

subroutine print_csr(self)
    class(SparseCSR), intent(in) :: self

    print *, "================ CSR (Compressed Sparse Row) ================"
    print *, "n   =", self%n
    print *, "nnz =", self%nnz

    if (self%nnz == 0) then
        print *, "<empty>"
        return
    end if

    !---------------------------------------------
    ! 打印 row_ptr
    !---------------------------------------------
    print *, "row_ptr:"
    write(*,'( *(I10,1X) )') self%row_ptr
    !---------------------------------------------
    ! 打印 col_ind
    !---------------------------------------------
    print *, "col_ind:"
    write(*,'( *(I10,1X) )') self%col_ind
    !---------------------------------------------
    ! 打印 val
    !---------------------------------------------
    print *, "val:"
    write(*,'( *(F12.4,1X) )') self%val
    print *, "============================================================="
end subroutine print_csr

end module SparseCSR_mod
