module MatrixDense_mod
  use MatrixBase_mod
  implicit none
  private
  public :: DenseMatrix

  type, extends(MatrixBase) :: DenseMatrix
  double precision, allocatable :: A(:,:)
contains
  procedure :: init
  procedure :: get_size
  procedure :: print
  procedure :: to_skyline
  end type DenseMatrix


contains

  ! ==============================
  ! Init dense matrix
  ! ==============================
  subroutine init(self, n, is_symmetric)
    class(DenseMatrix), intent(inout) :: self
    integer, intent(in) :: n
    logical, intent(in) :: is_symmetric

    self%n = n
    self%is_symmetric = is_symmetric
    allocate(self%A(n,n))
    self%A(:,:) = 0.0d0
  end subroutine init


  ! ==============================
  ! Get matrix size
  ! ==============================
  function get_size(self) result(n)
    class(DenseMatrix), intent(in) :: self
    integer :: n
    n = self%n
  end function get_size


  ! ==============================
  ! Print matrix
  ! ==============================
  subroutine print(self)
    class(DenseMatrix), intent(in) :: self
    integer :: i, j

    print *, "DenseMatrix size =", self%n
    ! 打印列号
    write(*,'(6X, *(I10,1x))') (j, j = 1, size(self%A,2))
    ! 打印矩阵内容（带行号）
    do i = 1, size(self%A,1)
      write(*,'(I4,2X, *(F10.4,1X))') i, self%A(i,:)
    end do

  end subroutine print


  ! ==============================
  ! Convert dense → skyline
  ! ==============================
  subroutine to_skyline(self, al, au, ad, jp)
    class(DenseMatrix), intent(in) :: self
    double precision, allocatable, intent(out) :: al(:), au(:), ad(:)
    integer, allocatable, intent(out) :: jp(:)

    integer :: n
    integer :: j, i, is_col, is_row, is, jh
    integer :: total, lastend, thisend, k

    n = self%n
    allocate(jp(n))

    ! ===== Step 1: Compute skyline profile =====
    total = 0
    do j = 1, n
      is_col = j
      is_row = j

      do i = 1, j-1
        if (self%A(i,j) /= 0.0d0) then
          is_col = i
          exit
        end if
      end do

      do i = 1, j-1
        if (self%A(j,i) /= 0.0d0) then
          is_row = i
          exit
        end if
      end do

      is = min(is_col, is_row)
      jh = j - is
      total = total + jh
      jp(j) = total
    end do

    ! Allocate skyline arrays
    allocate(ad(n))
    allocate(al(total))
    allocate(au(total))

    ! ===== Step 2: Fill skyline storage =====
    ad(1) = self%A(1,1)

    do j = 2, n
      ad(j) = self%A(j,j)
      lastend = jp(j-1)
      thisend = jp(j)
      jh = thisend - lastend
      if (jh <= 0) cycle

      is = j - jh
      k = lastend + 1

      ! Fill AU (column)
      do i = is, j-1
        au(k) = self%A(i,j)
        k = k + 1
      end do

      ! Fill AL (row)
      k = lastend + 1
      if (.not. self%is_symmetric) then
        do i = is, j-1
          al(k) = self%A(j,i)
          k = k + 1
        end do
      else
        do i = is, j-1
          al(k) = self%A(i,j)
          k = k + 1
        end do
      end if
    end do

  end subroutine to_skyline

end module MatrixDense_mod
