module MKLDirectSolver_mod
  use LinearSolverBase_mod, only: LinearSolver, dp
  use MatrixDense_mod,      only: DenseMatrix
  use SparseCSR_mod,        only: SparseCSR
  use PardisoDriver_mod,    only: PardisoContext
  implicit none
  private
  public :: MKLDirectSolver

  !===============================
  ! 外部 LAPACK 接口声明
  !===============================
  interface
    subroutine dgetrf(m, n, a, lda, ipiv, info)
      import dp
      integer :: m, n, lda, info
      integer :: ipiv(*)
      real(dp) :: a(lda,*)
    end subroutine dgetrf

    subroutine dgetrs(trans, n, nrhs, a, lda, ipiv, b, ldb, info)
      import dp
      character(len=1) :: trans
      integer :: n, nrhs, lda, ldb, info
      integer :: ipiv(*)
      real(dp) :: a(lda,*), b(ldb,*)
    end subroutine dgetrs
  end interface

  type, extends(LinearSolver) :: MKLDirectSolver
  type(DenseMatrix), pointer :: A_dense => null()
  type(SparseCSR),   pointer :: A_csr   => null()

  logical :: is_csr = .false.
  integer :: mtype  = 11

  character(len=16) :: profile = "robust"   ! matlab_like / robust / fast
  type(PardisoContext) :: ps

  ! ===== Dense LU cache =====
  real(dp), allocatable :: LU(:,:)      ! LU 分解结果
  integer,  allocatable :: ipiv(:)      ! pivot
  logical :: dense_factorized = .false. ! Dense 路径下的分解标志
contains
  procedure :: attach_matrix => mkl_attach
  procedure :: analyze       => mkl_analyze
  procedure :: factor        => mkl_factor
  procedure :: solve         => mkl_solve
  procedure :: free          => mkl_free
  end type MKLDirectSolver

contains

  subroutine mkl_attach(self, A)
    class(MKLDirectSolver), intent(inout) :: self
    class(*), intent(in), target :: A

    select type(pA => A)
      type is (DenseMatrix)
      self%A_dense => pA
      self%is_csr  = .false.
      type is (SparseCSR)
      self%A_csr => pA
      self%is_csr = .true.
      class default
      stop "MKLDirectSolver: attach_matrix expects DenseMatrix or SparseCSR."
    end select

    self%analyzed   = .false.
    self%factorized = .false.

    ! 切换矩阵后，dense LU 缓存失效
    self%dense_factorized = .false.
    if (allocated(self%LU))   deallocate(self%LU)
    if (allocated(self%ipiv)) deallocate(self%ipiv)

  end subroutine mkl_attach


  subroutine mkl_analyze(self)
    class(MKLDirectSolver), intent(inout) :: self
    integer :: n

    if (self%is_csr) then
      if (.not. associated(self%A_csr)) stop "MKLDirectSolver: A_csr not associated."
      n = self%A_csr%n

      call self%ps%init(self%mtype, verbosity=self%opts%verbosity, profile=self%profile)
      call self%ps%analyze(n, self%A_csr%row_ptr, self%A_csr%col_ind, self%A_csr%val)

      self%analyzed = .true.
    else
      if (.not. associated(self%A_dense)) stop "MKLDirectSolver: A_dense not associated."
      self%analyzed = .true.
    end if
  end subroutine mkl_analyze


  subroutine mkl_factor(self)
    class(MKLDirectSolver), intent(inout) :: self
    integer :: n

    if (self%is_csr) then
      if (.not. self%analyzed) stop "MKLDirectSolver: factor() before analyze()."
      n = self%A_csr%n
      call self%ps%factor(n, self%A_csr%row_ptr, self%A_csr%col_ind, self%A_csr%val)
      self%factorized = .true.
    else
      if (.not. associated(self%A_dense)) stop "MKLDirectSolver: factor() dense but A_dense not associated."
      if (.not. self%analyzed) stop "MKLDirectSolver: factor() before analyze() for dense."

      n = self%A_dense%n
      if (n <= 0) stop "MKLDirectSolver: dense matrix size n<=0."
      if (.not. allocated(self%A_dense%A)) stop "MKLDirectSolver: dense matrix A not allocated."
      if (size(self%A_dense%A,1) /= n .or. size(self%A_dense%A,2) /= n) stop "MKLDirectSolver: dense A size mismatch."

      ! allocate cache
      if (allocated(self%LU)) deallocate(self%LU)
      allocate(self%LU(n,n))
      self%LU = self%A_dense%A   ! copy, avoid overwriting original A

      if (allocated(self%ipiv)) deallocate(self%ipiv)
      allocate(self%ipiv(n))

      ! LU factorization
      block
        integer :: info
        call dgetrf(n, n, self%LU, n, self%ipiv, info)
        if (info /= 0) then
          write(*,*) "MKLDirectSolver(dense factor): dgetrf info=", info
          stop "MKLDirectSolver: dense LU factorization failed."
        end if
      end block

      self%dense_factorized = .true.
      self%factorized = .true.
    end if
  end subroutine mkl_factor


  subroutine mkl_solve(self, b, x)
    class(MKLDirectSolver), intent(inout) :: self
    real(dp), intent(inout)  :: b(:)
    real(dp), intent(out) :: x(:)

    integer :: n

    if (self%is_csr) then
      if (.not. self%factorized) stop "MKLDirectSolver: solve() before factor()."
      n = self%A_csr%n
      call self%ps%solve(n, self%A_csr%row_ptr, self%A_csr%col_ind, self%A_csr%val, b, x)
    else
      block
        integer :: nrhs, info, ldb
        if (.not. associated(self%A_dense)) stop "MKLDirectSolver: solve() dense but A_dense not associated."
        if (.not. self%dense_factorized) stop "MKLDirectSolver: solve() dense before factor()."

        n = self%A_dense%n
        if (size(b) /= n .or. size(x) /= n) stop "MKLDirectSolver: dense solve size mismatch."

        nrhs = 1
        ldb  = n

        ! dgetrs 会把 RHS 覆盖成解，所以我们先拷贝到 x
        x = b
        call dgetrs('N', n, nrhs, self%LU, n, self%ipiv, x, ldb, info)
        if (info /= 0) then
          write(*,*) "MKLDirectSolver(dense solve): dgetrs info=", info
          stop "MKLDirectSolver: dense solve failed."
        end if
      end block
    end if
  end subroutine mkl_solve


 subroutine mkl_free(self)
    class(MKLDirectSolver), intent(inout) :: self
    integer :: n

    if (.not. self%analyzed) return

    ! ==========================================
    ! 针对 MKL Pardiso phase=-1 在部分 Linux 系统上的释放 Bug：
    ! 我们注释掉显式的 finalize 调用，防止它触发 Double Free
    ! 内存将由操作系统在进程结束时 100% 安全地自动回收
    ! ==========================================
    ! if (self%is_csr .and. associated(self%A_csr)) then
    !   n = self%A_csr%n
    !   call self%ps%finalize(n)
    ! else
    !   call self%ps%finalize(0)
    ! end if

    nullify(self%A_dense, self%A_csr)
    self%analyzed   = .false.
    self%factorized = .false.

    self%dense_factorized = .false.
    if (allocated(self%LU))   deallocate(self%LU)
    if (allocated(self%ipiv)) deallocate(self%ipiv)
  end subroutine mkl_free
  end module MKLDirectSolver_mod