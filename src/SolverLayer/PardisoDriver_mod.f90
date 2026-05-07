include 'mkl_pardiso.f90'
module PardisoDriver_mod
  ! 统一封装 PARDISO：init/analyze/factor/solve/finalize
  !
  ! 要求：编译时能 include 到 mkl_pardiso.f90
  ! 例如 ifort ... -I${MKLROOT}/include!
  ! 

  use mkl_pardiso
  !include 'mkl_pardiso.f90'

  implicit none
  private
  public :: PardisoContext

  integer, parameter :: dp = kind(1.0d0)

  type :: PardisoContext
  type(MKL_PARDISO_HANDLE) :: pt(64)
  integer, allocatable :: iparm(:)

  integer :: mtype  = 0
  integer :: maxfct = 1
  integer :: mnum   = 1
  integer :: msglvl = 0
  integer :: nrhs_default = 1

  logical :: initialized = .false.
  logical :: analyzed    = .false.
  logical :: factorized  = .false.

contains
  procedure :: init        => pardiso_init
  procedure :: set_profile => pardiso_set_profile
  procedure :: analyze     => pardiso_analyze
  procedure :: factor      => pardiso_factor
  procedure :: solve       => pardiso_solve
  procedure :: finalize    => pardiso_finalize
  end type PardisoContext

contains

  subroutine pardiso_init(self, mtype, verbosity, profile)
    class(PardisoContext), intent(inout) :: self
    integer, intent(in) :: mtype
    integer, intent(in), optional :: verbosity
    character(len=*), intent(in), optional :: profile

    integer :: i

    if (.not. allocated(self%iparm)) allocate(self%iparm(64))
    self%iparm = 0

    self%mtype = mtype

    if (present(verbosity)) then
      self%msglvl = merge(1, 0, verbosity > 0)
    else
      self%msglvl = 0
    end if

    ! init pt
    do i = 1, 64
      self%pt(i)%DUMMY = 0
    end do

    ! default profile
    if (present(profile)) then
      call self%set_profile(profile)
    else
      call self%set_profile("robust")
    end if

    self%initialized = .true.
    self%analyzed    = .false.
    self%factorized  = .false.
  end subroutine pardiso_init


  subroutine pardiso_set_profile(self, profile)
    class(PardisoContext), intent(inout) :: self
    character(len=*), intent(in) :: profile

    if (.not. allocated(self%iparm)) allocate(self%iparm(64))
    self%iparm = 0

    ! 共同基础：不使用 solver default（我们自己设 iparm）
    self%iparm(1)  = 1   ! no solver default
    self%iparm(2)  = 2   ! METIS reorder
    self%iparm(3)  = 1   ! threads (1 means MKL decides or single; 你也可以外部设置 OMP)
    self%iparm(4)  = 0   ! no iterative-direct
    self%iparm(5)  = 0
    self%iparm(6)  = 0
    self%iparm(8)  = 2   ! iterative refinement steps (solve phase)

    self%iparm(18) = -1  ! report nnz in LU
    self%iparm(19) = -1  ! report MFLOPS
    self%iparm(20) = 0

    select case (adjustl(profile))
      case ("matlab_like", "matlab", "clean")
      ! 尽量贴近 MATLAB 的“朴素 LU”风格：关闭 scaling/matching/扰动
      self%iparm(10) = 0   ! no pivot perturbation
      self%iparm(11) = 0   ! no MPS scaling
      self%iparm(13) = 0   ! no weighted matching
      case ("robust", "stable")
      ! 更稳健：开启 scaling/matching/轻微扰动
      self%iparm(10) = 13  ! pivot perturbation 1e-13
      self%iparm(11) = 1   ! MPS scaling
      self%iparm(13) = 1   ! weighted matching
      case ("fast")
      ! 偏快：保留 reorder，关闭 matching，少做“花活”
      self%iparm(10) = 0
      self%iparm(11) = 1
      self%iparm(13) = 0
      case default
      stop "PardisoContext: unknown profile (use matlab_like/robust/fast)."
    end select
  end subroutine pardiso_set_profile


  subroutine pardiso_analyze(self, n, ia, ja, a)
    class(PardisoContext), intent(inout) :: self
    integer, intent(in) :: n
    integer, intent(in) :: ia(:), ja(:)
    real(dp), intent(inout) :: a(:)

    integer :: phase, error
    integer :: idum(1)
    real(dp) :: ddum(1)

    if (.not. self%initialized) stop "PardisoContext: call init() first."

    phase = 11
    error = 0
    idum(1) = 0
    ddum(1) = 0.0_dp

    call pardiso(self%pt, self%maxfct, self%mnum, self%mtype, phase, n, &
      a, ia, ja, idum, 0, self%iparm, self%msglvl, ddum, ddum, error)

    if (error /= 0) then
      write(*,*) "PARDISO analyze error=", error
      stop "PardisoContext: analyze failed."
    end if

    self%analyzed = .true.
  end subroutine pardiso_analyze


  subroutine pardiso_factor(self, n, ia, ja, a)
    class(PardisoContext), intent(inout) :: self
    integer, intent(in) :: n
    integer, intent(in) :: ia(:), ja(:)
    real(dp), intent(inout) :: a(:)

    integer :: phase, error
    integer :: idum(1)
    real(dp) :: ddum(1)

    if (.not. self%initialized) stop "PardisoContext: call init() first."
    if (.not. self%analyzed)    stop "PardisoContext: call analyze() first."

    phase = 22
    error = 0
    idum(1) = 0
    ddum(1) = 0.0_dp

    call pardiso(self%pt, self%maxfct, self%mnum, self%mtype, phase, n, &
      a, ia, ja, idum, 0, self%iparm, self%msglvl, ddum, ddum, error)

    if (error /= 0) then
      write(*,*) "PARDISO factor error=", error
      stop "PardisoContext: factor failed."
    end if

    self%factorized = .true.
  end subroutine pardiso_factor


  subroutine pardiso_solve(self, n, ia, ja, a, b, x, nrhs)
    class(PardisoContext), intent(inout) :: self
    integer, intent(in) :: n
    integer, intent(in) :: ia(:), ja(:)
    real(dp), intent(inout) :: a(:)
    real(dp), intent(inout) :: b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(in), optional :: nrhs

    integer :: phase, error, nrhs_
    integer :: idum(1)
    real(dp) :: ddum(1)

    if (.not. self%initialized) stop "PardisoContext: call init() first."
    if (.not. self%factorized)  stop "PardisoContext: call factor() first."

    if (size(b) /= n .or. size(x) /= n) stop "PardisoContext: solve size mismatch."

    nrhs_ = self%nrhs_default
    if (present(nrhs)) nrhs_ = nrhs

    phase = 33
    error = 0
    idum(1) = 0
    ddum(1) = 0.0_dp

    call pardiso(self%pt, self%maxfct, self%mnum, self%mtype, phase, n, &
      a, ia, ja, idum, nrhs_, self%iparm, self%msglvl, b, x, error)

    if (error /= 0) then
      write(*,*) "PARDISO solve error=", error
      stop "PardisoContext: solve failed."
    end if
  end subroutine pardiso_solve


  subroutine pardiso_finalize(self, n)
    class(PardisoContext), intent(inout) :: self
    integer, intent(in), optional :: n

    integer :: phase, error, n_
    integer :: idum(1)
    real(dp) :: ddum(1)

    if (.not. self%initialized) then
      if (allocated(self%iparm)) deallocate(self%iparm)
      return
    end if

    ! PARDISO 要一个 n；若外部没给，就用 0（一般也行，但保险起见你最好传）
    n_ = 0
    if (present(n)) n_ = n

    phase = -1
    error = 0
    idum(1) = 0
    ddum(1) = 0.0_dp

    call pardiso(self%pt, self%maxfct, self%mnum, self%mtype, phase, n_, &
      ddum, idum, idum, idum, 0, self%iparm, self%msglvl, ddum, ddum, error)

    if (allocated(self%iparm)) deallocate(self%iparm)

    self%initialized = .false.
    self%analyzed    = .false.
    self%factorized  = .false.
    self%mtype       = 0
  end subroutine pardiso_finalize
  
end module PardisoDriver_mod
