module SkylineLUSolverMod
  use Types
  use MatrixBase_mod
  use SolverBase_mod
  use SparseSkyline_mod
  implicit none
  private
  public :: SkylineLUSolver

  !===========================================================
  ! Skyline LU Solver: A ≈ L * U （无主元）
  ! 说明：
  !  - factorize: 直接在 A(SkylineMatrix) 上做原地 LU 分解
  !  - solve    : 做一次前代 + 对角缩放 + 回代
  !===========================================================
  type, extends(DirectSolverType) :: SkylineLUSolver
  logical :: is_factorized = .false.
contains
  procedure :: factorize => skyline_factorize
  procedure :: solve     => skyline_solve
  end type SkylineLUSolver

contains

  !===============================================================
  !  dredu —— 复制自你给的版本，用于更新对角线
  !===============================================================
  subroutine dredu(alow, aupp, diag, jh, ifl, dj)
    use Types
    implicit none
    integer,  intent(in)    :: jh
    logical,  intent(in)    :: ifl
    real(prec), intent(inout) :: alow(jh)
    real(prec), intent(inout) :: aupp(jh)
    real(prec), intent(in)    :: diag(jh)
    real(prec), intent(inout) :: dj

    real(prec) :: ud
    integer    :: k

    ! Reduce diagonal dj = dj - Σ alow(k) * aupp(k) * diag(k)
    do k = 1, jh
      ud = aupp(k)*diag(k)
      dj = dj - alow(k)*ud
      aupp(k) = ud
    end do

    ! 非对称情形：还要更新 alow
    if (ifl) then
      alow(1:jh) = alow(1:jh)*diag(1:jh)
    end if
  end subroutine dredu

  !===============================================================
  !  Skyline LU factorization（原地修改 A）
  !===============================================================
  subroutine skyline_factorize(self, A)
    class(SkylineLUSolver), intent(inout) :: self
    class(MatrixType),      intent(inout) :: A

    class(SkylineMatrix), pointer :: Ask
    integer :: n
    integer :: j, jd, jr, jh, i
    integer :: is, ie
    integer :: id, ih, jrh, idh
    real(prec) :: dd, daval, dimn, dimx, dfig
    logical :: unsym
    integer :: IOW=6

    !---- 类型选择：必须是 SkylineMatrix
    select type(A)
      type is (SkylineMatrix)
      Ask => A
      class default
      stop "SkylineLUSolver: A must be SkylineMatrix"
    end select

    ! 使用基类里的 nrows（方阵：nrows = ncols）
    n     = Ask%nrows
    unsym = Ask%unsymmetric   ! 或者：unsym = .not. Ask%is_symmetric

    !     --- Set initial values for conditioning check
    dimx = 0.0D0
    dfig = 0.0D0
    !  do i = 1, neq
    !    dimn = max(dimn, dabs(diag(i)))
    !  end do
    dimn = maxval(dabs(Ask%diag))

    !     --- Loop through columns to perform triangular decomposition
    jd = 1
    do j = 1, n
      jr = jd + 1
      jd = Ask%jpoin(j)
      jh = jd - jr
      if ( jh > 0 ) then
        is = j - jh
        ie = j - 1
        !     ---  If diagonal is zero compute a norm for singularity test
        if ( Ask%diag(j) == 0.0D0 ) daval =  sum(dabs(Ask%aupp(jr:jr+jh)))
        do i = is, ie
          jr = jr + 1
          id = Ask%jpoin(i)
          ih = min(id - Ask%jpoin(i - 1), i - is + 1)
          if ( ih > 0 ) then
            jrh = jr - ih
            idh = id - ih + 1
            Ask%aupp(jr) = Ask%aupp(jr) - dot_product(Ask%aupp(jrh:jrh+ih-1), Ask%alow(idh:idh+ih-1))
            if ( unsym ) Ask%alow(jr) = Ask%alow(jr) - dot_product(Ask%alow(jrh:jrh+ih-1), Ask%aupp(idh:idh+ih-1))
          end if
        end do
      end if
      !     ---   Reduce the diagonal
      if ( jh >= 0 ) then
        dd = Ask%diag(j)
        jr = jd - jh
        jrh = j - jh - 1
        call dredu(Ask%alow(jr:jr+jh), Ask%aupp(jr:jr+jh+1), Ask%diag(jrh:jrh+jh+1), jh + 1, unsym, Ask%diag(j) )
        !     ---   Check for conditioning errors and print warnings
        if ( dabs(Ask%diag(j)) < 0.5D-07*dabs(dd) ) write (IOW, 99001) j
        if ( Ask%diag(j) == 0.0D0 ) then
          if ( j/=n ) write (IOW, 99003) j
        end if
        if ( dd == 0.0D0 .and. jh > 0 ) then
          if ( dabs(Ask%diag(j)) < 0.5D-07*daval ) write (IOW, 99004) j
        end if
      end if

      !     ---    Store reciprocal of diagonal, compute condition checks
      if ( Ask%diag(j) /= 0.0D0 ) then
        dimx = dmax1(dimx, dabs(Ask%diag(j)))
        dimn = dmin1(dimn, dabs(Ask%diag(j)))
        dfig = dmax1(dfig, dabs(dd/Ask%diag(j)))
        Ask%diag(j) = 1.D0/Ask%diag(j)
      else
        Ask%diag(j) = 1.D0/(0.5D-07*dimn)
      end if
      self%is_factorized=.true.

    end do
  
99001 format (/' **** DIRECT SOLVER WARNING 1 **** '/  &
      '  Loss of at least 7 digits in reducing diagonal', ' of equation ', i5)
99003 format (/' **** DIRECT SOLVER WARNING 2 **** '/  &
      '  Reduced diagonal is zero for equation ', i5)
99004 format (/' **** DIRECT SOLVER WARNING 3 **** '/  &
      '  Rank failure for zero unreduced diagonal in ', 'equation', i5)
99005 format ( // ' Direct Solver has completed LU decomposition '/  &
      '    Conditioning information: '/  &
      '      Max diagonal in reduced matrix:   ',  &
      e11.4/'      Min diagonal in reduced matrix:   ',  &
      e11.4/'      Ratio:                            ',  &
      e11.4/'      Maximum no. diagonal digits lost: ', i3)  
  end subroutine skyline_factorize

  !===============================================================
  ! skyline_solve: 前代 + 对角缩放 + 回代
  !===============================================================
!===============================================================
! skyline_solve: FEAP-style forward + diag + backward
!===============================================================
subroutine skyline_solve(self, A, rhs, x)
  class(SkylineLUSolver), intent(inout) :: self
  class(MatrixType),      intent(in)    :: A
  real(prec),             intent(in)    :: rhs(:)
  real(prec),             intent(out)   :: x(:)

  class(SkylineMatrix), pointer :: Ask

  integer :: n
  integer :: is, j, jr, jh
  real(prec), parameter :: zero = 0.0_prec
  logical :: foundzerorhs

  if (.not. self%is_factorized) stop "SkylineLUSolver: factorize() must be called first"

  select type(A)
    type is (SkylineMatrix)
       Ask => A
    class default
       stop "SkylineLUSolver: A must be SkylineMatrix"
  end select

  n = Ask%nrows

  if (size(rhs) /= n .or. size(x) /= n) stop "skyline_solve: size mismatch"

  ! 拷贝 RHS
  x = rhs

  !===============================================================
  ! 1) Find first non-zero RHS entry: FEAP behavior
  !===============================================================
  foundzerorhs = .true.
  do is = 1, n
    if (x(is) /= zero) then
      foundzerorhs = .false.
      exit
    end if
  end do

  ! if RHS all zero, return zero solution
  if (foundzerorhs) then
    x = zero
    return
  end if

  !===============================================================
  ! 2) Forward substitution: L * y = rhs
  !    注意 FEAP 是从 is+1 开始！
  !===============================================================
  if (is < n) then
    do j = is+1, n
      jr = Ask%jpoin(j-1)
      jh = Ask%jpoin(j) - jr
      if (jh > 0) then
        x(j) = x(j) - dot_product( Ask%alow(jr+1:jr+jh), x(j-jh:j-1) )
      end if
    end do
  end if

  !===============================================================
  ! 3) Diagonal scaling: y = D^{-1} * y
  !===============================================================
  x(is:n) = x(is:n) * Ask%diag(is:n)

  !===============================================================
  ! 4) Backward substitution: U * x = y
  !===============================================================
  if (n > 1) then
    do j = n, 2, -1
      jr = Ask%jpoin(j-1)
      jh = Ask%jpoin(j) - jr
      if (jh > 0) then
        x(j-jh:j-1) = x(j-jh:j-1) - x(j) * Ask%aupp(jr+1:jr+jh)
      end if
    end do
  end if

end subroutine skyline_solve


end module SkylineLUSolverMod
