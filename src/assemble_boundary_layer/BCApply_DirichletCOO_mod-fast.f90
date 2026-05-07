module BCApply_DirichletCOO_mod
  !=============================================================
  ! Apply Dirichlet BCs on SparseCOO (high performance, correct RHS)
  !
  ! Enforce u(eq) = ubar for each eq in dbc_eq.
  !
  ! Correct elimination:
  !   For any row i != eq:
  !     rhs(i) := rhs(i) - K(i,eq) * ubar
  !   Remove all entries with row==eq OR col==eq
  !   Add (eq,eq)=1 and set rhs(eq)=ubar
  !
  ! Implementation (fast):
  !   - Build a boolean marker isD(1:n) and a value vector ubar(1:n)
  !   - One pass over nnz:
  !       * RHS correction for entries whose col is Dirichlet
  !       * Compact in-place keeping only entries with row/col NOT Dirichlet
  !   - Reset nnz to kept count, then append diagonal entries for Dirichlet eqs
  !
  ! Complexity:
  !   O(nnz + nbc) time, O(n) extra memory.
  !=============================================================
  use Types,        only : prec
  !use ParamIO,      only : IOW
  use SparseCOO_mod, only : SparseCOO
  implicit none
  private
  public :: apply_dirichlet_coo

  integer:: IOW=6

contains

  subroutine apply_dirichlet_coo(Kcoo, rhs, dbc_eq, dbc_val)
    type(SparseCOO), intent(inout) :: Kcoo
    real(prec),     intent(inout) :: rhs(:)
    integer,        intent(in)    :: dbc_eq(:)
    real(prec),     intent(in)    :: dbc_val(:)

    integer :: n, nnz_old, nbc
    integer :: k, w, i, j, p, eq
    logical, allocatable :: isD(:)
    real(prec), allocatable :: ubar(:)

    if (size(dbc_eq) /= size(dbc_val)) stop "apply_dirichlet_coo: size mismatch"

    n = size(rhs)
    if (Kcoo%n /= n) then
      write(IOW,*) "apply_dirichlet_coo: size(rhs) != Kcoo%n"
      write(IOW,*) "  rhs size = ", n, "  Kcoo%n = ", Kcoo%n
      stop
    end if

    nbc = size(dbc_eq)
    if (nbc == 0) return

    allocate(isD(n));  isD  = .false.
    allocate(ubar(n)); ubar = 0.0_prec

    ! Mark Dirichlet dofs and store prescribed values
    do p = 1, nbc
      eq = dbc_eq(p)
      if (eq < 1 .or. eq > n) then
        write(IOW,*) "apply_dirichlet_coo: eq out of range:", eq
        stop
      end if
      isD(eq)  = .true.
      ubar(eq) = dbc_val(p)
    end do

    !-----------------------------------------------------------
    ! Pass 1: RHS correction + in-place compaction
    !   - Read pointer: k = 1..nnz_old
    !   - Write pointer: w = number of kept entries
    !-----------------------------------------------------------
    nnz_old = Kcoo%nnz
    w = 0

    do k = 1, nnz_old
      i = Kcoo%row(k)
      j = Kcoo%col(k)

      ! RHS correction: if column is Dirichlet and row is NOT Dirichlet
      if (isD(j) .and. .not. isD(i)) then
        rhs(i) = rhs(i) - Kcoo%val(k) * ubar(j)
      end if

      ! Keep only non-Dirichlet rows and columns
      if (.not. isD(i) .and. .not. isD(j)) then
        w = w + 1
        if (w /= k) then
          Kcoo%row(w) = Kcoo%row(k)
          Kcoo%col(w) = Kcoo%col(k)
          Kcoo%val(w) = Kcoo%val(k)
        end if
      end if
    end do

    Kcoo%nnz = w

    !-----------------------------------------------------------
    ! Pass 2: Add diagonal (eq,eq)=1 and set rhs(eq)=ubar(eq)
    !   - Use your existing add_entry (handles grow)
    !-----------------------------------------------------------
    do p = 1, nbc
      eq = dbc_eq(p)
      call Kcoo%add_entry(eq, eq, 1.0_prec)
      rhs(eq) = ubar(eq)
    end do

    deallocate(isD, ubar)
  end subroutine apply_dirichlet_coo

end module BCApply_DirichletCOO_mod
