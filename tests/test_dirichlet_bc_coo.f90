program test_dirichlet_bc_coo
  use Types,              only : prec
  use SparseCOO_mod
  use BCApply_DirichletCOO_mod, only : apply_dirichlet_coo
  implicit none

  type(SparseCOO) :: Kcoo
  real(prec), allocatable :: rhs(:)
  real(prec), allocatable :: Kdense(:,:)

  integer :: n
  integer, allocatable :: dbc_eq(:)
  real(prec), allocatable :: dbc_val(:)

  !------------------------------------------------------------
  ! Problem size
  !------------------------------------------------------------
  n = 3
  allocate(rhs(n))
  rhs = 0.0_prec

  call Kcoo%init(n)

  !------------------------------------------------------------
  ! Assemble 3x3 tridiagonal matrix:
  !
  ! [ 2 -1  0 ]
  ! [ -1 2 -1 ]
  ! [ 0 -1  2 ]
  !------------------------------------------------------------
  call Kcoo%add_entry(1,1,  2.0_prec)
  call Kcoo%add_entry(1,2, -1.0_prec)

  call Kcoo%add_entry(2,1, -1.0_prec)
  call Kcoo%add_entry(2,2,  2.0_prec)
  call Kcoo%add_entry(2,3, -1.0_prec)

  call Kcoo%add_entry(3,2, -1.0_prec)
  call Kcoo%add_entry(3,3,  2.0_prec)

  !------------------------------------------------------------
  ! Dirichlet BC: u(1) = 5
  !------------------------------------------------------------
  allocate(dbc_eq(1), dbc_val(1))
  dbc_eq(1)  = 1
  dbc_val(1) = 5.0_prec

  call apply_dirichlet_coo(Kcoo, rhs, dbc_eq, dbc_val)

  !------------------------------------------------------------
  ! Preprocess (merge/sort)
  !------------------------------------------------------------
  call Kcoo%preprocess()

  !------------------------------------------------------------
  ! Check RHS correctness
  !------------------------------------------------------------
  print *, "RHS after Dirichlet:"
  print *, rhs

  if (abs(rhs(1) - 5.0_prec) > 1.0e-12_prec) then
    print *, "ERROR: rhs(1) wrong"
    stop 1
  end if

  if (abs(rhs(2) - 5.0_prec) > 1.0e-12_prec) then
    print *, "ERROR: rhs(2) should be +5 due to elimination"
    stop 1
  end if

  if (abs(rhs(3)) > 1.0e-12_prec) then
    print *, "ERROR: rhs(3) should remain 0"
    stop 1
  end if

  !------------------------------------------------------------
  ! Inspect resulting matrix
  !------------------------------------------------------------
  print *
  print *, "COO matrix after Dirichlet (row, col, val):"
  call Kcoo%print()

  Kdense = Kcoo%to_dense()

  print *
  print *, "Dense matrix after Dirichlet:"
  print '(3F10.3)', Kdense(1,:)
  print '(3F10.3)', Kdense(2,:)
  print '(3F10.3)', Kdense(3,:)

  ! Expected:
  ! [ 1  0  0 ]
  ! [ 0  2 -1 ]
  ! [ 0 -1  2 ]

  if (abs(Kdense(1,1) - 1.0_prec) > 1.0e-12_prec) stop "ERROR: K(1,1) wrong"
  if (abs(Kdense(1,2)) > 1.0e-12_prec) stop "ERROR: K(1,2) not zero"
  if (abs(Kdense(2,1)) > 1.0e-12_prec) stop "ERROR: K(2,1) not zero"

  print *
  print *, "test_dirichlet_bc_coo_fast PASSED"

end program test_dirichlet_bc_coo