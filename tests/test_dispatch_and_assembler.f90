program test_dispatch_and_assembler
  use Types, only : prec
  use SparseCOO_mod, only : SparseCOO
  use AssemblerCOO_mod, only : assemble_all_elements_coo
  implicit none

  type(SparseCOO) :: Kcoo
  real(prec), allocatable :: rhs(:)
  logical :: fail

  integer, parameter :: n_elements = 1
  ! 修改点：Q4 单元有 4 个节点，每个节点 2 个自由度，所以总共有 8 个未知量
  integer, parameter :: n_unknowns = 8

  ! ---------------------------
  ! Run assembly
  ! ---------------------------
  call assemble_all_elements_coo(n_elements, n_unknowns, Kcoo, rhs, fail)

  if (fail) then
    print *, 'TEST FAILED: assembler returned fail=.true.'
    stop
  end if

  ! ---------------------------
  ! Check RHS
  ! ---------------------------
  print *, 'rhs = ', rhs

  ! ---------------------------
  ! Inspect COO entries
  ! ---------------------------
  print *, 'COO entries (first 20):'
  print *, 'nnz = ', Kcoo%nnz
  print *, 'First few entries:'
  call print_partial_coo(Kcoo, 20)

  print *, 'TEST PASSED'

contains

  subroutine print_partial_coo(A, n_limit)
    type(SparseCOO), intent(in) :: A
    integer, intent(in) :: n_limit
    integer :: k
    
    do k = 1, min(A%nnz, n_limit)
      print '(I4,2X,I4,2X,I4,2X,E12.4)', k, A%row(k), A%col(k), A%val(k)
    end do
  end subroutine

end program test_dispatch_and_assembler
