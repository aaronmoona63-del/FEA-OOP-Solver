program test_multi_element_assembly
  use Types, only : prec
  use SparseCOO_mod
  use AssemblerCOO_mod
  implicit none
  type(SparseCOO) :: Kcoo
  real(prec), allocatable :: rhs(:)
  logical :: fail
  integer :: i
  print *, "--- FORCED REBUILD: Running Multi-Element Test ---"
  call assemble_all_elements_coo(2, 12, Kcoo, rhs, fail)
  if (fail) stop "assembly failed"
  if (Kcoo%nnz /= 112) then
    print *, "Matrix NNZ mismatch! Expected 112, got ", Kcoo%nnz
    stop 1
  end if
  if (any(abs(rhs) > 1e-10)) then
    print *, "RHS should be zero, but is not!"
    stop 1
  end if
  print *, "test_multi_element_assembly PASSED"
end program
