program TestCSR
  use SparseCOO_mod
  use SparseCSR_mod
  use SparseConvert_mod
  implicit none

  type(SparseCOO) :: K, K_after
  type(SparseCSR) :: Kcsr
  real(8), allocatable :: Kdense(:,:)
  integer :: i
  real(8) :: eps

  !--------------------------------------
  ! Step 1: 构造矩阵 K 的 COO 格式
  !--------------------------------------
  call K%init(5)

  ! Row 1
  call K%add_entry(1,1,1.0d0)
  call K%add_entry(1,2,0.2d0)
  call K%add_entry(1,3,0.4d0)
  call K%add_entry(1,4,0.1d0)

  ! Row 2
  call K%add_entry(2,1,0.3d0)
  call K%add_entry(2,2,2.0d0)
  call K%add_entry(2,3,0.3d0)

  ! Row 3
  call K%add_entry(3,1,0.4d0)
  call K%add_entry(3,2,0.3d0)
  call K%add_entry(3,3,3.0d0)
  call K%add_entry(3,5,0.5d0)

  ! Row 4
  call K%add_entry(4,1,0.1d0)
  call K%add_entry(4,2,3.0d0)
  call K%add_entry(4,4,4.0d0)
  call K%add_entry(4,5,0.7d0)

  ! Row 5
  call K%add_entry(5,2,3.0d0)
  call K%add_entry(5,3,0.5d0)
  call K%add_entry(5,4,0.7d0)
  call K%add_entry(5,5,5.0d0)

  print *, "====================================================="
  print *, "              COO <--> CSR  Conversion Test"
  print *, "====================================================="

  !--------------------------------------
  ! Step 2: 打印 K 的 dense 格式
  !--------------------------------------
  allocate(Kdense(5,5))
  Kdense = K%to_dense()
  print *, "K (dense format):"
  do i = 1, 5
    print "(5F10.4)", Kdense(i, :)
  end do

  !--------------------------------------
  ! Step 3: COO -> CSR
  !--------------------------------------
  print *, "-----------------------------------------------------"
  print *, "Converting COO → CSR"
  print *, "-----------------------------------------------------"
  Kcsr = COO_to_CSR(K)

  print *, "row_ptr:"
  print "( *(I6,1X) )", Kcsr%row_ptr

  print *, "col_ind:"
  print "( 5(I6,1X) )", Kcsr%col_ind

  print *, "val:"
  print "( 5(F10.4,1X) )", Kcsr%val

  !--------------------------------------
  ! Step 4: CSR -> COO
  !--------------------------------------
  print *, "-----------------------------------------------------"
  print *, "Converting CSR → COO"
  print *, "-----------------------------------------------------"
  K_after = CSR_to_COO(Kcsr)

  print *, "K_after (COO format):"
  call K_after%print()

  !--------------------------------------
  ! Step 5: 打印 dense 格式
  !--------------------------------------
  Kdense = K_after%to_dense()
  print *, "K_after (dense format):"
  do i = 1, 5
    print "(5F10.4)", Kdense(i, :)
  end do

  !--------------------------------------
  ! Step 6: COO preprocess（排序 + 去重）
  !--------------------------------------
  call K%preprocess()
  !call K_after%preprocess()

  !--------------------------------------
  ! Step 7: 比较 K 和 K_after（必须一致）
  !--------------------------------------
  eps = 1.0d-10

  print *, "----------------------------------------------------"
  print *, "Checking: original K  <-->  K_after (CSR loop)"
  print *, "----------------------------------------------------"

  if (K%nnz /= K_after%nnz) then
    print *, "❌ ERROR: nnz mismatch: ", K%nnz, " vs ", K_after%nnz
    stop "CSR mismatch: nnz differ"
  end if

  do i = 1, K%nnz
    if (K%row(i) /= K_after%row(i) .or. K%col(i) /= K_after%col(i)) then
      print *, "❌ ERROR: location mismatch at index", i
      print *, "K: ", K%row(i), K%col(i)
      print *, "K_after: ", K_after%row(i), K_after%col(i)
      stop "CSR mismatch: row/col differ"
    end if

    if (abs(K%val(i) - K_after%val(i)) > eps) then
      print *, "❌ ERROR: value mismatch at index", i
      print *, "K:", K%val(i), "   K_after:", K_after%val(i)
      stop "CSR mismatch: values differ"
    end if
  end do

  print *, "✅ SUCCESS: K and K_after match after CSR conversion!"
  print *, "====================================================="

end program TestCSR
