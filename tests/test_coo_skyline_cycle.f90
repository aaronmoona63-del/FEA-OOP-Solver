program TestSkyline
  use SparseCOO_mod
  use SparseSkyline_mod
  use SparseConvert_mod
  implicit none

  type(SparseCOO) :: K, K_after
  type(SparseSkyline) :: Ksky
  logical :: symmetric
  real(8),allocatable :: Kdense(:,:)
  integer:: i
  real(8) :: eps

  !--------------------------------------
  ! Step 1: 构造矩阵 K 的 COO 格式
  !--------------------------------------
  call K%init(5)

  ! K =[1.0000    0.2000    0.4000    0.1000         0
  !     0.3000    2.0000    0.3000         0         0
  !     0.4000    0.3000    3.0000         0    0.5000
  !     0.1000         3         0    4.0000    0.7000
  !          0         3    0.5000    0.7000    5.0000];

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
  print *, "       测试矩阵 K 的 Skyline 闭环转换"
  print *, "====================================================="

  !--------------------------------------
  ! 测试 1：非对称模式（应当通过）
  !--------------------------------------
  symmetric = .false.
  print *, "👉 测试：symmetric = .false."
  allocate(Kdense(5,5))
  Kdense = K%to_dense()
  print *, "K (dense format):"
  do i = 1, 5
    print "(5F10.4)", Kdense(i, :)
  end do

  Ksky = COO_to_Skyline(K, .false.)

  print *, "================ SKYLINE STORAGE ================"

!-----------------------------------------
! 输出 diagonal (ad)
!-----------------------------------------
  print *, "ad (diagonal):"
  print "( *(F10.4,1X) )", Ksky%ad

!-----------------------------------------
! 输出 al（下三角列存储）
!-----------------------------------------
  print *, "al (lower skyline block):"
  if (size(Ksky%al) > 0) then
    print "( *(F10.4,1X) )", Ksky%al
  else
    print *, "<empty>"
  end if

!-----------------------------------------
! 输出 au（上三角列存储：仅在 unsymmetric 时存在）
!-----------------------------------------
  print *, "au (upper skyline block):"
  if (allocated(Ksky%au)) then
    if (size(Ksky%au) > 0) then
      print "( *(F10.4,1X) )", Ksky%au
    else
      print *, "<empty>"
    end if
  else
    print *, "<not allocated>"
  end if
!-----------------------------------------
! 输出 skyline profile jp
!-----------------------------------------
  print *, "jp pointer array (skyline profile):"
  if (size(Ksky%jp) > 0) then
    print "( *(I10,1X) )", Ksky%jp
  else
    print *, "<empty>"
  end if
!-----------------------------------------
! 输出 n 与 nnz
!-----------------------------------------
  print *, "n   =", Ksky%n
  print *, "nnz =", Ksky%nnz

  print *, "================================================="
!--------------------------------------
!     skyline to coo, then print
!--------------------------------------
  K_after = Skyline_to_COO(Ksky, .false.)
  call K_after%print()
  Kdense = K_after%to_dense()
  print *, "K_after (dense format):"
  do i = 1, 5
    print "(5F10.4)", Kdense(i, :)
  end do
  print *, "after prepocess:"
  call K_after%preprocess()
  call K_after%print()

  call K%preprocess()

  eps = 1.0d-10

  print *, "----------------------------------------------------"
  print *, "Checking: original K  <-->  K_after (Skyline loop)"
  print *, "----------------------------------------------------"

  if (K%nnz /= K_after%nnz) then
    print *, "❌ ERROR: nnz mismatch: ", K%nnz, " vs ", K_after%nnz
    stop "COO mismatch: nnz differ"
  end if

  do i = 1, K%nnz
    if (K%row(i) /= K_after%row(i) .or. K%col(i) /= K_after%col(i)) then
      print *, "❌ ERROR: location mismatch at index", i
      print *, "K: ", K%row(i), K%col(i)
      print *, "K_after: ", K_after%row(i), K_after%col(i)
      stop "COO mismatch: row/col differ"
    end if

    if (abs(K%val(i) - K_after%val(i)) > eps) then
      print *, "❌ ERROR: value mismatch at index", i
      print *, "K:", K%val(i), "   K_after:", K_after%val(i)
      stop "COO mismatch: values differ"
    end if
  end do

  print *, "✅ SUCCESS: K and K_after are identical!"

end program TestSkyline
