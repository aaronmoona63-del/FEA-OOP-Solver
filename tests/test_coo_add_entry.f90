program TestAddEntry
  use SparseCOO_mod
  implicit none

  type(SparseCOO) :: A
  real(8), allocatable :: D(:,:)
  integer :: i,j

  !--------------------------------------
  ! 1. 初始化，矩阵大小 n = 5
  !--------------------------------------
  call A%init(5)

  print *, "Initial capacity =", A%capacity
  print *, "Initial nnz      =", A%nnz
  print *, "------------------------------------------------"

  !--------------------------------------
  ! 2. 添加一些三元组
  !--------------------------------------
  print *, "Adding entries..."
  call A%add_entry(1,1,1.0d0)
  call A%add_entry(1,2,2.0d0)
  call A%add_entry(3,3,3.0d0)
  call A%add_entry(5,4,4.0d0)
  call A%add_entry(2,1,5.0d0)

  print *, "After 5 entries:"
  print *, "nnz      =", A%nnz
  print *, "capacity =", A%capacity
  print *, "------------------------------------------------"

  !--------------------------------------
  ! 3. 添加更多条目以触发自动扩容
  !--------------------------------------
  print *, "Adding many entries to test auto-grow..."

  do i = 1, 200
    call A%add_entry(mod(i,5)+1, mod(i,5)+1, 0.1d0 * i)
  end do

  print *, "After adding 200 entries:"
  print *, "nnz      =", A%nnz
  print *, "capacity =", A%capacity
  print *, "------------------------------------------------"

  !--------------------------------------
  ! 4. 打印 COO 内容（前 20 个）
  !--------------------------------------
  print *, "Print first 20 entries:"
  A%nnz = min(A%nnz, 20)    ! 防止打印太多
  call A%print()

  ! 恢复 nnz（如果你不想修改 nnz，可以复制数组打印）
  A%nnz = 200 + 5

  !--------------------------------------
  ! 5. COO → Dense 验证是否正确
  !--------------------------------------
  D = A%to_dense()

  print *, "Dense matrix:"

! 打印列号
  write(*,'(6X, *(I10,1x))') (j, j = 1, size(D,2))

! 打印矩阵内容（带行号）
  do i = 1, size(D,1)
    write(*,'(I4,2X, *(F10.4,1X))') i, D(i,:)
  end do


end program TestAddEntry
