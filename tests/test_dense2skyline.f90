program test_skyline
  use MatrixDense
  implicit none

  type(DenseMatrix) :: K1,K2 
  integer :: n
  double precision, allocatable :: al(:), au(:), ad(:)
  integer, allocatable :: jp(:)

  n = 5
  call K1%init(n, .true.)   !  symmetric
  call K2%init(n, .false.)   ! not symmetric

  ! Fill matrix (example)
  K1%A = reshape( (/ &
  1.0d0, 0.2d0, 0.4d0, 0.1d0, 0.0d0, &
  0.2d0, 2.0d0, 0.3d0, 0.0d0, 0.0d0, &
  0.4d0, 0.3d0, 3.0d0, 0.0d0, 0.5d0, &
  0.1d0, 0.0d0, 0.0d0, 4.0d0, 0.7d0, &
  0.0d0, 0.0d0, 0.5d0, 0.7d0, 5.0d0  /), [n,n] )


  ! Fill matrix (example)
  K2%A = reshape( (/ &
  1.0d0, 0.3d0, 0.4d0, 0.1d0, 0.0d0, &   ! 第 1 列
  0.2d0, 2.0d0, 0.3d0, 3.0d0, 3.0d0, &   ! 第 2 列
  0.4d0, 0.3d0, 3.0d0, 0.0d0, 0.5d0, &   ! 第 3 列
  0.1d0, 0.0d0, 0.0d0, 4.0d0, 0.7d0, &   ! 第 4 列
  0.0d0, 0.0d0, 0.5d0, 0.7d0, 5.0d0  /), [n,n] )

  print *, "Dense matrix K1:"
  call K1%print()

  ! Convert to skyline
  call K1%to_skyline(al, au, ad, jp)

  print *, "ad = ", ad
  print *, "jp = ", jp
  print *, "al = ", al
  print *, "au = ", au


  print *, "Dense matrix K2:"
  call K2%print()

  ! Convert to skyline
  call K2%to_skyline(al, au, ad, jp)

  print *, "ad = ", ad
  print *, "jp = ", jp
  print *, "al = ", al
  print *, "au = ", au

end program test_skyline
