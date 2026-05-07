module MeshIO_mod
  use Types, only : prec
  implicit none
  private
  public :: load_mesh_from_txt

contains

  ! 完全解耦：直接通过参数传出纯数组，不再污染任何全局变量
  subroutine load_mesh_from_txt(filename, n_nodes, n_elements, coords, conn, props)
      character(len=*), intent(in) :: filename
      integer, intent(out) :: n_nodes, n_elements
      real(prec), allocatable, intent(out) :: coords(:,:)
      integer, allocatable, intent(out) :: conn(:,:)
      real(prec), intent(out) :: props(2)

      integer :: file_unit, i, id
      integer :: n1, n2, n3, n4, n5, n6, n7, n8
      real(prec) :: x, y, z, E, nu

      open(newunit=file_unit, file=filename, status='old', action='read')
      read(file_unit, *) n_nodes, n_elements

      allocate(coords(3, n_nodes))
      allocate(conn(8, n_elements))

      do i = 1, n_nodes
          read(file_unit, *) id, x, y, z
          coords(1, i) = x
          coords(2, i) = y
          coords(3, i) = z
      end do

      do i = 1, n_elements
          read(file_unit, *) id, n1, n2, n3, n4, n5, n6, n7, n8
          conn(:, i) = [n1, n2, n3, n4, n5, n6, n7, n8]
      end do

      read(file_unit, *) E, nu
      props(1) = E
      props(2) = nu

      close(file_unit)
      print *, "=> Successfully loaded mesh from: ", trim(filename)
      print *, "   [Nodes]: ", n_nodes, " | [Elements]: ", n_elements
  end subroutine load_mesh_from_txt

end module MeshIO_mod