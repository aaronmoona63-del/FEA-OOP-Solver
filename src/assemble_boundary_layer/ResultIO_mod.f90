module ResultIO_mod
  use Types, only : prec
  implicit none
  private
  public :: export_vtk_hex8

contains

  subroutine export_vtk_hex8(filename, nnode, nelem, coords, conn, u, von_mises)
      character(len=*), intent(in) :: filename
      integer, intent(in)          :: nnode, nelem
      real(prec), intent(in)       :: coords(3, nnode)
      integer, intent(in)          :: conn(8, nelem)
      real(prec), intent(in)       :: u(nnode * 3)
      real(prec), intent(in)       :: von_mises(nelem)  ! 新增的应力数组

      integer :: iunit, id

      open(newunit=iunit, file=filename, status='replace')
      write(iunit, '(A)') '# vtk DataFile Version 3.0'
      write(iunit, '(A)') 'FEM 3D Hex8 Results (OOP Version)'
      write(iunit, '(A)') 'ASCII'
      write(iunit, '(A)') 'DATASET UNSTRUCTURED_GRID'

      write(iunit, '(A, I8, A)') 'POINTS ', nnode, ' float'
      do id = 1, nnode
          write(iunit, '(3(E14.6, 1X))') coords(1, id), coords(2, id), coords(3, id)
      end do

      write(iunit, '(A, I8, I8)') 'CELLS ', nelem, nelem * 9
      do id = 1, nelem
          write(iunit, '(I2, 8(1X, I8))') 8, (conn(1,id)-1), (conn(2,id)-1), (conn(3,id)-1), &
                                             (conn(4,id)-1), (conn(5,id)-1), (conn(6,id)-1), &
                                             (conn(7,id)-1), (conn(8,id)-1)
      end do

      write(iunit, '(A, I8)') 'CELL_TYPES ', nelem
      do id = 1, nelem
          write(iunit, '(I2)') 12  
      end do

      ! 写入节点位移 (POINT_DATA)
      write(iunit, '(A, I8)') 'POINT_DATA ', nnode
      write(iunit, '(A)') 'VECTORS Displacement float'
      do id = 1, nnode
          write(iunit, '(3(E14.6, 1X))') u(3*id-2), u(3*id-1), u(3*id)
      end do
      
      ! 写入单元应力 (CELL_DATA) - 这里是关键！
      write(iunit, '(A, I8)') 'CELL_DATA ', nelem
      write(iunit, '(A)') 'SCALARS Von_Mises float 1'
      write(iunit, '(A)') 'LOOKUP_TABLE default'
      do id = 1, nelem
          write(iunit, '(E14.6)') von_mises(id)
      end do
      
      close(iunit)
      print *, "=> VTK Results (with Stress) exported to: ", trim(filename)
  end subroutine export_vtk_hex8

end module ResultIO_mod