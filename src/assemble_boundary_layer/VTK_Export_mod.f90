module VTK_Export_mod
    implicit none
    private
    public :: export_to_vtk

    integer, parameter :: dp = kind(1.0d0)

contains

    !====================================================================
    ! 导出二维有限元结果到 VTK 格式，供 ParaView 读取
    !====================================================================
    subroutine export_to_vtk(filename, nnode, nelem, coords, conn, x)
        character(len=*), intent(in) :: filename
        integer, intent(in)          :: nnode, nelem
        real(dp), intent(in)         :: coords(:, :) ! 节点坐标 (ndims, nnode)
        integer, intent(in)          :: conn(:, :)   ! 单元连结矩阵 (nodes_per_elem, nelem)
        real(dp), intent(in)         :: x(:)         ! 求解出的全局位移向量
        
        integer :: iunit, i, j, nodes_per_elem
        
        nodes_per_elem = size(conn, 1) ! 自动推断单元节点数（如4代表四边形）
        
        open(newunit=iunit, file=filename, status='replace')
        
        ! 1. VTK 头部要求
        write(iunit, '(A)') '# vtk DataFile Version 3.0'
        write(iunit, '(A)') 'FEM 2D Results'
        write(iunit, '(A)') 'ASCII'
        write(iunit, '(A)') 'DATASET UNSTRUCTURED_GRID'
        
        ! 2. 写入节点坐标
        write(iunit, '(A, I8, A)') 'POINTS ', nnode, ' float'
        do i = 1, nnode
            write(iunit, '(3(E14.6, 1X))') coords(1, i), coords(2, i), 0.0_dp
        end do
        
        ! 3. 写入单元拓扑 (CELLS)
        write(iunit, '(A, I8, I8)') 'CELLS ', nelem, nelem * (nodes_per_elem + 1)
        do i = 1, nelem
            ! 注意：VTK节点索引从0开始，所以要减1
            write(iunit, '(I2, 10(1X, I8))') nodes_per_elem, (conn(j,i)-1, j=1,nodes_per_elem)
        end do
        
        ! 4. 写入单元类型 (CELL_TYPES)
        write(iunit, '(A, I8)') 'CELL_TYPES ', nelem
        do i = 1, nelem
            if (nodes_per_elem == 4) then
                write(iunit, '(I2)') 9 ! 9 = VTK_QUAD (四边形)
            else if (nodes_per_elem == 3) then
                write(iunit, '(I2)') 5 ! 5 = VTK_TRIANGLE (三角形)
            else
                write(iunit, '(I2)') 9 ! 默认
            end if
        end do
        
        ! 5. 写入节点位移场 (POINT_DATA)
        write(iunit, '(A, I8)') 'POINT_DATA ', nnode
        write(iunit, '(A)') 'VECTORS Displacement float'
        do i = 1, nnode
            ! 假设 2D 框架，每个节点有 ux, uy 两个自由度
            write(iunit, '(3(E14.6, 1X))') x(2*i-1), x(2*i), 0.0_dp
        end do
        
        close(iunit)
    end subroutine export_to_vtk

end module VTK_Export_mod