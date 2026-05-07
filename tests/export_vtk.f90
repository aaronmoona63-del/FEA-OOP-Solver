! ======================================================================
! 将 Hex8 网格和位移结果输出为 ParaView 可读的 VTK 文件
! ======================================================================
SUBROUTINE ExportToVTK(filename, numNodes, numElems, nodes, elements, U)
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: filename
    INTEGER, INTENT(in)          :: numNodes, numElems
    ! 假设 nodes 是 (numNodes, 3)，elements 是 (numElems, 8)，U 是位移 (numNodes, 3)
    REAL*8, INTENT(in)           :: nodes(numNodes, 3)
    INTEGER, INTENT(in)          :: elements(numElems, 8)
    REAL*8, INTENT(in)           :: U(numNodes, 3)
    
    INTEGER :: i, fileUnit
    
    fileUnit = 99
    OPEN(unit=fileUnit, file=filename, status='replace', action='write')
    
    ! 1. VTK 文件头 (Header)
    WRITE(fileUnit, '(A)') '# vtk DataFile Version 3.0'
    WRITE(fileUnit, '(A)') 'FEM 3D Hex8 Results (OOP Version)'
    WRITE(fileUnit, '(A)') 'ASCII'
    WRITE(fileUnit, '(A)') 'DATASET UNSTRUCTURED_GRID'
    
    ! 2. 写入节点坐标 (POINTS)
    WRITE(fileUnit, '(A, I0, A)') 'POINTS ', numNodes, ' float'
    DO i = 1, numNodes
        WRITE(fileUnit, '(3(E15.6, 1X))') nodes(i, 1), nodes(i, 2), nodes(i, 3)
    END DO
    
    ! 3. 写入单元连接关系 (CELLS)
    ! 每个 Hex8 单元由 8 个节点组成，在文件里第一个数字是表示节点数量的 8。
    ! 所以对于 n 个单元，共有 n * 9 个整数需要写入。
    WRITE(fileUnit, '(A, I0, 1X, I0)') 'CELLS ', numElems, numElems * 9
    DO i = 1, numElems
        ! ⭐ 极其重要的坑：VTK 的节点编号是从 0 开始的！
        ! 你的 Fortran 代码里的编号通常是从 1 开始的，所以这里必须全部减 1！
        WRITE(fileUnit, '(I2, 8(1X, I0))') 8, elements(i, 1)-1, elements(i, 2)-1, &
                                              elements(i, 3)-1, elements(i, 4)-1, &
                                              elements(i, 5)-1, elements(i, 6)-1, &
                                              elements(i, 7)-1, elements(i, 8)-1
    END DO
    
    ! 4. 写入单元类型 (CELL_TYPES)
    ! VTK_HEXAHEDRON 对应的数字编号是 12
    WRITE(fileUnit, '(A, I0)') 'CELL_TYPES ', numElems
    DO i = 1, numElems
        WRITE(fileUnit, '(I2)') 12
    END DO
    
    ! 5. 写入节点位移数据 (POINT_DATA)
    ! 这里将每个节点的 X, Y, Z 位移作为一个矢量写入
    WRITE(fileUnit, '(A, I0)') 'POINT_DATA ', numNodes
    WRITE(fileUnit, '(A)') 'VECTORS Displacement float'
    DO i = 1, numNodes
        WRITE(fileUnit, '(3(E15.6, 1X))') U(i, 1), U(i, 2), U(i, 3)
    END DO
    
    CLOSE(fileUnit)
    PRINT *, "✅ 成功生成 VTK 结果文件: ", TRIM(filename)

END SUBROUTINE ExportToVTK