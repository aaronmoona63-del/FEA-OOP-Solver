# MatrixLayer - 稀疏矩阵基础设施
add_library(MatrixLayer STATIC
    Types.f90
    MatrixBase_mod.f90
    MatrixDense_mod.f90
    MatrixUtilities_mod.f90
    SparseMatrixBase_mod.f90
    SparseCOO_mod.f90
    SparseCSR_mod.f90
    SparseSkyline_mod.f90
    SparseConvert_mod.f90
)

# 指定 .mod 文件输出目录
 set_target_properties(MatrixLayer PROPERTIES
     Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
 )
