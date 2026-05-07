# add_library(SolverLayer STATIC
#     LinearSolverBase_mod.f90
#     LinearSolverFactory_mod.f90
#     # datri.f90
#     # dasol.f90
#     MKLDirectSolver_mod.f90
#     MKLIterativeSolver_mod.f90
#     SkylineDirectSolver_mod.f90
#     PardisoDriver_mod.f90
# )

# # 依赖 MatrixLayer
# target_link_libraries(SolverLayer PUBLIC MatrixLayer)

# set_target_properties(SolverLayer PROPERTIES
#     Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
# )

# SolverLayer library
add_library(SolverLayer STATIC
    LinearSolverBase_mod.f90
    MKLDirectSolver_mod.f90
    CGSolver_mod.f90
    PCGSolver_mod.f90
    LinearSolverFactory_mod.f90
    SkylineDirectSolver_mod.f90
    PardisoDriver_mod.f90
)

# 依赖 MatrixLayer（保持你原来的）
target_link_libraries(SolverLayer
    PUBLIC
        MatrixLayer
)

# 👉 关键 1：让 ifx 能找到 MKL 的 .mod 文件
target_include_directories(SolverLayer
    PRIVATE
        ${MKL_INCLUDE_DIRS}
)

# 👉 关键 2：链接 MKL（PARDISO 在这里）
target_link_libraries(SolverLayer
    PRIVATE
        MKL::MKL
)

# Fortran module 输出目录（你原来的，保持）
 set_target_properties(SolverLayer PROPERTIES
     Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules
 )
