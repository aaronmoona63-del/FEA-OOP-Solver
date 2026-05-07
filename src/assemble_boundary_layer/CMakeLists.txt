add_library(assemble_boundary_layer STATIC
    ParamIO.f90
    Boundaryconditions.f90
    Element_Utilities.f90
    element_dispatch.f90
    Mesh_sub.f90
    Assembler_coo.f90
    BCPlan_mod.f90
    BCApply_DistributedLoad_mod.f90
    BCApply_DirichletCOO_mod-fast.f90
    BoundaryLayer_mod_new.f90
    BCApply_DirichletCSR_mod.f90
    VTK_Export_mod.f90
    ElementBase_mod.f90
    Element_Hex8_mod.f90
    Assembler_OOP_mod.f90
    MeshIO_mod.f90
    ResultIO_mod.f90
    PostProcessor_mod.f90
    BoundaryCondition_mod.f90
)

target_link_libraries(assemble_boundary_layer
    PUBLIC
    MatrixLayer
    SolverLayer
)

target_include_directories(assemble_boundary_layer
    PUBLIC
    ${CMAKE_BINARY_DIR}/modules
)
