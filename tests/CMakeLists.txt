enable_testing()

function(add_fortran_test TEST_NAME TEST_FILE)
    add_executable(${TEST_NAME} ${TEST_FILE})
    target_link_libraries(${TEST_NAME} PRIVATE MatrixLayer SolverLayer assemble_boundary_layer)
    target_include_directories(${TEST_NAME} PRIVATE ${CMAKE_BINARY_DIR}/modules)
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endfunction()

add_fortran_test(test_coo_skyline_cycle     test_coo_skyline_cycle.f90)
add_fortran_test(test_coo_crs_cycle         test_coo_crs_cycle.f90)
add_fortran_test(test_coo_add_entry         test_coo_add_entry.f90)
add_fortran_test(test_directSolver          test_directSolver.f90)
add_fortran_test(test_denseSolver           test_denseSolver.f90)
add_fortran_test(test_dispatch_and_assembler test_dispatch_and_assembler.f90)
add_fortran_test(test_multi_element_assembly test_multi_element_assembly.f90)
add_fortran_test(test_dirichlet_bc_coo      test_dirichlet_bc_coo.f90)
add_fortran_test(test_full_solver_mkl       test_full_solver_mkl.f90)
add_fortran_test(test_bc_csr test_bc_csr.f90)
add_fortran_test(test_tension test_tension.f90)
add_fortran_test(test_pure_bending test_pure_bending.f90)
add_fortran_test(test_four_elements_matlab test_four_elements_matlab.f90)

add_executable(test_eccentric_tension test_eccentric_tension.f90)
target_link_libraries(test_eccentric_tension PRIVATE MatrixLayer SolverLayer assemble_boundary_layer)

add_executable(test_two_elements_3d test_two_elements_3d.f90)
target_link_libraries(test_two_elements_3d PRIVATE MatrixLayer SolverLayer assemble_boundary_layer)

add_executable(test_oop_assembly test_oop_assembly.f90)
target_link_libraries(test_oop_assembly PRIVATE MatrixLayer SolverLayer assemble_boundary_layer)

add_executable(test_benchmark test_benchmark.f90)
target_link_libraries(test_benchmark assemble_boundary_layer)

add_executable(test_benchmark_iter test_benchmark_iter.f90)
target_link_libraries(test_benchmark_iter assemble_boundary_layer)

add_executable(test_benchmark_cg test_benchmark_cg.f90)
target_link_libraries(test_benchmark_cg assemble_boundary_layer)

add_executable(test_benchmark_skyline test_benchmark_skyline.f90)
target_link_libraries(test_benchmark_skyline assemble_boundary_layer)

add_executable(test_hole_plate_compare test_hole_plate_compare.f90)
target_link_libraries(test_hole_plate_compare PRIVATE assemble_boundary_layer SolverLayer MatrixLayer)

add_executable(main_hole_plate_analysis main_hole_plate_analysis.f90)
target_link_libraries(main_hole_plate_analysis PRIVATE assemble_boundary_layer SolverLayer MatrixLayer)

add_executable(main_regular_plate_analysis main_regular_plate_analysis.f90)
target_link_libraries(main_regular_plate_analysis PRIVATE assemble_boundary_layer SolverLayer MatrixLayer)

# =============== 核心：新增的大比武程序 ===============
add_executable(main_hole_plate_benchmark main_hole_plate_benchmark.f90)
target_link_libraries(main_hole_plate_benchmark PRIVATE assemble_boundary_layer SolverLayer MatrixLayer)
add_executable(main_hole_medium_benchmark main_hole_medium_benchmark.f90)
target_link_libraries(main_hole_medium_benchmark PRIVATE assemble_boundary_layer SolverLayer MatrixLayer)
add_executable(main_hole_400k main_hole_400k.f90)
target_link_libraries(main_hole_400k PRIVATE assemble_boundary_layer SolverLayer MatrixLayer)
add_executable(main_hole_800k main_hole_800k.f90)
target_link_libraries(main_hole_800k PRIVATE assemble_boundary_layer SolverLayer MatrixLayer)
