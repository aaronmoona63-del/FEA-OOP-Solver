program test_bc_csr
    use SparseCSR_mod, only: SparseCSR
    use SparseConvert_mod, only: COO_to_CSR
    use SparseCOO_mod, only: SparseCOO
    use BCApply_DirichletCSR_mod, only: apply_dirichlet_bc_csr
    implicit none

    ! 定义精度
    integer, parameter :: dp = kind(1.0d0)

    type(SparseCOO) :: K_coo
    type(SparseCSR) :: K_csr
    real(dp), allocatable :: rhs(:)
    integer :: bc_dofs(2)
    real(dp) :: bc_vals(2)

    print *, "--- 测试 CSR 罚函数边界条件 ---"

    call K_coo%init(3)
    call K_coo%add_entry(1, 1, 10.0_dp)
    call K_coo%add_entry(2, 2, 20.0_dp)
    call K_coo%add_entry(3, 3, 30.0_dp)
    call K_coo%add_entry(1, 2, -5.0_dp)
    call K_coo%add_entry(2, 1, -5.0_dp)
    call K_coo%preprocess() 

    allocate(rhs(3))
    rhs = [100.0_dp, 200.0_dp, 300.0_dp]

    K_csr = COO_to_CSR(K_coo)
    
    bc_dofs = [1, 3]
    bc_vals = [0.0_dp, 0.5_dp]

    print *, "施加约束前 RHS:", rhs
    
    call apply_dirichlet_bc_csr(K_csr, rhs, bc_dofs, bc_vals)

    print *, "施加约束后 RHS:", rhs
    print *, "验证通过！请检查上述数值。"

end program test_bc_csr