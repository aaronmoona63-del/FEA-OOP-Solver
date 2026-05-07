module BCApply_DirichletCSR_mod
    use SparseCSR_mod, only: SparseCSR
    implicit none
    private
    public :: apply_dirichlet_bc_csr

    ! 直接在这里定义双精度参数 dp，与你的矩阵模块保持绝对一致
    integer, parameter :: dp = kind(1.0d0)

contains

    !====================================================================
    ! 采用“罚函数法”对 CSR 格式矩阵施加本质 (Dirichlet) 边界条件
    !====================================================================
    subroutine apply_dirichlet_bc_csr(K_csr, rhs, bc_dofs, bc_vals)
        type(SparseCSR), intent(inout) :: K_csr
        real(dp), intent(inout)        :: rhs(:)
        integer, intent(in)            :: bc_dofs(:)
        real(dp), intent(in)           :: bc_vals(:)
        
        integer :: i, dof, k, start_idx, end_idx
        real(dp), parameter :: PENALTY = 1.0e15_dp ! 极大的惩罚因子
        
        if (size(bc_dofs) /= size(bc_vals)) then
            print *, "Error: bc_dofs and bc_vals size mismatch!"
            return
        end if
        
        do i = 1, size(bc_dofs)
            dof = bc_dofs(i)
            
            start_idx = K_csr%row_ptr(dof)
            end_idx   = K_csr%row_ptr(dof+1) - 1
            
            do k = start_idx, end_idx
                if (K_csr%col_ind(k) == dof) then
                    K_csr%val(k) = K_csr%val(k) * PENALTY
                    rhs(dof) = K_csr%val(k) * bc_vals(i)
                    exit 
                end if
            end do
        end do
        
    end subroutine apply_dirichlet_bc_csr

end module BCApply_DirichletCSR_mod