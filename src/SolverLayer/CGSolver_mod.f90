module CGSolver_mod
  use LinearSolverBase_mod, only : LinearSolver, dp
  use SparseCSR_mod,        only : SparseCSR
  implicit none
  private
  public :: CGSolver

  type, extends(LinearSolver) :: CGSolver
      type(SparseCSR), pointer :: A_csr => null()
  contains
      procedure :: attach_matrix => cg_attach_matrix
      procedure :: analyze       => cg_analyze
      procedure :: factor        => cg_factor
      procedure :: solve         => cg_solve
      procedure :: free          => cg_free
  end type CGSolver

contains

  subroutine cg_attach_matrix(self, A)
      class(CGSolver), intent(inout) :: self
      class(*), intent(in), target    :: A
      select type(pA => A)
      type is (SparseCSR)
          self%A_csr => pA
      class default
          stop "CGSolver: attach_matrix expects SparseCSR."
      end select
  end subroutine cg_attach_matrix

  subroutine cg_analyze(self)
      class(CGSolver), intent(inout) :: self
      self%analyzed = .true.
  end subroutine cg_analyze

  subroutine cg_factor(self)
      class(CGSolver), intent(inout) :: self
      ! 标准 CG 无需预处理，直接跳过
      self%factorized = .true.
  end subroutine cg_factor

  subroutine cg_solve(self, b, x)
      class(CGSolver), intent(inout) :: self
      real(dp), intent(inout)         :: b(:)
      real(dp), intent(out)           :: x(:)
      
      integer :: n, max_iter, iter, i, k
      real(dp) :: tol, b_norm, r_norm, alpha, beta, rho, rho_new, p_dot_q
      real(dp), allocatable :: r(:), p(:), q(:)
      
      if (.not. associated(self%A_csr)) stop "CGSolver: A_csr not associated!"
      n = self%A_csr%n
      allocate(r(n), p(n), q(n))
      
      max_iter = self%opts%max_iter
      if (max_iter <= 0) max_iter = 10000
      tol = self%opts%tol
      if (tol <= 0.0_dp) tol = 1.0e-6_dp
      
      x = 0.0_dp
      r = b  
      b_norm = sqrt(sum(b**2))
      if (b_norm < 1.0e-14_dp) b_norm = 1.0_dp
      
      p = r
      rho = dot_product(r, r)
      
      do iter = 1, max_iter
          ! 稀疏矩阵向量乘法 (SpMV)
          q = 0.0_dp
          do i = 1, n
              do k = self%A_csr%row_ptr(i), self%A_csr%row_ptr(i+1)-1
                  q(i) = q(i) + self%A_csr%val(k) * p(self%A_csr%col_ind(k))
              end do
          end do
          
          p_dot_q = dot_product(p, q)
          alpha = rho / p_dot_q
          x = x + alpha * p
          r = r - alpha * q
          r_norm = sqrt(sum(r**2))
          
          if (r_norm / b_norm < tol) then
              print *, "  => [标准 CG] 成功收敛! 迭代次数: ", iter
              exit
          end if
          
          rho_new = dot_product(r, r)
          beta = rho_new / rho
          p = r + beta * p
          rho = rho_new
      end do
      
      if (iter > max_iter) print *, "  => [标准 CG 警告] 达到最大迭代次数未收敛！"
      deallocate(r, p, q)
  end subroutine cg_solve

  subroutine cg_free(self)
      class(CGSolver), intent(inout) :: self
      nullify(self%A_csr)
      self%analyzed = .false.
      self%factorized = .false.
  end subroutine cg_free
end module CGSolver_mod
