module PCGSolver_mod
  use LinearSolverBase_mod, only : LinearSolver, dp
  use SparseCSR_mod,        only : SparseCSR
  implicit none
  private
  public :: PCGSolver

  type, extends(LinearSolver) :: PCGSolver
      ! 👉 新增：用于存放绑定进来的稀疏矩阵指针
      type(SparseCSR), pointer :: A_csr => null()
      ! 用于存放对角线预处理数据
      real(dp), allocatable :: inv_diag(:)
  contains
      ! 严格对应基类的五大接口！
      procedure :: attach_matrix => pcg_attach_matrix
      procedure :: analyze       => pcg_analyze
      procedure :: factor        => pcg_factor
      procedure :: solve         => pcg_solve
      procedure :: free          => pcg_free
  end type PCGSolver

contains

  ! ========================================================
  ! 1. 绑定矩阵 (解决 error #8322 和 error #8247)
  ! ========================================================
  subroutine pcg_attach_matrix(self, A)
      class(PCGSolver), intent(inout) :: self
      class(*), intent(in), target    :: A
      
      select type(pA => A)
      type is (SparseCSR)
          self%A_csr => pA
      class default
          stop "PCGSolver: attach_matrix expects SparseCSR."
      end select
  end subroutine pcg_attach_matrix

  ! ========================================================
  ! 2. 分析阶段 (PCG不需要额外分析结构，直接标记即可)
  ! ========================================================
  subroutine pcg_analyze(self)
      class(PCGSolver), intent(inout) :: self
      self%analyzed = .true.
  end subroutine pcg_analyze

  ! ========================================================
  ! 3. 预处理提取阶段 (完美化解大数惩罚法)
  ! ========================================================
  subroutine pcg_factor(self)
      class(PCGSolver), intent(inout) :: self
      integer :: i, k, n
      
      if (.not. associated(self%A_csr)) stop "PCGSolver: A_csr not associated!"

      n = self%A_csr%n  ! 修正 error #6460 (你的类里叫 n 不是 n_rows)
      
      if (allocated(self%inv_diag)) deallocate(self%inv_diag)
      allocate(self%inv_diag(n))
      self%inv_diag = 1.0_dp

      ! 提取对角线并求倒数
      do i = 1, n
          do k = self%A_csr%row_ptr(i), self%A_csr%row_ptr(i+1)-1
              if (self%A_csr%col_ind(k) == i) then
                  if (abs(self%A_csr%val(k)) > 1.0e-20_dp) then
                      self%inv_diag(i) = 1.0_dp / self%A_csr%val(k)
                  end if
                  exit
              end if
          end do
      end do
      self%factorized = .true.
  end subroutine pcg_factor

  ! ========================================================
  ! 4. 核心求解阶段 (严格按照基类的 intent 定义参数)
  ! ========================================================
  subroutine pcg_solve(self, b, x)
      class(PCGSolver), intent(inout) :: self
      real(dp), intent(inout)         :: b(:)
      real(dp), intent(out)           :: x(:)
      
      integer :: n, max_iter, iter, i, k
      real(dp) :: tol, b_norm, r_norm, alpha, beta, rho, rho_new, p_dot_q, current_err
      real(dp), allocatable :: r(:), z(:), p(:), q(:)
      
      if (.not. associated(self%A_csr)) stop "PCGSolver: A_csr not associated!"
      n = self%A_csr%n
      allocate(r(n), z(n), p(n), q(n))
      
      ! 动态读取基类中的设定选项
      max_iter = self%opts%max_iter
      if (max_iter <= 0) max_iter = 10000
      tol = self%opts%tol
      if (tol <= 0.0_dp) tol = 1.0e-6_dp
      
      x = 0.0_dp
      r = b  ! 初始残差
      b_norm = sqrt(sum(b**2))
      if (b_norm < 1.0e-14_dp) b_norm = 1.0_dp
      
      ! 第一次预处理
      z = r * self%inv_diag
      p = z
      rho = dot_product(r, z)

      ! ====================================================
      ! 📊 实验一：打开文件用于记录残差收敛历史
      ! ====================================================
      open(unit=105, file='pcg_residual.txt', status='replace')
      write(105, *) "Iteration    Relative_Residual"
      
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
          
          ! 计算当前的相对误差
          r_norm = sqrt(sum(r**2))
          current_err = r_norm / b_norm
          
          ! ================================================
          ! 📊 实验一：实时将当前迭代次数和残差写入文件
          ! ================================================
          write(105, *) iter, current_err
          
          if (current_err < tol) then
              print *, "  => [自研 PCG] 成功收敛! 迭代次数: ", iter
              exit
          end if
          
          ! 应用预处理
          z = r * self%inv_diag
          rho_new = dot_product(r, z)
          beta = rho_new / rho
          p = z + beta * p
          rho = rho_new
      end do
      
      ! ====================================================
      ! 📊 实验一：关闭记录文件
      ! ====================================================
      close(105)
      
      if (iter > max_iter) print *, "  => [PCG 警告] 达到最大迭代次数未收敛！"
      deallocate(r, z, p, q)
  end subroutine pcg_solve

  ! ========================================================
  ! 5. 资源释放阶段
  ! ========================================================
  subroutine pcg_free(self)
      class(PCGSolver), intent(inout) :: self
      if (allocated(self%inv_diag)) deallocate(self%inv_diag)
      nullify(self%A_csr)
      self%analyzed = .false.
      self%factorized = .false.
  end subroutine pcg_free

end module PCGSolver_mod