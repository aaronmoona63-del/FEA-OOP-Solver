module LinearSolverBase_mod
    implicit none
    private
    public :: LinearSolver, LinearSolverOptions, dp

    integer, parameter :: dp = kind(1.0d0)

    ! ==========================
    ! 选项结构：控制求解策略
    ! ==========================
    type :: LinearSolverOptions
        character(len=16) :: solver_family = "auto"
        ! "auto" / "direct" / "iter"

        character(len=16) :: storage_pref  = "auto"
        ! "auto" / "skyline" / "csr" / "dense"

        real(dp) :: tol      = 1.0d-8     ! 迭代法容差
        integer  :: max_iter = 1000       ! 最大迭代次数
        integer  :: max_refinement = 0    ! 迭代改进次数（直接法可用）

        integer  :: verbosity = 0         ! 输出级别 0/1/2
    end type LinearSolverOptions

    ! ==========================
    ! 抽象基类：LinearSolver
    ! ==========================
    type, abstract :: LinearSolver
        type(LinearSolverOptions) :: opts
        logical :: analyzed  = .false.
        logical :: factorized = .false.
    contains
        procedure(ls_attach_iface), deferred :: attach_matrix
        procedure(ls_analyze_iface), deferred :: analyze
        procedure(ls_factor_iface),  deferred :: factor
        procedure(ls_solve_iface),   deferred :: solve
        procedure(ls_free_iface),    deferred :: free
    end type LinearSolver

    ! ==========================
    ! 抽象接口定义
    ! ==========================
    abstract interface
        ! 绑定矩阵（不同 solver 自己 downcast）
        subroutine ls_attach_iface(self, A)
            import :: LinearSolver
            class(LinearSolver), intent(inout) :: self
            class(*), intent(in), target       :: A
        end subroutine ls_attach_iface

        ! 可选分析阶段（PARDISO 之类用）
        subroutine ls_analyze_iface(self)
            import :: LinearSolver
            class(LinearSolver), intent(inout) :: self
        end subroutine ls_analyze_iface

        ! 分解（直接法）或预处理（迭代法）
        subroutine ls_factor_iface(self)
            import :: LinearSolver
            class(LinearSolver), intent(inout) :: self
        end subroutine ls_factor_iface

        ! 求解：给定右端 b，返回 x
        subroutine ls_solve_iface(self, b, x)
            import :: LinearSolver, dp
            class(LinearSolver), intent(inout) :: self
            real(dp), intent(inout)  :: b(:)
            real(dp), intent(out) :: x(:)
        end subroutine ls_solve_iface

        ! 释放内部资源（工作区、MKL handler 等）
        subroutine ls_free_iface(self)
            import :: LinearSolver
            class(LinearSolver), intent(inout) :: self
        end subroutine ls_free_iface
    end interface

end module LinearSolverBase_mod
