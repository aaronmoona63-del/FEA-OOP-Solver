module LinearSolverFactory_mod
    use LinearSolverBase_mod
    use SkylineDirectSolver_mod
    use MKLDirectSolver_mod
    use PCGSolver_mod
    use CGSolver_mod     ! <--- 新增引入标准 CG

    use SparseSkyline_mod, only: SparseSkyline
    use SparseCSR_mod,     only: SparseCSR
    use MatrixDense_mod,   only: DenseMatrix

    implicit none
    private
    public :: create_linear_solver

contains

    subroutine create_linear_solver(A, opts, solver)
        class(*), intent(in), target                  :: A
        type(LinearSolverOptions), intent(in)         :: opts
        class(LinearSolver), allocatable, intent(out) :: solver

        select type(pA => A)
        type is (SparseSkyline)
            allocate(SkylineDirectSolver :: solver)
            solver%opts = opts
            call solver%attach_matrix(pA)

        type is (SparseCSR)
            ! 👇 根据传入的名字，智能分发求解器
            if (opts%solver_family == "pcg" .or. opts%solver_family == "iter") then
                allocate(PCGSolver :: solver)
            elseif (opts%solver_family == "cg") then
                allocate(CGSolver :: solver)  ! <--- 调用新写的标准 CG
            else
                allocate(MKLDirectSolver :: solver)
            end if
            solver%opts = opts
            call solver%attach_matrix(pA)

        type is (DenseMatrix)
            allocate(MKLDirectSolver :: solver)
            solver%opts = opts
            call solver%attach_matrix(pA)
        class default
            stop "create_linear_solver: unsupported matrix type."
        end select
    end subroutine create_linear_solver
end module LinearSolverFactory_mod