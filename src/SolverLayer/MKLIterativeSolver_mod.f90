module LinearSolverFactory_mod
    use LinearSolverBase_mod
    use SkylineDirectSolver_mod
    use MKLDirectSolver_mod
    ! 👇 换成我们自主研发的 PCG
    use PCGSolver_mod

    use SparseSkyline_mod, only: SparseSkyline
    use SparseCSR_mod,      only: SparseCSR
    use MatrixDense_mod,    only: DenseMatrix

    implicit none
    private
    public :: create_linear_solver

contains

    subroutine create_linear_solver(A, opts, solver)
        class(*), intent(in), target              :: A
        type(LinearSolverOptions), intent(in)     :: opts
        class(LinearSolver), allocatable, intent(out) :: solver

        select type(pA => A)
        type is (SparseSkyline)
            allocate(SkylineDirectSolver :: solver)
            solver%opts = opts
            call solver%attach_matrix(pA)

        type is (SparseCSR)
            if (opts%solver_family == "iter") then
                ! 👇 接入我们的自主研发 PCG 求解器
                allocate(PCGSolver :: solver)
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