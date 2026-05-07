module SkylineDirectSolver_mod
    use LinearSolverBase_mod, only: LinearSolver, dp
    use SparseSkyline_mod,   only: SparseSkyline
    implicit none
    private
    public :: SkylineDirectSolver


    !===========================================================
    ! SkylineDirectSolver : 对称 Skyline 矩阵的 LDL^T 直接解法
    ! 说明：
    !   * 只支持对称矩阵（A%is_symmetric = .true.）
    !   * 在 A 上原地做 LDL^T 分解：
    !       A = L * D * L^T
    !     - A%ad(j) 存 D(j)
    !     - A%al 存 L(j,i) (j>i) 的下三角，按列 Skyline 存储
    !     - A%au 在对称情况下简单设置为 A%al（保持格式）
    !   * solve(b,x)：
    !       1) 前代   L y = b
    !       2) 对角   D z = y
    !       3) 回代   L^T x = z
    !===========================================================
    type, extends(LinearSolver) :: SkylineDirectSolver
        type(SparseSkyline), pointer :: A => null()
        logical :: is_factorized = .false.
    contains
        procedure :: attach_matrix => sky_attach
        procedure :: analyze       => sky_analyze
        procedure :: factor        => sky_factor
        procedure :: solve         => sky_solve
        procedure :: free          => sky_free
    end type SkylineDirectSolver

contains

    !-----------------------------------------------------------
    !  绑定矩阵：必须是 SparseSkyline
    !-----------------------------------------------------------
    subroutine sky_attach(self, A)
        class(SkylineDirectSolver), intent(inout) :: self
        class(*),                  intent(in), target :: A

        select type(pA => A)
        type is (SparseSkyline)
            self%A => pA
        class default
            stop "SkylineDirectSolver: attach_matrix expects SparseSkyline."
        end select

        self%is_factorized = .false.
        self%analyzed      = .false.
        self%factorized    = .false.
    end subroutine sky_attach

    !-----------------------------------------------------------
    !  对 Skyline 来说，analyze 可以只是做一些检查
    !-----------------------------------------------------------
    subroutine sky_analyze(self)
        class(SkylineDirectSolver), intent(inout) :: self

        if (.not. associated(self%A)) then
            stop "SkylineDirectSolver: analyze() called before attach_matrix()."
        end if

        if (.not. self%A%is_symmetric) then
            stop "SkylineDirectSolver: only symmetric Skyline matrices are supported in this LDL^T version."
        end if

        self%analyzed = .true.
    end subroutine sky_analyze

    !-----------------------------------------------------------
    !  Skyline LDL^T 分解（原地修改 self%A）
    !  格式说明（SparseSkyline_mod 中）：
    !
    !   对于列 j：
    !     if j == 1:
    !         base = 0
    !     else
    !         base = jp(j-1)
    !     jh   = jp(j) - base          ! 这一列的“高度”（非零数）
    !     is   = j - jh                ! 这一列最顶端行号
    !
    !   对应元素：
    !     AL(base + k)   对应 (j, is + k - 1)
    !     AU(base + k)   对应 (is + k - 1, j)
    !
    !  LDL^T 算法大意：
    !    for j = 1..n:
    !      for 每个 i < j 且在 skyline 中:
    !         L(j,i) = ( A(j,i) - sum_{k in intersection} L(j,k)*D(k)*L(i,k) ) / D(i)
    !      D(j) = A(j,j) - sum_{i<j} L(j,i)^2 * D(i)
    !-----------------------------------------------------------
    subroutine sky_factor(self)
        class(SkylineDirectSolver), intent(inout) :: self

        type(SparseSkyline), pointer :: Ask
        integer :: n
        integer :: j, i, k
        integer :: base_j, h_j, is_j
        integer :: base_i, h_i, is_i
        integer :: loc
        integer :: k_start, k_end
        real(dp) :: sum, sumd, lij

        if (.not. associated(self%A)) then
            stop "SkylineDirectSolver: factor() called before attach_matrix()."
        end if

        Ask => self%A
        n   = Ask%n

        if (.not. Ask%is_symmetric) then
            stop "SkylineDirectSolver: factor() currently supports only symmetric Skyline matrices (LDL^T)."
        end if

        if (.not. allocated(Ask%jp) .or. .not. allocated(Ask%ad) .or. .not. allocated(Ask%al)) then
            stop "SkylineDirectSolver: SparseSkyline storage not allocated properly."
        end if

        !---------------- LDL^T 主循环 ----------------
        do j = 1, n

            ! 本列 j 的 skyline 结构
            if (j == 1) then
                base_j = 0
            else
                base_j = Ask%jp(j-1)
            end if
            h_j = Ask%jp(j) - base_j          ! 高度
            is_j = j - h_j                    ! 这一列最顶部行号

            !------------------------------------------------------
            ! 1) 计算列 j 的下三角 L(j,i), i = is_j..j-1
            !------------------------------------------------------
            do loc = 1, h_j
                i = is_j + loc - 1           ! 行号 i < j
                ! 当前位置在 AL 里的索引
                sum = Ask%al(base_j + loc)

                ! 行 i 的 skyline 信息（来自列 i）
                if (i == 1) then
                    base_i = 0
                else
                    base_i = Ask%jp(i-1)
                end if
                h_i  = Ask%jp(i) - base_i
                is_i = i - h_i

                ! 与行 i 的交集：k in [max(is_j,is_i) .. i-1]
                k_start = max(is_j, is_i)
                k_end   = i - 1

                if (k_start <= k_end) then
                    do k = k_start, k_end
                        ! L(j,k) 在列 j 的位置：
                        !   只有当 k >= is_j 时才在 skyline 中：
                        !   index_jk = base_j + (k - is_j + 1)
                        ! L(i,k) 在列 i 的位置：
                        !   index_ik = base_i + (k - is_i + 1)
                        ! 这两个位置一定在存储内，因为 k 已经在交集区间里

                        sum = sum - Ask%al( base_j + (k - is_j + 1) ) * &
                                    Ask%ad( k )                       * &
                                    Ask%al( base_i + (k - is_i + 1) )
                    end do
                end if

                ! L(j,i) = sum / D(i)
                Ask%al(base_j + loc) = sum / Ask%ad(i)
            end do

            !------------------------------------------------------
            ! 2) 更新对角 D(j)
            !    D(j) = A(j,j) - sum_{i<j} L(j,i)^2 * D(i)
            !------------------------------------------------------
            sumd = 0.0_dp
            do loc = 1, h_j
                i   = is_j + loc - 1
                lij = Ask%al(base_j + loc)
                sumd = sumd + lij * lij * Ask%ad(i)
            end do

            Ask%ad(j) = Ask%ad(j) - sumd

            ! 简单的奇异性检测（可根据需要改成 warn 而不是 stop）
            if (Ask%ad(j) == 0.0_dp) then
                stop "SkylineDirectSolver: zero pivot encountered in LDL^T factorization."
            end if
        end do

        ! 对称矩阵：通常不需要显式存 AU，但为了保持结构，我们把 AU 拷贝为 AL
        if (allocated(Ask%au)) then
            Ask%au = Ask%al
        end if

        self%is_factorized = .true.
        self%factorized    = .true.
    end subroutine sky_factor

    !-----------------------------------------------------------
    !  求解：L D L^T x = b
    !   1) forward:  L y = b
    !   2) diag:     D z = y
    !   3) backward: L^T x = z
    !-----------------------------------------------------------
    subroutine sky_solve(self, b, x)
        class(SkylineDirectSolver), intent(inout) :: self
        real(dp), intent(inout)  :: b(:)
        real(dp), intent(out) :: x(:)

        type(SparseSkyline), pointer :: Ask
        integer :: n
        integer :: j
        integer :: base_j, h_j, is_j
        real(dp), allocatable :: y(:)

        if (.not. associated(self%A)) then
            stop "SkylineDirectSolver: solve() called before attach_matrix()."
        end if
        if (.not. self%is_factorized) then
            stop "SkylineDirectSolver: solve() called before factor()."
        end if

        Ask => self%A
        n   = Ask%n

        if (size(b) /= n .or. size(x) /= n) then
            stop "SkylineDirectSolver: solve() size mismatch between b/x and matrix."
        end if

        allocate(y(n))
        y = b

        !---------------- 1) Forward: L y = b ----------------
        do j = 1, n
            if (j == 1) then
                base_j = 0
            else
                base_j = Ask%jp(j-1)
            end if
            h_j  = Ask%jp(j) - base_j
            is_j = j - h_j

            if (h_j > 0) then
                ! y(j) = y(j) - sum_{k=is_j}^{j-1} L(j,k) * y(k)
                y(j) = y(j) - dot_product( Ask%al(base_j+1:base_j+h_j), y(is_j:j-1) )
            end if
        end do

        !---------------- 2) Diagonal scaling: D z = y ----------------
        do j = 1, n
            y(j) = y(j) / Ask%ad(j)
        end do

        !---------------- 3) Backward: L^T x = z ----------------
        x = y
        do j = n, 1, -1
            if (j == 1) then
                base_j = 0
            else
                base_j = Ask%jp(j-1)
            end if
            h_j  = Ask%jp(j) - base_j
            is_j = j - h_j

            if (h_j > 0) then
                ! 对于每个 i = is_j..j-1:
                !   x(i) = x(i) - L(j,i) * x(j)
                x(is_j:j-1) = x(is_j:j-1) - x(j) * Ask%al(base_j+1:base_j+h_j)
            end if
        end do

        deallocate(y)
    end subroutine sky_solve

    !-----------------------------------------------------------
    !  释放：这里只是断开指针 & 状态置位
    !-----------------------------------------------------------
    !-----------------------------------------------------------
    ! 释放：这里只是断开指针 & 状态置位
    !-----------------------------------------------------------
    subroutine sky_free(self)
        class(SkylineDirectSolver), intent(inout) :: self

        nullify(self%A)
        self%is_factorized = .false.
        self%analyzed      = .false.
        self%factorized    = .false.
    end subroutine sky_free

end module SkylineDirectSolver_mod