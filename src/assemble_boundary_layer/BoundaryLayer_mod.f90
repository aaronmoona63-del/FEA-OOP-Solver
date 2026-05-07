module BoundaryLayer_mod
  !=============================================================
  ! BoundaryLayer_mod
  !
  ! 目标：
  !   - 继续使用 legacy 的 module Boundaryconditions 作为“BC 数据库”
  !   - 提供一个可扩展的 BoundaryLayer：
  !       1) build_plan():  将 BC 数据“编译”为可执行的 BCPlan（扁平列表）
  !       2) apply_loads():  对 rhs / (K,rhs) 施加 Neumann/traction 等（可扩展）
  !       3) apply_dirichlet(): 对 (K,rhs) 施加强制位移（推荐在 CSR 上做）
  !
  ! 重要说明：
  !   - 本模块不绑定具体 MeshBase / DofMap 实现，避免耦合。
  !   - 映射 node+dof -> eq_id 通过回调 get_eq() 注入（你工程里怎么编号都行）。
  !   - Dirichlet 的“矩阵修改策略”通过回调 apply_dirichlet_op() 注入：
  !       例如：你可以在 CSR 上实现行列消元，然后把过程指针传进来。
  !
  ! 这样设计的好处：
  !   - BoundaryLayer 与矩阵格式/装配层解耦（以后 COO/CSR/Skyline 都能用）
  !   - Boundaryconditions 数据结构不动
  !   - 后续拓展：history / user subroutine / MPC / traction-by-marker 都有位置可挂
  !=============================================================
  use Types,  only : prec
  use ParamIO, only : IOW
  use Boundaryconditions
  implicit none
  private

  public :: BCPlan, BoundaryLayer
  public :: BC_SRC_CONST, BC_SRC_HISTORY, BC_SRC_USER
  public :: BC_OBJ_NODE, BC_OBJ_NODESET

  !------------------------
  ! Source types (values)
  !------------------------
  integer, parameter :: BC_SRC_CONST   = 1
  integer, parameter :: BC_SRC_HISTORY = 2
  integer, parameter :: BC_SRC_USER    = 3

  !------------------------
  ! Object types (targets)
  !------------------------
  integer, parameter :: BC_OBJ_NODE    = 1
  integer, parameter :: BC_OBJ_NODESET = 2

  !=============================================================
  ! Callback interfaces (依赖注入)
  !=============================================================
  abstract interface
    ! 将 (node_id, dof_id) 映射为全局方程号 eq（1-based）
    integer function get_eq_iface(node_id, dof_id) result(eq)
      import :: prec
      integer, intent(in) :: node_id
      integer, intent(in) :: dof_id
    end function get_eq_iface

    ! 对给定的 Dirichlet 列表施加到线性系统上：
    ! 你可以让它操作 CSR / Skyline / Dense 等具体结构。
    subroutine apply_dirichlet_op_iface(dbc_eq, dbc_val, rhs)
      import :: prec
      integer,   intent(in)    :: dbc_eq(:)
      real(prec),intent(in)    :: dbc_val(:)
      real(prec),intent(inout) :: rhs(:)
    end subroutine apply_dirichlet_op_iface
  end interface

  !=============================================================
  ! BCPlan: 扁平化“执行计划”
  !=============================================================
  type :: BCPlan
    !----------------------------
    ! Dirichlet (prescribed DOF)
    !----------------------------
    integer, allocatable :: dbc_eq(:)         ! 全局方程号
    integer, allocatable :: dbc_src(:)        ! BC_SRC_*
    real(prec), allocatable :: dbc_const(:)   ! 常值（若 src=CONST）
    integer, allocatable :: dbc_hist(:)       ! history id（若 src=HISTORY）
    integer, allocatable :: dbc_user(:)       ! 预留（若 src=USER）
    integer, allocatable :: dbc_rate(:)       ! rate_flag（预留）

    !----------------------------
    ! Nodal forces (prescribed forces)
    !----------------------------
    integer, allocatable :: fbc_eq(:)
    integer, allocatable :: fbc_src(:)
    real(prec), allocatable :: fbc_const(:)
    integer, allocatable :: fbc_hist(:)
    integer, allocatable :: fbc_user(:)

    !----------------------------
    ! Traction by marker (预留：推荐用 marker 绑定)
    !----------------------------
    integer, allocatable :: tr_marker(:)      ! unique markers
    integer, allocatable :: tr_src(:)
    integer, allocatable :: tr_hist(:)
    integer, allocatable :: tr_val_index(:)   ! index into dload_values
    integer, allocatable :: tr_ncomp(:)

    logical :: built = .false.
  contains
    procedure :: clear => plan_clear
  end type BCPlan

  !=============================================================
  ! BoundaryLayer: 对外使用的边界层对象
  !=============================================================
  type :: BoundaryLayer
    type(BCPlan) :: plan
    procedure(get_eq_iface), pointer, nopass :: get_eq => null()
    procedure(apply_dirichlet_op_iface), pointer, nopass :: apply_dirichlet_op => null()
  contains
    procedure :: attach_get_eq
    procedure :: attach_apply_dirichlet_op
    procedure :: build_plan
    procedure :: apply_loads
    procedure :: apply_dirichlet
  end type BoundaryLayer

contains

  !=============================================================
  ! BCPlan utilities
  !=============================================================
  subroutine plan_clear(self)
    class(BCPlan), intent(inout) :: self
    if (allocated(self%dbc_eq))       deallocate(self%dbc_eq)
    if (allocated(self%dbc_src))      deallocate(self%dbc_src)
    if (allocated(self%dbc_const))    deallocate(self%dbc_const)
    if (allocated(self%dbc_hist))     deallocate(self%dbc_hist)
    if (allocated(self%dbc_user))     deallocate(self%dbc_user)
    if (allocated(self%dbc_rate))     deallocate(self%dbc_rate)

    if (allocated(self%fbc_eq))       deallocate(self%fbc_eq)
    if (allocated(self%fbc_src))      deallocate(self%fbc_src)
    if (allocated(self%fbc_const))    deallocate(self%fbc_const)
    if (allocated(self%fbc_hist))     deallocate(self%fbc_hist)
    if (allocated(self%fbc_user))     deallocate(self%fbc_user)

    if (allocated(self%tr_marker))    deallocate(self%tr_marker)
    if (allocated(self%tr_src))       deallocate(self%tr_src)
    if (allocated(self%tr_hist))      deallocate(self%tr_hist)
    if (allocated(self%tr_val_index)) deallocate(self%tr_val_index)
    if (allocated(self%tr_ncomp))     deallocate(self%tr_ncomp)

    self%built = .false.
  end subroutine plan_clear

  !=============================================================
  ! Attach callbacks
  !=============================================================
  subroutine attach_get_eq(self, f)
    class(BoundaryLayer), intent(inout) :: self
    procedure(get_eq_iface), pointer, intent(in) :: f
    self%get_eq => f
  end subroutine attach_get_eq

  subroutine attach_apply_dirichlet_op(self, f)
    class(BoundaryLayer), intent(inout) :: self
    procedure(apply_dirichlet_op_iface), pointer, intent(in) :: f
    self%apply_dirichlet_op => f
  end subroutine attach_apply_dirichlet_op

  !=============================================================
  ! Build plan: compile Boundaryconditions DB into BCPlan
  !
  ! 目前实现：
  !   - prescribeddof_list: node_number 单点（node_set 将来可扩展）
  !   - prescribedforce_list: node_number 单点
  !
  ! 你之后扩展：
  !   - node_set: 展开 nodeset_list + node_lists
  !   - history: interpolate_history_table
  !   - user: 回调/插件
  !   - traction: 按 marker 或 elementset+face 展开
  !=============================================================
  subroutine build_plan(self)
    class(BoundaryLayer), intent(inout) :: self

    integer :: i, cnt
    integer :: node_id, dof_id, eq
    integer :: src
    real(prec) :: vconst

    if (.not. associated(self%get_eq)) then
      write(IOW,*) "BoundaryLayer: get_eq callback not attached."
      stop
    end if

    call self%plan%clear()

    !----------------------------
    ! 1) Dirichlet: prescribeddof_list
    !----------------------------
    cnt = 0
    do i = 1, n_prescribeddof
      if (prescribeddof_list(i)%node_number > 0) cnt = cnt + 1
      ! node_set 扩展：这里以后把 nodeset 展开成多个 node
    end do

    if (cnt > 0) then
      allocate(self%plan%dbc_eq(cnt))
      allocate(self%plan%dbc_src(cnt))
      allocate(self%plan%dbc_const(cnt))
      allocate(self%plan%dbc_hist(cnt))
      allocate(self%plan%dbc_user(cnt))
      allocate(self%plan%dbc_rate(cnt))

      cnt = 0
      do i = 1, n_prescribeddof
        node_id = prescribeddof_list(i)%node_number
        if (node_id <= 0) cycle

        dof_id = prescribeddof_list(i)%dof
        eq = self%get_eq(node_id, dof_id)

        src = prescribeddof_list(i)%flag
        select case (src)
        case (1) ! direct value
          vconst = dof_values(prescribeddof_list(i)%index_dof_values)
          self%plan%dbc_src(cnt+1)   = BC_SRC_CONST
          self%plan%dbc_const(cnt+1) = vconst
          self%plan%dbc_hist(cnt+1)  = 0
          self%plan%dbc_user(cnt+1)  = 0
        case (2) ! history
          self%plan%dbc_src(cnt+1)   = BC_SRC_HISTORY
          self%plan%dbc_const(cnt+1) = 0.0_prec
          self%plan%dbc_hist(cnt+1)  = prescribeddof_list(i)%history_number
          self%plan%dbc_user(cnt+1)  = 0
        case (3) ! user subroutine
          self%plan%dbc_src(cnt+1)   = BC_SRC_USER
          self%plan%dbc_const(cnt+1) = 0.0_prec
          self%plan%dbc_hist(cnt+1)  = 0
          self%plan%dbc_user(cnt+1)  = prescribeddof_list(i)%subroutine_parameter_number
        case default
          write(IOW,*) "BoundaryLayer: unsupported prescribeddof flag=", src
          stop
        end select

        self%plan%dbc_eq(cnt+1)   = eq
        self%plan%dbc_rate(cnt+1) = prescribeddof_list(i)%rate_flag
        cnt = cnt + 1
      end do
    end if

    !----------------------------
    ! 2) Nodal force: prescribedforce_list
    !----------------------------
    cnt = 0
    do i = 1, n_prescribedforces
      if (prescribedforce_list(i)%node_number > 0) cnt = cnt + 1
    end do

    if (cnt > 0) then
      allocate(self%plan%fbc_eq(cnt))
      allocate(self%plan%fbc_src(cnt))
      allocate(self%plan%fbc_const(cnt))
      allocate(self%plan%fbc_hist(cnt))
      allocate(self%plan%fbc_user(cnt))

      cnt = 0
      do i = 1, n_prescribedforces
        node_id = prescribedforce_list(i)%node_number
        if (node_id <= 0) cycle

        dof_id = prescribedforce_list(i)%dof
        eq = self%get_eq(node_id, dof_id)

        src = prescribedforce_list(i)%flag
        select case (src)
        case (1) ! direct value
          vconst = dof_values(prescribedforce_list(i)%index_dof_values)
          self%plan%fbc_src(cnt+1)   = BC_SRC_CONST
          self%plan%fbc_const(cnt+1) = vconst
          self%plan%fbc_hist(cnt+1)  = 0
          self%plan%fbc_user(cnt+1)  = 0
        case (2) ! history
          self%plan%fbc_src(cnt+1)   = BC_SRC_HISTORY
          self%plan%fbc_const(cnt+1) = 0.0_prec
          self%plan%fbc_hist(cnt+1)  = prescribedforce_list(i)%history_number
          self%plan%fbc_user(cnt+1)  = 0
        case (3) ! user subroutine
          self%plan%fbc_src(cnt+1)   = BC_SRC_USER
          self%plan%fbc_const(cnt+1) = 0.0_prec
          self%plan%fbc_hist(cnt+1)  = 0
          self%plan%fbc_user(cnt+1)  = prescribedforce_list(i)%subroutine_parameter_number
        case default
          write(IOW,*) "BoundaryLayer: unsupported prescribedforce flag=", src
          stop
        end select

        self%plan%fbc_eq(cnt+1) = eq
        cnt = cnt + 1
      end do
    end if

    !----------------------------
    ! 3) Traction (预留：你后续可以在这里编译)
    !    推荐：先按 marker 绑定；或者 elementset+face 展开到 bfaces
    !----------------------------
    ! TODO

    self%plan%built = .true.
  end subroutine build_plan

  !=============================================================
  ! Apply loads (Neumann-like): add to rhs and/or add boundary contributions
  !
  ! 当前实现：
  !   - nodal force：rhs(eq) += value(time)
  !
  ! 后续扩展：
  !   - traction：遍历 mesh bfaces -> marker -> boundary kernel -> add to COO/rhs
  !=============================================================
  subroutine apply_loads(self, time, rhs)
    class(BoundaryLayer), intent(inout) :: self
    real(prec), intent(in) :: time
    real(prec), intent(inout) :: rhs(:)

    integer :: i, eq
    real(prec) :: val

    if (.not. self%plan%built) then
      write(IOW,*) "BoundaryLayer: plan not built. Call build_plan() first."
      stop
    end if

    if (allocated(self%plan%fbc_eq)) then
      do i = 1, size(self%plan%fbc_eq)
        eq = self%plan%fbc_eq(i)
        val = eval_value(self%plan%fbc_src(i), self%plan%fbc_const(i), self%plan%fbc_hist(i), time)
        rhs(eq) = rhs(eq) + val
      end do
    end if

  end subroutine apply_loads

  !=============================================================
  ! Apply Dirichlet: compute dbc values at time and delegate to apply_dirichlet_op
  !
  ! 注意：
  !   - 这里不直接改 K（因为 K 的格式可能是 CSR/Skyline/Dense）
  !   - 你把具体实现（例如 CSR 行列消元）作为回调 attach 进来
  !=============================================================
  subroutine apply_dirichlet(self, time, rhs)
    class(BoundaryLayer), intent(inout) :: self
    real(prec), intent(in) :: time
    real(prec), intent(inout) :: rhs(:)

    real(prec), allocatable :: dbc_val(:)
    integer :: i

    if (.not. self%plan%built) then
      write(IOW,*) "BoundaryLayer: plan not built. Call build_plan() first."
      stop
    end if

    if (.not. associated(self%apply_dirichlet_op)) then
      write(IOW,*) "BoundaryLayer: apply_dirichlet_op callback not attached."
      stop
    end if

    if (.not. allocated(self%plan%dbc_eq)) return
    allocate(dbc_val(size(self%plan%dbc_eq)))

    do i = 1, size(self%plan%dbc_eq)
      dbc_val(i) = eval_value(self%plan%dbc_src(i), self%plan%dbc_const(i), self%plan%dbc_hist(i), time)
      ! rate_flag（预留）：如果你做瞬态并用 rate_flag=1，这里可改成 eval_rate(...)
    end do

    call self%apply_dirichlet_op(self%plan%dbc_eq, dbc_val, rhs)
    deallocate(dbc_val)

  end subroutine apply_dirichlet

  !=============================================================
  ! eval_value: evaluate const/history/user source at given time
  !=============================================================
  real(prec) function eval_value(src, cval, hist_id, time) result(v)
    integer, intent(in) :: src
    real(prec), intent(in) :: cval
    integer, intent(in) :: hist_id
    real(prec), intent(in) :: time

    real(prec) :: tmp
    integer :: idx, nh

    select case (src)
    case (BC_SRC_CONST)
      v = cval

    case (BC_SRC_HISTORY)
      if (hist_id <= 0 .or. hist_id > n_histories) then
        write(IOW,*) "BoundaryLayer: invalid history id=", hist_id
        stop
      end if
      idx = history_list(hist_id)%index
      nh  = history_list(hist_id)%n_timevals
      ! history_data is (2, length_history_data) in legacy, your snippet showed (:,:) anyway.
      ! We assume history_data(:, idx:idx+nh-1) stores [time; value] pairs.
      call interpolate_history_table(history_data(:, idx:idx+nh-1), nh, time, tmp)
      v = tmp

    case (BC_SRC_USER)
      ! 预留：user subroutine
      ! 你可以在这里调用 user_bc_value(hist_id, time, ...) 或通过 procedure pointer 注入
      write(IOW,*) "BoundaryLayer: USER source not implemented yet."
      stop

    case default
      write(IOW,*) "BoundaryLayer: unknown src=", src
      stop
    end select
  end function eval_value

end module BoundaryLayer_mod
