module BCPlan_mod
  !=============================================================
  ! BCPlan_mod
  !
  ! 定义 Boundary Condition 的“执行计划（Compiled Plan）”
  !
  ! 设计目标：
  !   - 承载从 Boundaryconditions(DB) 编译后的结果
  !   - 数据应当是：
  !       * 扁平的
  !       * 可缓存的
  !       * 执行期零判断/零解析
  !
  ! 本模块只定义数据结构和最基本的管理操作
  ! 不做任何：
  !   - Mesh 查询
  !   - 矩阵修改
  !   - 时间积分
  !
  !=============================================================
  use Types, only : prec
  implicit none
  private

  public :: BCPlan
  public :: BC_SRC_CONST, BC_SRC_HISTORY, BC_SRC_USER

  !-------------------------------------------------------------
  ! Source type flags（值的来源）
  !-------------------------------------------------------------
  integer, parameter :: BC_SRC_CONST   = 1
  integer, parameter :: BC_SRC_HISTORY = 2
  integer, parameter :: BC_SRC_USER    = 3

  !=============================================================
  ! BCPlan type
  !
  ! 含义：
  !   每个数组元素代表“一条已经展开好的 BC”
  !
  !   例如 Dirichlet：
  !     dbc_eq(i)    = 第 i 条 BC 作用的全局方程号
  !     dbc_src(i)   = 值来源（CONST/HISTORY/USER）
  !=============================================================
  type :: BCPlan

    !----------------------------
    ! Dirichlet BC（强制位移）
    !----------------------------
    integer, allocatable :: dbc_eq(:)         ! 全局方程号
    integer, allocatable :: dbc_src(:)        ! BC_SRC_*
    real(prec), allocatable :: dbc_const(:)   ! 常值（src=CONST）
    integer, allocatable :: dbc_hist(:)       ! history id（src=HISTORY）
    integer, allocatable :: dbc_user(:)       ! user parameter id（src=USER）
    integer, allocatable :: dbc_rate(:)       ! rate_flag（预留：瞬态）

    !----------------------------
    ! Nodal force（节点力）
    !----------------------------
    integer, allocatable :: fbc_eq(:)
    integer, allocatable :: fbc_src(:)
    real(prec), allocatable :: fbc_const(:)
    integer, allocatable :: fbc_hist(:)
    integer, allocatable :: fbc_user(:)

    !----------------------------
    ! Distributed load / traction
    ! （推荐 marker-based，便于和 MeshBase 对齐）
    !----------------------------
    integer, allocatable :: tr_marker(:)      ! boundary marker
    integer, allocatable :: tr_src(:)
    integer, allocatable :: tr_hist(:)
    integer, allocatable :: tr_val_index(:)   ! index into dload_values
    integer, allocatable :: tr_ncomp(:)       ! 分量数（2D=2, 3D=3）

    logical :: built = .false.

  contains
    procedure :: clear => bcplan_clear
  end type BCPlan

contains

  !=============================================================
  ! 清空 BCPlan
  !
  ! 说明：
  !   - 用于重新 build_plan
  !   - 不保留任何旧的 BC
  !=============================================================
  subroutine bcplan_clear(self)
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
  end subroutine bcplan_clear

end module BCPlan_mod
