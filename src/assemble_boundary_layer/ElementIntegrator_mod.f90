!> ElementIntegrator_mod.f90
!> ------------------------------------------------------------
!> 角色定位：
!>   - 负责“算局部贡献”：遍历单元/边界面，调用 ElementKernel 生成 (Ke, Fe)
!>   - 不负责全局稀疏矩阵的 scatter-add（那是 Assembler 的职责）
!>
!> 依赖关系（有意为之）：
!>   ElementIntegrator 依赖 MeshBase / FESpace / DofMap / ElementKernel
!>   Assembler 只依赖 Sparse 矩阵类型（已在你项目里拆分）
!>
!> 说明（结合你的 legacy Boundaryconditions）：
!>   - 建议把 distributedloads / traction / marker 预处理成 mesh 的 bfaces + marker
!>     然后在 kernel%compute_bface(...) 里根据 marker 和输入参数决定 Kb/fbv。
!>   - 本模块“只遍历 mesh%get_num_bfaces()”，不直接读 legacy Boundaryconditions。
!>     （这样 Integrator 干净、通用、可并行；legacy BC 的复杂语义放到“预处理/编译器”层）
!> ------------------------------------------------------------

module ElementIntegrator_mod
  use LinearSolverBase_mod,  only: dp
  use MeshBase_mod,          only: MeshBase
  use FESpace_mod,           only: FESpace
  use DofMap_mod,            only: DofMap
  use ElementKernelBase_mod, only: ElementKernel
  use Assembler_mod,         only: ElementContribution
  implicit none
  private
  public :: ElementIntegrator

  type :: ElementIntegrator
    class(ElementKernel), pointer :: kernel => null()
    type(DofMap),         pointer :: dofmap => null()

    ! 控制项
    logical :: do_domain   = .true.
    logical :: do_boundary = .true.
  contains
    procedure :: attach_kernel
    procedure :: attach_dofmap

    ! 生成贡献（局部）：
    procedure :: integrate_domain
    procedure :: integrate_boundary

    ! 便捷封装：一次性拿到域+边界贡献
    procedure :: integrate_all
  end type ElementIntegrator

contains

  subroutine attach_kernel(self, ker)
    class(ElementIntegrator), intent(inout) :: self
    class(ElementKernel), target, intent(inout) :: ker
    self%kernel => ker
  end subroutine attach_kernel

  subroutine attach_dofmap(self, dm)
    class(ElementIntegrator), intent(inout) :: self
    type(DofMap), target, intent(inout) :: dm
    self%dofmap => dm
  end subroutine attach_dofmap

  !------------------------------------------------------------
  ! integrate_domain: ∫Ω 逐单元调用 kernel%compute_element
  ! 输出：econtribs(:) 每个元素一个 ElementContribution
  !------------------------------------------------------------
  subroutine integrate_domain(self, mesh, fe, econtribs)
    class(ElementIntegrator), intent(inout) :: self
    class(MeshBase), intent(in) :: mesh
    type(FESpace),   intent(in) :: fe
    type(ElementContribution), allocatable, intent(out) :: econtribs(:)

    integer :: ne, e
    integer, allocatable :: gdofs(:)
    real(dp), allocatable :: k_el(:,:), f_el(:)

    if (.not. associated(self%kernel)) stop "ElementIntegrator: kernel not attached."
    if (.not. associated(self%dofmap)) stop "ElementIntegrator: dofmap not attached."

    if (.not. self%do_domain) then
      allocate(econtribs(0))
      return
    end if

    ne = mesh%get_num_elems()
    allocate(econtribs(ne))

    do e = 1, ne
      call self%dofmap%get_elem_dofs(mesh, e, gdofs)
      call self%kernel%compute_element(mesh, fe, e, k_el, f_el)

      call set_contrib(econtribs(e), gdofs, k_el, f_el)

      call safe_dealloc_int(gdofs)
      call safe_dealloc_mat(k_el)
      call safe_dealloc_vec(f_el)
    end do
  end subroutine integrate_domain

  !------------------------------------------------------------
  ! integrate_boundary: ∫Γ 逐边界面调用 kernel%compute_bface
  ! 输出：bcontribs(:) 每个边界面一个 ElementContribution（或空贡献）
  !
  ! 约定：
  !   - marker = mesh%get_bface_marker(bf)
  !   - kernel%compute_bface(mesh, fe, bf, marker, Kb, fbv)
  !     若该 marker 不需要贡献，可令 Kb/fbv 不分配或分配为 size=0。
  !------------------------------------------------------------
  subroutine integrate_boundary(self, mesh, fe, bcontribs)
    class(ElementIntegrator), intent(inout) :: self
    class(MeshBase), intent(in) :: mesh
    type(FESpace),   intent(in) :: fe
    type(ElementContribution), allocatable, intent(out) :: bcontribs(:)

    integer :: nbf, bf, marker
    integer :: n_keep
    integer, allocatable :: gdofs(:)
    real(dp), allocatable :: Kb(:,:), fbv(:)
    type(ElementContribution), allocatable :: tmp(:)

    if (.not. associated(self%kernel)) stop "ElementIntegrator: kernel not attached."
    if (.not. associated(self%dofmap)) stop "ElementIntegrator: dofmap not attached."

    if (.not. self%do_boundary) then
      allocate(bcontribs(0))
      return
    end if

    nbf = mesh%get_num_bfaces()
    if (nbf <= 0) then
      allocate(bcontribs(0))
      return
    end if

    ! 先按最大可能分配，再压缩
    allocate(tmp(nbf))
    n_keep = 0

    do bf = 1, nbf
      marker = mesh%get_bface_marker(bf)

      call self%dofmap%get_bface_dofs(mesh, bf, gdofs)
      call self%kernel%compute_bface(mesh, fe, bf, marker, Kb, fbv)

      if (has_nonempty_contrib(Kb, fbv, gdofs)) then
        n_keep = n_keep + 1
        call set_contrib(tmp(n_keep), gdofs, Kb, fbv)
      end if

      call safe_dealloc_int(gdofs)
      call safe_dealloc_mat(Kb)
      call safe_dealloc_vec(fbv)
    end do

    if (n_keep == 0) then
      allocate(bcontribs(0))
      deallocate(tmp)
    else
      allocate(bcontribs(n_keep))
      bcontribs = tmp(1:n_keep)
      deallocate(tmp)
    end if
  end subroutine integrate_boundary

  !------------------------------------------------------------
  ! integrate_all: 返回域贡献 + 边界贡献
  !------------------------------------------------------------
  subroutine integrate_all(self, mesh, fe, econtribs, bcontribs)
    class(ElementIntegrator), intent(inout) :: self
    class(MeshBase), intent(in) :: mesh
    type(FESpace),   intent(in) :: fe
    type(ElementContribution), allocatable, intent(out) :: econtribs(:)
    type(ElementContribution), allocatable, intent(out) :: bcontribs(:)

    call self%integrate_domain(mesh, fe, econtribs)
    call self%integrate_boundary(mesh, fe, bcontribs)
  end subroutine integrate_all

  !============================================================
  ! Helpers
  !============================================================

  subroutine set_contrib(ec, gdofs, ke, fev)
    type(ElementContribution), intent(inout) :: ec
    integer, intent(in) :: gdofs(:)
    real(dp), intent(in) :: ke(:,:), fev(:)

    ! 深拷贝（allocatable assignment 会自动分配并复制）
    ec%gdofs = gdofs
    ec%ke    = ke
    ec%fe    = fev
  end subroutine set_contrib

  logical function has_nonempty_contrib(Kb, fbv, gdofs)
    real(dp), allocatable, intent(in) :: Kb(:,:), fbv(:)
    integer, allocatable, intent(in) :: gdofs(:)

    has_nonempty_contrib = .false.

    if (.not. allocated(gdofs)) return
    if (size(gdofs) == 0) return

    if (allocated(Kb)) then
      if (size(Kb,1) > 0 .and. size(Kb,2) > 0) then
        has_nonempty_contrib = .true.
        return
      end if
    end if

    if (allocated(fbv)) then
      if (size(fbv) > 0) then
        has_nonempty_contrib = .true.
        return
      end if
    end if
  end function has_nonempty_contrib

  subroutine safe_dealloc_int(a)
    integer, allocatable, intent(inout) :: a(:)
    if (allocated(a)) deallocate(a)
  end subroutine safe_dealloc_int

  subroutine safe_dealloc_vec(a)
    real(dp), allocatable, intent(inout) :: a(:)
    if (allocated(a)) deallocate(a)
  end subroutine safe_dealloc_vec

  subroutine safe_dealloc_mat(a)
    real(dp), allocatable, intent(inout) :: a(:,:)
    if (allocated(a)) deallocate(a)
  end subroutine safe_dealloc_mat

end module ElementIntegrator_mod
