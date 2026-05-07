module AssemblerCOO_mod
  use Types,        only : prec
  !use ParamIO,      only : IOW
  use Mesh,         only : mesh_query_element_sizes, mesh_get_element_for_assembly
  use SparseCOO_mod,only : SparseCOO
  use ElementDispatch_mod, only : element_dispatch

  implicit none
  private
  public :: assemble_all_elements_coo

contains

  subroutine assemble_all_elements_coo(n_elements, n_unknowns, Kcoo, rhs, fail)
    !==========================================================
    ! Assemble global stiffness and residual using:
    !   Mesh query -> allocate -> Mesh get -> Kernel -> COO add
    !
    ! Inputs:
    !   n_elements : number of elements (loop count)
    !   n_unknowns : global DOF count (matrix size)
    !
    ! In/Out:
    !   Kcoo       : SparseCOO global matrix
    !
    ! Outputs:
    !   rhs(:)     : global RHS/residual vector
    !   fail       : error flag
    !==========================================================
    integer, intent(in) :: n_elements
    integer, intent(in) :: n_unknowns
    type(SparseCOO), intent(out) :: Kcoo
    real(prec), allocatable, intent(out) :: rhs(:)
    logical, intent(out) :: fail
    
    ! local variables
    integer :: lmn
    integer :: flag, nnode, ndims, ndofpn, n_props, n_svars
    integer :: ndofs_elem

    ! element-local buffers
    real(prec), allocatable :: coords(:,:)        ! (ndims, nnode)
    real(prec), allocatable :: u_tot(:,:)         ! (ndofpn, nnode)
    real(prec), allocatable :: u_inc(:,:)         ! (ndofpn, nnode)
    integer,    allocatable :: gdofs(:)           ! (ndofs_elem=ndofpn*nnode)
    real(prec), allocatable :: props(:)
    real(prec), allocatable :: svars0(:)
    real(prec), allocatable :: svars(:)
    real(prec), allocatable :: Ke(:,:)            ! (ndofs_elem, ndofs_elem)
    real(prec), allocatable :: fe(:)              ! (ndofs_elem)

    ! init
    fail = .false.
    allocate(rhs(n_unknowns))
    rhs = 0.0_prec
    call Kcoo%init(n_unknowns)

    ! loop over elements
    do lmn = 1, n_elements

      !--------------------------------------------------------
      ! (1) Query sizes (NO direct access to element_list etc.)
      !--------------------------------------------------------
      call mesh_query_element_sizes(lmn, flag, nnode, ndims, ndofpn, n_props, n_svars)
      ndofs_elem = nnode * ndofpn

      !--------------------------------------------------------
      ! (2) Allocate (or reuse) element-local buffers
      !--------------------------------------------------------
      call prepare_element_workspace(nnode, ndims, ndofpn, ndofs_elem, n_props, n_svars, &
                               coords, u_tot, u_inc, gdofs, props, svars0, svars, Ke, fe)

      !--------------------------------------------------------
      ! (3) Get all assembler-ready data in ONE call
      !--------------------------------------------------------
      call mesh_get_element_for_assembly(lmn, flag, nnode, ndims, ndofpn, &
                                         coords, u_tot, u_inc, gdofs, &
                                         n_props, props, n_svars, svars0, svars)

      !--------------------------------------------------------
      ! (4) Call element_dispatch to compute -> Ke, fe
      !--------------------------------------------------------
      Ke = 0.0_prec
      fe = 0.0_prec

      call element_dispatch(flag, coords, u_tot, u_inc, &
                            props(1:n_props), svars0(1:n_svars), svars(1:n_svars), &
                            Ke, fe, fail)
      
      !--------------------------------------------------------
      ! (5) COO scatter-add
      !--------------------------------------------------------
      call add_vec(rhs, gdofs(1:ndofs_elem), fe(1:ndofs_elem))
      call add_block_to_coo(Kcoo, gdofs(1:ndofs_elem), Ke(1:ndofs_elem,1:ndofs_elem))

    end do

    if (fail) return

    !----------------------------------------------------------
    ! (6) Sort/merge duplicates
    !----------------------------------------------------------
    call Kcoo%preprocess()


    ! ---------------------------
    ! Check RHS
    ! ---------------------------
    print *, 'svars0 = ', svars0
    print *, 'svars = ', svars

    return

  end subroutine assemble_all_elements_coo


  !============================================================
  ! Helpers
  !============================================================

  subroutine add_vec(rhs, idx, v)
    real(prec), intent(inout) :: rhs(:)
    integer, intent(in) :: idx(:)
    real(prec), intent(in) :: v(:)
    integer :: i
    do i = 1, size(idx)
      rhs(idx(i)) = rhs(idx(i)) + v(i)
    end do
  end subroutine add_vec


  subroutine add_block_to_coo(A, gd, blk)
    type(SparseCOO), intent(inout) :: A
    integer, intent(in) :: gd(:)
    real(prec), intent(in) :: blk(:,:)
    integer :: i, j
    do i = 1, size(gd)
      do j = 1, size(gd)
        if (blk(i,j) /= 0.0_prec) call A%add_entry(gd(i), gd(j), blk(i,j))
      end do
    end do
  end subroutine add_block_to_coo


  subroutine prepare_element_workspace(nnode, ndims, ndofpn, ndofs_elem, n_props, n_svars, &
                                 coords, u_tot, u_inc, gdofs, props, svars0, svars, Ke, fe)
    ! allocate or reuse arrays
    integer, intent(in) :: nnode, ndims, ndofpn, ndofs_elem, n_props, n_svars
    real(prec), allocatable, intent(inout) :: coords(:,:), u_tot(:,:), u_inc(:,:)
    integer,    allocatable, intent(inout) :: gdofs(:)
    real(prec), allocatable, intent(inout) :: props(:), svars0(:), svars(:)
    real(prec), allocatable, intent(inout) :: Ke(:,:), fe(:)

    if (allocated(coords)) then
      if (size(coords,1) /= ndims .or. size(coords,2) /= nnode) deallocate(coords)
    end if
    if (.not. allocated(coords)) allocate(coords(ndims, nnode))

    if (allocated(u_tot)) then
      if (size(u_tot,1) /= ndofpn .or. size(u_tot,2) /= nnode) deallocate(u_tot)
    end if
    if (.not. allocated(u_tot)) allocate(u_tot(ndofpn, nnode))

    if (allocated(u_inc)) then
      if (size(u_inc,1) /= ndofpn .or. size(u_inc,2) /= nnode) deallocate(u_inc)
    end if
    if (.not. allocated(u_inc)) allocate(u_inc(ndofpn, nnode))

    if (allocated(gdofs)) then
      if (size(gdofs) /= ndofs_elem) deallocate(gdofs)
    end if
    if (.not. allocated(gdofs)) allocate(gdofs(ndofs_elem))

    if (allocated(Ke)) then
      if (size(Ke,1) /= ndofs_elem .or. size(Ke,2) /= ndofs_elem) deallocate(Ke)
    end if
    if (.not. allocated(Ke)) allocate(Ke(ndofs_elem, ndofs_elem))

    if (allocated(fe)) then
      if (size(fe) /= ndofs_elem) deallocate(fe)
    end if
    if (.not. allocated(fe)) allocate(fe(ndofs_elem))

    ! props/svars can be zero-length in theory; to keep code simple,
    ! allocate at least 1 and only use first n_props/n_svars entries.
    if (allocated(props)) then
      if (size(props) < max(1,n_props)) deallocate(props)
    end if
    if (.not. allocated(props)) allocate(props(max(1,n_props)))

    if (allocated(svars0)) then
      if (size(svars0) < max(1,n_svars)) deallocate(svars0)
    end if
    if (.not. allocated(svars0)) allocate(svars0(max(1,n_svars)))

    if (allocated(svars)) then
      if (size(svars) < max(1,n_svars)) deallocate(svars)
    end if
    if (.not. allocated(svars)) allocate(svars(max(1,n_svars)))

  end subroutine prepare_element_workspace

  subroutine assemble_constraints_coo(A, rhs)
    type(SparseCOO), intent(inout) :: A
    real(prec), intent(inout) :: rhs(:)

    ! 这里照抄你原来约束组装的数学意义，
    ! 只是把 skyline 写入改成 COO add_entry
    !
    ! 你原来做的事情包括：
    ! - tie constraint: 在 (icol,irow) 位置加 +/-1
    ! - rhs(icol) 写 constraint residual
    ! - diag(icol) 写 penalty-like 项
    !
    ! 这里给你“二节点 tie(flag<3)”的示例模板：
    integer :: nc, node1, node2, dof1, dof2, iof1, iof2, irow, icol
    real(prec) :: diagnorm, lmult
    diagnorm = 1.d0   ! 你原来用 diag 的范数做缩放；COO 下可改成固定或另算

    ! do nc = 1, n_constraints
    !   if (constraint_list(nc)%flag < 3) then
    !     icol = ieqs(length_dofs + nc)

    !     node1 = constraint_list(nc)%node1
    !     dof1  = constraint_list(nc)%dof1
    !     iof1  = dof1 + node_list(node1)%dof_index - 1
    !     irow  = ieqs(iof1)

    !     node2 = constraint_list(nc)%node2
    !     dof2  = constraint_list(nc)%dof2
    !     iof2  = dof2 + node_list(node2)%dof_index - 1

    !     ! rhs 对应原 skyline 写法
    !     rhs(irow) = rhs(irow) - lagrange_multipliers(nc)
    !     call A%add_entry(icol, irow,  1.d0)
    !     call A%add_entry(irow, icol,  1.d0)   ! 若你要保持对称；非对称可按需要删掉

    !     irow = ieqs(iof2)
    !     rhs(irow) = rhs(irow) + lagrange_multipliers(nc)
    !     call A%add_entry(icol, irow, -1.d0)
    !     call A%add_entry(irow, icol, -1.d0)

    !     rhs(icol) = dof_total(iof2)+dof_increment(iof2) - &
    !                 dof_total(iof1)-dof_increment(iof1) - lagrange_multipliers(nc)*1.d-12*diagnorm
    !     call A%add_entry(icol, icol, 1.d-12*diagnorm)

    !   else
    !     ! 多节点约束(flag==3)：
    !     ! 你原来调用 user_constraint 得到 element_stiffness/element_residual，
    !     ! 然后用 skyline scatter；这里照单元一样：
    !     ! 1) 构造 gdofs(1:iu)（不过节点来自 nodeset + dof list）
    !     ! 2) add_vec(rhs, gdofs, element_residual)
    !     ! 3) add_block_to_coo(A, gdofs, element_stiffness)
    !     !
    !     ! （这块我建议你把“构造 gdofs”的逻辑单独写一个 build_constraint_gdofs）
    !   end if
    ! end do
  end subroutine
end module AssemblerCOO_mod