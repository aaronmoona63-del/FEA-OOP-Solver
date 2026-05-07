module Boundaryconditions
    use Types
    !  Data defining a time history of prescribed DOF, generalized forces, or generalized traction
    type history
        sequence
        integer :: index                         ! Index of start of time, value pairs in history_data(:)
        integer :: n_timevals                    ! No. time/value pairs
    end type history
    !  Parameters passed to user-subroutine controlled boundary conditions
    type subroutineparameters
        sequence
        integer :: index                         ! Index of start of parameter list in subroutine_parameters(:)
        integer :: n_parameters                  ! No. subroutine parameters
    end type subroutineparameters
    !  Parameters for constraints
    type constraintparameters
        sequence
        integer :: index                        ! Index of start of parameter list in constraint_parameters(:)
        integer :: n_parameters                 ! No. constraint parameters
    end type constraintparameters
    !  Parameters for node sets
    type nodeset
        sequence
        integer :: index                        ! Index of start of list in node_lists(:)
        integer :: n_nodes                      ! No. nodes in set
    end type nodeset
    !  Parameters for element sets
    type elementset
        sequence
        integer :: index                        ! Index of start of list in element_lists(:)
        integer :: n_elements                   ! No. elements in set
    end type elementset

    type prescribeddof
        sequence
        integer :: flag                          ! Flag specifying nature of DOF - flag=1 prescribed value; flag=2 history, flag=3 user subroutine
        integer :: dof                           ! DOF to prescribe
        integer :: node_number                   ! Node number if single node is constrained, zero otherwise
        integer :: node_set                      ! Node set number if a set of nodes is constrained, zero otherwise
        integer :: history_number                ! Index of DOF history table
        integer :: subroutine_parameter_number   ! Index of user subroutine parameters
        integer :: index_dof_values              ! Index for value of prescribed DOF
        integer :: rate_flag                     ! Set to 1 if DOF rate is to be prescribed; 0 if value is to be prescribed.
    end type prescribeddof

    type prescribedforce
        sequence
        integer :: flag                          ! Flag specifying nature of force - flag=1 prescribed value; flag=2 history, flag=3 user subroutine
        integer :: dof                           ! DOF to apply force to
        integer :: node_number                   ! Node number if single node is loaded, zero otherwise
        integer :: node_set                      ! Node set number if a set of nodes is loaded, zero otherwise
        integer :: history_number                ! Index of DOF history table
        integer :: subroutine_parameter_number   ! Index of user subroutine parameters
        integer :: index_dof_values              ! Index for value of prescribed force
    end type prescribedforce

    type distributedloads
        sequence
        integer :: flag                          ! Flag specifying distributed load type. flag=1, traction value given, flag=2 history and direction give, flag=3 history and normal to element, flag=4 user subroutine
        integer :: element_set                   ! Element set to be loaded
        integer :: face                          ! Face to be loaded
        integer :: history_number                ! History number specifying variation of traction
        integer :: subroutine_parameter_number   !
        integer :: index_dload_values
        integer :: n_dload_values
    end type distributedloads
   
    type constraint
        sequence
        integer :: flag                       ! Flag identifying type of constraint.  flag=1 tie a node pair; flag=2, tie a node set, flag=3 general MPC
        integer :: node1                      ! First node (or node set) in constraint
        integer :: dof1                       ! Dof for first node (not used if flag=3)
        integer :: node2                      ! Second node (or node set, or DOF list) in constraint
        integer :: dof2                       ! Dof for second node (not used if flag=3)
        integer :: index_parameters           ! Index of parameters associated with constraint in parameter list
    end type constraint


   
    integer :: n_histories                   ! No. load histories
    integer :: n_subroutine_parameters       ! No. lists of subroutine parameters for distributed loads or DOF
    integer :: n_constraint_parameters       ! No. lists of parameters for constraints
    integer :: n_nodesets                    ! No. node sets
    integer :: n_elementsets                 ! No. element sets
    integer :: n_prescribeddof               ! No. prescribed DOF
    integer :: n_prescribedforces            ! No. prescribed forces
    integer :: n_distributedloads            ! No. distributed loads or fluxes
    integer :: n_constraints                 ! No. constraints
   
    integer :: length_node_lists
    integer :: length_element_lists
    integer :: length_history_data
    integer :: length_subroutine_parameters
    integer :: length_dof_values
    integer :: length_dload_values
    integer :: length_constraint_parameters
   
    integer, allocatable :: node_lists(:)
    integer, allocatable :: element_lists(:)
   
    real (prec), allocatable :: history_data(:,:)
    real (prec), allocatable :: subroutine_parameters(:)
    real (prec), allocatable :: dload_values(:)
    real (prec), allocatable :: dof_values(:)
    real (prec), allocatable :: constraint_parameters(:)
   
    real (prec), allocatable :: lagrange_multipliers(:)
   
    character (len=100), allocatable :: elementset_namelist(:)
    character (len=100), allocatable :: nodeset_namelist(:)
    character (len=100), allocatable :: history_namelist(:)
    character (len=100), allocatable :: subroutineparameter_namelist(:)
    character (len=100), allocatable :: constraintparameter_namelist(:)
     
    type (history), allocatable :: history_list(:)
    type (subroutineparameters), allocatable :: subroutineparameter_list(:)
    type (constraintparameters), allocatable :: constraintparameter_list(:)
    type (elementset), allocatable :: elementset_list(:)
    type (nodeset), allocatable :: nodeset_list(:)
    type (prescribeddof), allocatable :: prescribeddof_list(:)
    type (prescribedforce), allocatable :: prescribedforce_list(:)
    type (distributedloads), allocatable :: distributedload_list(:)
    type (constraint), allocatable :: constraint_list(:)
   
contains
  subroutine initialize_Boundaryconditions_vars()
    n_histories = 0
    n_subroutine_parameters = 0
    n_constraint_parameters = 0
    n_nodesets = 0
    n_elementsets = 0
      
    n_prescribeddof = 0
    n_prescribedforces = 0
    n_distributedloads = 0
    n_constraints = 0
  
    length_element_lists = 0
    length_node_lists = 0
    length_history_data = 0
    length_subroutine_parameters = 0
    length_dof_values = 0
    length_dload_values = 0
    length_constraint_parameters = 0
    
    ! Solution data
    if (allocated(lagrange_multipliers)) lagrange_multipliers = 0.d0  
  end subroutine initialize_Boundaryconditions_vars


  function find_name_index(name,len,name_list,n_names)
    !通过在给定的名称列表中查找特定名称来返回其索引

    use Types
    !use ParamIO
    implicit none

    character ( len = 100 )     :: name   !待查找的名称，最大长度为100。
    integer                     :: len   !name的实际长度
    integer                     :: n_names   !名称列表中的名称总数。
    character ( len = 100 )     :: name_list(n_names) !一个包含多个名称的数组，其大小由 n_names 确定

    integer                     :: find_name_index

    integer:: IOW=6

    ! Local Variables
    integer :: n
    
    do n = 1, n_names
      if (len == len_trim(name_list(n)) ) then
        if (strcmp(name, name_list(n), len)) then
          find_name_index = n
          return
        end if
      endif
    end do

    write(IOW,*) ' *** Error detected in input file ***'
    write(IOW,*) ' The named history, node set, element set, or parameter list ',name
    write(IOW,*) ' was not defined '
    stop

  end function find_name_index
 
  subroutine interpolate_history_table(history,nhist,time,dofvalue)
    use Types

    integer,  intent( in )         :: nhist

    real( prec ), intent( in )     :: history(2,nhist)
    real( prec ), intent( in )     :: time

    real( prec ), intent( out )    :: dofvalue

    !Local variables
    integer :: klo,khi,k

    !   Subroutine to interpolate a load history table
    !
    !  Find positions in history table within which to interpolate

    if (nhist == 1) then
      dofvalue = history(2,1)
      return
    endif

    if (time <= history(1,1) ) then
      dofvalue = history(2,1)
    else if (time >= history(1,nhist) ) then
      dofvalue = history(2,nhist)
    else

      klo = 1
      khi = nhist

      do while ( .true. )
        if ( khi - klo>1 ) then
            k = (khi + klo)/2
            if ( history(1, k)>time ) then
                khi = k
            else
                klo = k
            end if
            cycle
        else
            exit
        end if
      end do

      if (history(1,khi) == history(1,klo)) then
        dofvalue = 0.5D0*(history(2,khi) + history(2,klo))
        return
      endif
      dofvalue =( (history(1,khi)-time)*history(2,klo) +  &
        (time-history(1,klo))*history(2,khi)  )/ &
        (history(1,khi)-history(1,klo))
    endif
    return
  end subroutine interpolate_history_table

   logical function strcmp(a, b, n)
    implicit none
    character(len=*), intent(in) :: a, b
    integer, intent(in) :: n
    strcmp = (a(1:n) == b(1:n))
  end function strcmp

  subroutine get_elem_list(elset_id, elems)
    integer, intent(in) :: elset_id
    integer, allocatable, intent(out) :: elems(:)

    ! 从 elementset_list / element_lists 里取
    allocate(elems(elementset_list(elset_id)%n_elements))
    elems = element_lists( elementset_list(elset_id)%index : &
                            elementset_list(elset_id)%index + size(elems) - 1 )
  end subroutine
  
  !=======================================================================
  ! Subroutine: get_gdofs_from_face_continuous
  !
  ! Purpose:
  !   Construct the global degree-of-freedom (DOF) index list corresponding
  !   to a face of an element, assuming a *continuous global DOF numbering*
  !   scheme.
  !
  ! Description:
  !   This routine maps local face nodes and their per-node DOFs to a flat
  !   list of global DOF indices (gdofs), suitable for assembling element-
  !   level boundary contributions (e.g. traction forces) into the global
  !   system.
  !
  !   The global DOF numbering is assumed to be:
  !
  !       global_dof = (node_id - 1) * ndofpn + local_dof
  !
  !   where:
  !     - node_id   : global node number (1-based)
  !     - ndofpn    : number of DOFs per node
  !     - local_dof : local DOF index at the node (1-based)
  !
  !   This mapping corresponds to a *dense, continuous* DOF layout without
  !   equation reordering, constraints elimination, or DOF compression.
  !
  ! Typical Use Case:
  !   - Assembly of boundary condition contributions (Neumann / traction)
  !   - Assembly of element stiffness/residual for simple test problems
  !   - Early-stage solver development before introducing equation maps
  !
  ! Inputs:
  !   face_nodes(:)
  !     Global node numbers belonging to the element face, ordered according
  !     to the element's face connectivity.
  !
  !   ndofpn
  !     Number of degrees of freedom per node.
  !
  ! Outputs:
  !   gdofs(:)
  !     Allocated array of length size(face_nodes) * ndofpn containing the
  !     global DOF indices corresponding to all DOFs on the face.
  !
  !     The ordering is:
  !       [ node1_dof1, node1_dof2, ..., node1_dof_ndofpn,
  !         node2_dof1, node2_dof2, ..., node2_dof_ndofpn,
  !         ... ]
  !
  ! Notes:
  !   - This routine assumes *no constrained DOFs* and *no equation mapping*.
  !   - For constrained systems or reordered equations, this routine should
  !     be replaced by a mapping through an equation/DOF map.
  !   - Intended to be injected via BCExecutionContext or assembler callbacks,
  !     not hard-coded inside BoundaryLayer logic.
  !
  !=======================================================================
  subroutine get_gdofs_from_face_continuous(face_nodes, ndofpn, gdofs)
    integer, intent(in) :: face_nodes(:), ndofpn
    integer, allocatable, intent(out) :: gdofs(:)

    integer :: a, d, pos
    allocate(gdofs(size(face_nodes)*ndofpn))

    pos = 0
    do a = 1, size(face_nodes)
      do d = 1, ndofpn
        pos = pos + 1
        gdofs(pos) = (face_nodes(a)-1)*ndofpn + d
      end do
    end do
  end subroutine




    !------------------------------------------------------------
  ! 下面这些 bc_* 例程：你应该放进 Boundaryconditions 的“封装层”
  ! 我这里只给接口名，避免你现在把 legacy 数组直接暴露给 apply
  !------------------------------------------------------------
  integer function get_n_distributedloads(bc) result(n)
    class(*), intent(in) :: bc
    stop "Implement get_n_distributedloads in your BC wrapper"
  end function

  subroutine bc_get_distributedload(bc, load, flag, elset, ifac, ntract, traction, scale)
    class(*), intent(in) :: bc
    integer, intent(in)  :: load
    integer, intent(out) :: flag, elset, ifac, ntract
    real(prec), intent(out) :: traction(3)
    real(prec), intent(out) :: scale
    stop "Implement bc_get_distributedload in your BC wrapper"
  end subroutine

  logical function bc_is_uel_element(bc, mesh, lmn) result(isuel)
    class(*), intent(in) :: bc, mesh
    integer, intent(in) :: lmn
    isuel = .false.
  end function



  subroutine bc_eval_distributedload(load, time, dtime, flag, elset, ifac, ntract, traction)
  use Types, only : prec
  !use Boundaryconditions   ! 直接用你现有全局数组
  implicit none
  integer, intent(in) :: load
  real(prec), intent(in) :: time, dtime
  integer, intent(out) :: flag, elset, ifac, ntract
  real(prec), intent(out) :: traction(3)

  integer :: iof2, nhist
  real(prec) :: scale, norm2

  traction = 0.0_prec

  flag  = distributedload_list(load)%flag
  elset = distributedload_list(load)%element_set
  ifac  = distributedload_list(load)%face

  if (flag < 3) then
    ntract = distributedload_list(load)%n_dload_values
    if (ntract > 3) stop "bc_eval_distributedload: ntract>3 not supported here"
    traction(1:ntract) = dload_values(distributedload_list(load)%index_dload_values : &
                                      distributedload_list(load)%index_dload_values + ntract - 1)
  else
    ! flag==3: normal traction, only magnitude provided later by history
    ntract = 1
    traction(1) = 1.0_prec
  end if

  if (flag == 2) then
    norm2 = dot_product(traction(1:ntract), traction(1:ntract))
    if (norm2 <= 0.0_prec) stop "bc_eval_distributedload: direction norm = 0"
    traction(1:ntract) = traction(1:ntract) / sqrt(norm2)
  end if

  if (flag > 1 .and. distributedload_list(load)%history_number > 0) then
    iof2  = history_list(distributedload_list(load)%history_number)%index
    nhist = history_list(distributedload_list(load)%history_number)%n_timevals
    call interpolate_history_table(history_data(1,iof2), nhist, time + dtime, scale)
    traction(1:ntract) = traction(1:ntract) * scale
  end if
  end subroutine



end module