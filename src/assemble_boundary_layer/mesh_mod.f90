module Mesh
  use Types

  !  Data type for nodes
  type node
  sequence
  integer :: flag                          ! Integer identifier
  integer :: coord_index                   ! Index of first coordinate in coordinate array
  integer :: n_coords                      ! Total no. coordinates for the node
  integer :: dof_index                     ! Index of first DOF in dof array
  integer :: n_dof                         ! Total no. of DOF for node
  integer :: displacement_map_index        ! Index of displacement node map
  integer :: n_displacements               ! No. displacement DOFs
  end type node

  !  Data type for elements
  type element
  sequence
  integer :: flag                          ! Integer identifier for element
  integer :: connect_index                 ! Index of first node on element in connectivity(:)
  integer :: n_nodes                       ! No. nodes on the element
  integer :: state_index                   ! Index of first state variable in element_state_variables(:)
  integer :: n_states                      ! No. state variables
  integer :: element_property_index        ! Index of first element property in element_properties(:)
  integer :: n_element_properties          ! No. element properties
  integer :: density_index                 ! Index of density value
  integer :: int_element_property_index    ! Index of integer element properties (used for ABAQUS subroutines)
  integer :: n_int_element_properties      ! No. integer valued element properties (used for ABAQUS subroutines)
  integer :: material_index                ! Index of material assigned to element
  end type element

  !  Data type for material
  type material
  sequence
  integer :: prop_index                    ! Index of first material property in property array
  integer :: n_properties                  ! No. properties for this material
  integer :: n_states                      ! No. history dependent state variables for this material
  end type material

  !  Data type for zone
  type zone
  sequence
  integer :: start_element                 ! First element in a zone
  integer :: end_element                   ! Last element in a zone
  end type zone

  !  Data type for storing ABAQUS UEL distributed loads
  type abq_uel_bc
  sequence
  integer :: mdload                        ! Number of BCs applied to abaqus UEL
  integer :: mag_index                     ! Index to arrays storing magnitude and type of BC applied to abaqus UEL
  end type abq_uel_bc


  integer :: n_zones                          ! 区域总数
  integer :: n_nodes                          ! 节点总数
  integer :: n_elements                       ! 单元总数
  integer :: n_materials                      ! Total number of materials (used in ABAQUS UMAT and VUMAT)
  integer :: length_coords                    ! Length of coordinate array
  integer :: length_dofs                      ! Length of nodal DOF array
  integer :: length_connectivity              ! Length of connectivity array
  integer :: length_element_properties        ! Length of element property array
  integer :: length_int_element_properties    ! Length of integer valued element property array
  integer :: length_material_properties       ! Length of material property array
  integer :: length_densities                 ! Length of density array
  integer :: length_state_variables           ! Length of state variable array
  integer :: length_displacement_map          ! Length of array mapping nodal DOF to displacements
  integer :: n_mesh_parameters                ! Parameeters controlling a user-subrouine generated mesh
  integer :: length_abq_dlmag_array           ! Array dimension for abaqus uel boundary conditions

  integer, allocatable :: displacement_map(:)       ! Array storing mapping of nodal DOF to displacements
  integer, allocatable :: connectivity(:)           ! Array storing element connectivity
  integer, allocatable :: int_element_properties(:) ! Array storing integer valued element properties

  integer, allocatable :: abq_uel_bc_typ(:)          ! Array storing boundary condition type flags for abaqus UEL bcs
  integer, allocatable :: abq_MCRD(:)                         ! Abaqus uel MCRD parameter

  real (prec) :: nodal_force_norm             ! Norm of nodal generalized force
  real (prec) :: unbalanced_force_norm        ! Norm of out-of-balance forces
  real (prec) :: correction_norm              ! Norm of solution correction

  real (prec), allocatable :: element_properties(:)            ! List of element properties
  real (prec), allocatable :: material_properties(:)           ! List of material properties
  real (prec), allocatable :: densities(:)                     ! List of density values for zones
  real (prec), allocatable :: initial_state_variables(:)       ! Element state variables at the start of a time increment
  real (prec), allocatable :: updated_state_variables(:)       ! Element state variables at the end of a time increment
  real (prec), allocatable :: coords(:)                        ! List of nodal coordinates
  real (prec), allocatable :: dof_total(:)                     ! List of accumulated DOF
  real (prec), allocatable :: dof_increment(:)                 ! List of increment in DOF
  real (prec), allocatable :: energy(:)                        ! Energy (used by ABAQUS UMAT and UEL)
  real (prec), allocatable :: velocity(:)                      ! Velocity (for explicit dynamics)
  real (prec), allocatable :: acceleration(:)                  ! Acceleration (for explicit dynamics)
  real (prec), allocatable :: lumped_mass(:)                   ! Lumped mass matrix (for explicit dynamics)
  real (prec), allocatable :: rforce(:)
  real (prec), allocatable :: mesh_subroutine_parameters(:)    ! List of parameters controlling a user-subroutine generated mesh


  real (prec), allocatable :: abq_uel_bc_mag(:)                      ! Array storing magnitudes of BCs for abaqus UEL
  real (prec), allocatable :: abq_uel_bc_dmag(:)                     ! Array storing magnitudes of increment in BCs for abaqus UEL

  type (node),       save, allocatable :: node_list(:)
  type (element),    save, allocatable :: element_list(:)
  type (zone),       save, allocatable :: zone_list(:)
  type (material),   save, allocatable :: material_list(:)
  type (abq_uel_bc), allocatable :: abq_uel_bc_list(:)

  character (len=100), allocatable :: zone_namelist(:)         ! Names of zones in mesh
  character (len=100), allocatable :: material_namelist(:)     ! Names of materials

  logical, allocatable :: element_deleted(:)             ! Flag listing deleted elements in an explicit dynamic simulation

contains

  subroutine initialize_Mesh_vars()
    n_zones = 0
    n_nodes = 0
    n_elements = 0
    n_materials = 0
    length_coords = 0
    length_dofs = 0
    length_displacement_map = 0

    length_element_properties = 0
    length_int_element_properties = 0
    length_material_properties = 0

    length_densities = 0
    length_state_variables = 0
    length_connectivity = 0

    n_mesh_parameters = 0

    if (allocated(zone_namelist)) zone_namelist = 'Unnamed Zone'
    if (allocated(node_list)) node_list(1:n_nodes)%displacement_map_index = 0
    if (allocated(node_list)) node_list(1:n_nodes)%n_displacements = 0


    ! Solution data
    if (allocated(initial_state_variables) ) initial_state_variables = 0.d0
    if (allocated(updated_state_variables) ) updated_state_variables = 0.d0
    if (allocated(dof_total) ) dof_total = 0.d0
    if (allocated(dof_increment) )dof_increment = 0.d0
    if (allocated(rforce) ) rforce = 0.d0
    if (allocated(element_properties) ) element_properties = 0.d0
    if (allocated(int_element_properties) ) int_element_properties = 0

  end subroutine initialize_Mesh_vars
    
  subroutine extract_element_data(lmn,flag,n_nodes,node_list,n_properties,properties,n_state_variables,initial_svars,updated_svars)

    use Types
    use ParamIO

    implicit none

    integer, intent( in )  :: lmn

    integer, intent( out ) :: node_list(*)
    integer, intent( out ) :: flag
    integer, intent( out ) :: n_nodes
    integer, intent( out ) :: n_properties
    integer, intent( out ) :: n_state_variables

    real (prec), intent( out ) ::  properties(*)
    real (prec), intent( out ) ::  initial_svars(*)
    real (prec), intent( out ) ::  updated_svars(*)

    !  Function to extract data for element number lmn
    !
    if (lmn>n_elements) then
      write(IOW,*) ' Error in subroutine extract_element_data '
      write(IOW,*) ' Element number ',lmn,' exceeds the number of elements in the mesh '
      stop
    endif

    flag = element_list(lmn)%flag
    n_nodes = element_list(lmn)%n_nodes
    n_state_variables = element_list(lmn)%n_states
    n_properties = element_list(lmn)%n_element_properties

    node_list(1:n_nodes) = connectivity(element_list(lmn)%connect_index:element_list(lmn)%connect_index+n_nodes-1)
    if (n_properties>0) then
      properties(1:n_properties) = element_properties(element_list(lmn)%element_property_index: &
        element_list(lmn)%element_property_index+n_properties-1)
    endif
    if (n_state_variables>0) then
      initial_svars(1:n_state_variables) = initial_state_variables(element_list(lmn)%state_index: &
        element_list(lmn)%state_index+n_state_variables-1)
      updated_svars(1:n_state_variables) = updated_state_variables(element_list(lmn)%state_index: &
        element_list(lmn)%state_index+n_state_variables-1)
    endif

  end subroutine extract_element_data

  subroutine extract_node_data(nn,flag,n_coords,nodal_coords,n_dof,nodal_dof_increment,nodal_dof_total)

    use Types
    use ParamIO

    implicit none

    integer, intent( in )  ::  nn

    integer, intent( out ) :: n_coords
    integer, intent( out ) :: flag
    integer, intent( out ) :: n_dof

    real (prec), intent( out ) ::  nodal_coords(*)
    real (prec), intent( out ) ::  nodal_dof_increment(*)
    real (prec), intent( out ) ::  nodal_dof_total(*)

    !  Function to extract data for node number nn
    !
    if (nn>n_nodes) then
      write(IOW,*) ' Error in subroutine exatract_node_data '
      write(IOW,*) ' Node number ',nn,' exceeds the number of nodes in the mesh '
      stop
    endif

    flag = node_list(nn)%flag
    n_coords = node_list(nn)%n_coords
    n_dof = node_list(nn)%n_dof


    nodal_coords(1:n_coords) = coords(node_list(nn)%coord_index:node_list(nn)%coord_index+n_coords-1)
    nodal_dof_increment(1:n_dof) = dof_increment(node_list(nn)%dof_index: &
      node_list(nn)%dof_index + n_dof-1)
    nodal_dof_total(1:n_dof) = dof_total(node_list(nn)%dof_index: &
      node_list(nn)%dof_index + n_dof-1)

  end subroutine extract_node_data

  subroutine get_element_data(lmn, flag, n_nodes, nodes, &
                              n_properties, properties, &
                              n_state_variables, initial_svars, updated_svars)

    integer, intent(in)  :: lmn
    integer, intent(out) :: flag, n_nodes, n_properties, n_state_variables
    integer, intent(out) :: nodes(:)
    real(prec), intent(out) :: properties(:)
    real(prec), intent(out) :: initial_svars(:), updated_svars(:)

    if (lmn > n_elements) then
      write(IOW,*) 'Error: element index out of range:', lmn
      stop
    endif

    flag = element_list(lmn)%flag
    n_nodes = element_list(lmn)%n_nodes
    n_state_variables = element_list(lmn)%n_states
    n_properties = element_list(lmn)%n_element_properties

    nodes(1:n_nodes) = connectivity( &
      element_list(lmn)%connect_index : &
      element_list(lmn)%connect_index + n_nodes - 1 )

    if (n_properties > 0) then
      properties(1:n_properties) = element_properties( &
        element_list(lmn)%element_property_index : &
        element_list(lmn)%element_property_index + n_properties - 1 )
    endif

    if (n_state_variables > 0) then
      initial_svars(1:n_state_variables) = initial_state_variables( &
        element_list(lmn)%state_index : &
        element_list(lmn)%state_index + n_state_variables - 1 )

      updated_svars(1:n_state_variables) = updated_state_variables( &
        element_list(lmn)%state_index : &
        element_list(lmn)%state_index + n_state_variables - 1 )
    endif

  end subroutine get_element_data

  subroutine get_node_data(nn, flag, n_coords, nodal_coords, &
                           n_dof, nodal_dof_increment, nodal_dof_total)

    integer, intent(in)  :: nn
    integer, intent(out) :: flag, n_coords, n_dof
    real(prec), intent(out) :: nodal_coords(:)
    real(prec), intent(out) :: nodal_dof_increment(:)
    real(prec), intent(out) :: nodal_dof_total(:)

    if (nn > n_nodes) then
      write(IOW,*) 'Error: node index out of range:', nn
      stop
    endif

    flag = node_list(nn)%flag
    n_coords = node_list(nn)%n_coords
    n_dof = node_list(nn)%n_dof

    nodal_coords(1:n_coords) = coords( &
      node_list(nn)%coord_index : &
      node_list(nn)%coord_index + n_coords - 1 )

    nodal_dof_increment(1:n_dof) = dof_increment( &
      node_list(nn)%dof_index : &
      node_list(nn)%dof_index + n_dof - 1 )

    nodal_dof_total(1:n_dof) = dof_total( &
      node_list(nn)%dof_index : &
      node_list(nn)%dof_index + n_dof - 1 )

  end subroutine get_node_data

  subroutine mesh_get_element_for_assembly( lmn,          &
      flag,                                               &
      nnode, ndims, ndof_per_node,                        &
      coords_elem,                                        &
      u_elem_total,                                       &
      u_elem_increment,                                   &
      gdofs,                                              &
      n_props, props,                                     &
      n_svars, svars0, svars )

    use Types
    use ParamIO
    implicit none

    !------------------------
    ! Input
    !------------------------
    integer, intent(in) :: lmn

    !------------------------
    ! Output: topology
    !------------------------
    integer, intent(out) :: flag
    integer, intent(out) :: nnode
    integer, intent(out) :: ndims
    integer, intent(out) :: ndof_per_node

    !------------------------
    ! Output: geometry & dof
    !------------------------
    real(prec), intent(out) :: coords_elem(:,:)        ! (ndims, nnode)
    real(prec), intent(out) :: u_elem_total(:,:)       ! (ndof_per_node, nnode)
    real(prec), intent(out) :: u_elem_increment(:,:)   ! (ndof_per_node, nnode)
    integer,    intent(out) :: gdofs(:)                 ! (ndof_per_node*nnode)

    !------------------------
    ! Output: material & state
    !------------------------
    integer, intent(out) :: n_props, n_svars
    real(prec), intent(out) :: props(:)
    real(prec), intent(out) :: svars0(:)
    real(prec), intent(out) :: svars(:)

    !------------------------
    ! Local variables
    !------------------------
    integer :: j, d, pos
    integer :: node_id
    integer :: prop_idx, svar_idx

    !------------------------
    ! Sanity check
    !------------------------
    if (lmn > n_elements) then
      write(IOW,*) 'Error: element index out of range:', lmn
      stop
    end if

    !------------------------
    ! Basic element info
    !------------------------
    flag  = element_list(lmn)%flag
    nnode = element_list(lmn)%n_nodes

    call get_element_type_info(flag, ndims, ndof_per_node)

    !------------------------
    ! Fill node-based data
    !------------------------
    pos = 0
    do j = 1, nnode
      node_id = connectivity( element_list(lmn)%connect_index + j - 1 )

      ! Coordinates
      coords_elem(:, j) = coords( &
        node_list(node_id)%coord_index : &
        node_list(node_id)%coord_index + ndims - 1 )

      ! DOF total
      u_elem_total(:, j) = dof_total( &
        node_list(node_id)%dof_index : &
        node_list(node_id)%dof_index + ndof_per_node - 1 )

      ! DOF increment
      u_elem_increment(:, j) = dof_increment( &
        node_list(node_id)%dof_index : &
        node_list(node_id)%dof_index + ndof_per_node - 1 )

      ! Global DOF indices (COO row/col numbers)
      do d = 1, ndof_per_node
        pos = pos + 1
        gdofs(pos) = node_list(node_id)%dof_index + d - 1
      end do
    end do

    !------------------------
    ! Element properties
    !------------------------
    n_props = element_list(lmn)%n_element_properties
    if (n_props > 0) then
      prop_idx = element_list(lmn)%element_property_index
      props(1:n_props) = element_properties( prop_idx : prop_idx + n_props - 1 )
    end if

    !------------------------
    ! State variables
    !------------------------
    n_svars = element_list(lmn)%n_states
    if (n_svars > 0) then
      svar_idx = element_list(lmn)%state_index
      svars0(1:n_svars) = initial_state_variables( svar_idx : svar_idx + n_svars - 1 )
      svars (1:n_svars) = updated_state_variables( svar_idx : svar_idx + n_svars - 1 )
    end if

  end subroutine mesh_get_element_for_assembly



  subroutine get_element_type_info(flag, ndims, ndof_per_node)
    implicit none
    integer, intent(in)  :: flag
    integer, intent(out) :: ndims
    integer, intent(out) :: ndof_per_node

    select case (flag)
    case (10002)     ! 2D continuum
      ndims = 2
      ndof_per_node = 2

    case (10003)     ! 3D continuum
      ndims = 3
      ndof_per_node = 3

    case default
      write(*,*) 'Unknown element flag:', flag
      stop
    end select
  end subroutine


  subroutine mesh_query_element_sizes( lmn,      &
      flag,                                      &
      nnode, ndims, ndof_per_node,               &
      n_props, n_svars )

    use Types
    use ParamIO
    implicit none

    integer, intent(in)  :: lmn

    integer, intent(out) :: flag
    integer, intent(out) :: nnode
    integer, intent(out) :: ndims
    integer, intent(out) :: ndof_per_node
    integer, intent(out) :: n_props
    integer, intent(out) :: n_svars

    ! --- bounds check ---
    if (lmn > n_elements) then
      write(IOW,*) 'mesh_query_element_sizes: element index out of range:', lmn
      stop
    end if

    ! --- topology ---
    flag  = element_list(lmn)%flag
    nnode = element_list(lmn)%n_nodes

    ! --- element-type metadata ---
    call get_element_type_info(flag, ndims, ndof_per_node)

    ! --- material / state sizes ---
    n_props = element_list(lmn)%n_element_properties
    n_svars = element_list(lmn)%n_states

  end subroutine mesh_query_element_sizes

end module
end module Mesh