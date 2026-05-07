module Mesh
  use Types, only : prec
  implicit none
  private
  public :: mesh_query_element_sizes, mesh_get_element_for_assembly

contains

  subroutine mesh_query_element_sizes(lmn, flag, nnode, ndims, ndofpn, n_props, n_svars)
    integer, intent(in)  :: lmn
    integer, intent(out) :: flag, nnode, ndims, ndofpn, n_props, n_svars
    flag   = 10002
    nnode  = 4
    ndims  = 2
    ndofpn = 2
    n_props = 2
    n_svars = 0
  end subroutine

  subroutine mesh_get_element_for_assembly( lmn, flag, nnode, ndims, ndofpn, &
      coords, u_tot, u_inc, gdofs, n_props, props, n_svars, svars0, svars )
    integer, intent(in) :: lmn
    integer, intent(out) :: flag, nnode, ndims, ndofpn
    real(prec), intent(out) :: coords(:,:), u_tot(:,:), u_inc(:,:)
    integer, intent(out) :: gdofs(:)
    integer, intent(out) :: n_props, n_svars
    real(prec), intent(out) :: props(:), svars0(:), svars(:)

    flag=10002; nnode=4; ndims=2; ndofpn=2; n_props=2; n_svars=0
    if (lmn == 1) then
        coords(1,1)=0.; coords(2,1)=0.; coords(1,2)=1.; coords(2,2)=0.
        coords(1,3)=1.; coords(2,3)=1.; coords(1,4)=0.; coords(2,4)=1.
        gdofs = (/1,2, 3,4, 9,10, 7,8/)
    else
        coords(1,1)=1.; coords(2,1)=0.; coords(1,2)=2.; coords(2,2)=0.
        coords(1,3)=2.; coords(2,3)=1.; coords(1,4)=1.; coords(2,4)=1.
        gdofs = (/3,4, 5,6, 11,12, 9,10/)
    end if
    u_tot=0.; u_inc=0.; props(1)=1000.; props(2)=0.3; svars0=0.; svars=0.
  end subroutine
end module Mesh
