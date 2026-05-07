module ElementDispatch_mod
  use Types, only : prec
  use Element_Utilities
  implicit none
  private
  public :: element_dispatch

contains

  subroutine element_dispatch(flag, coords, u_tot, u_inc, &
                              props, svars0, svars, Ke, fe, fail)

    integer, intent(in) :: flag
    real(prec), intent(in)  :: coords(:,:), u_tot(:,:), u_inc(:,:)
    real(prec), intent(in)  :: props(:)
    real(prec), intent(in)  :: svars0(:)
    real(prec), intent(inout) :: svars(:)
    real(prec), intent(out) :: Ke(:,:), fe(:)
    logical, intent(out) :: fail

    fail = .false.
    Ke = 0.0_prec
    fe = 0.0_prec

    select case (flag)
    case (10002)   ! 2D continuum
      call continuum_kernel_2D(coords, u_tot, u_inc, props, svars0, svars, Ke, fe)
    case default
      write(*,*) 'element_dispatch: unsupported element flag = ', flag
      fail = .true.
    end select
  end subroutine element_dispatch

  subroutine continuum_kernel_2D(coords, u_tot, u_inc, props, svars0, svars, Ke, fe)
    use Types, only : prec
    implicit none

    real(prec), intent(in)    :: coords(:,:), u_tot(:,:), u_inc(:,:), props(:)
    real(prec), intent(in)    :: svars0(:)
    real(prec), intent(inout) :: svars(:)
    real(prec), intent(out)   :: Ke(:,:), fe(:)

    ! Local variables
    integer :: i, k, a
    real(prec) :: E, nu
    real(prec) :: D(3,3)
    real(prec) :: xi_list(2,4), w_list(4)
    real(prec) :: N(4), dN_dxi(4,2), dN_dx(4,2)
    real(prec) :: J(2,2), J_inv(2,2), detJ
    real(prec) :: B(3,8)
    real(prec) :: xi_p(2), w_p
    real(prec) :: DB(3,8) 

    ! 1. Material properties
    E  = props(1)
    nu = props(2)

    ! 2. Construct D matrix (Plane Stress)
    D = 0.0_prec
    D(1,1) = 1.0_prec
    D(1,2) = nu
    D(2,1) = nu
    D(2,2) = 1.0_prec
    D(3,3) = (1.0_prec - nu) / 2.0_prec
    D = D * (E / (1.0_prec - nu*nu))

    ! 3. Initialize integration points 
    call initialize_integration_points(4, 4, xi_list, w_list)

    Ke = 0.0_prec
    fe = 0.0_prec 

    ! 4. Integration loop
    do k = 1, 4
        xi_p = xi_list(:,k)
        w_p  = w_list(k)
        call calculate_shapefunctions(xi_p, 4, N, dN_dxi)
        J = matmul(coords, dN_dxi)
        call invert_small(J, J_inv, detJ)
        dN_dx = matmul(dN_dxi, J_inv)
        B = 0.0_prec
        do a = 1, 4
            B(1, 2*a-1) = dN_dx(a,1)
            B(2, 2*a)   = dN_dx(a,2)
            B(3, 2*a-1) = dN_dx(a,2)
            B(3, 2*a)   = dN_dx(a,1)
        end do
        DB = matmul(D, B)
        Ke = Ke + matmul(transpose(B), DB) * detJ * w_p
    end do
  end subroutine
end module ElementDispatch_mod
