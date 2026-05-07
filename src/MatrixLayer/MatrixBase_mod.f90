module MatrixBase_mod
  implicit none
  private
  public :: MatrixBase

  type, abstract :: MatrixBase
     integer :: n = 0
     logical :: is_symmetric = .false.
   contains
     procedure(init_iface),          deferred :: init
     procedure(get_size_iface),      deferred :: get_size
     procedure(print_iface),         deferred :: print
     procedure(to_skyline_iface),    deferred :: to_skyline
  end type MatrixBase

  ! ============================
  ! Abstract Interfaces
  ! ============================
  abstract interface

    subroutine init_iface(self, n, is_symmetric)
      import :: MatrixBase
      class(MatrixBase), intent(inout) :: self
      integer, intent(in) :: n
      logical, intent(in) :: is_symmetric
    end subroutine init_iface

    function get_size_iface(self) result(n)
      import :: MatrixBase
      class(MatrixBase), intent(in) :: self
      integer :: n
    end function get_size_iface

    subroutine print_iface(self)
      import :: MatrixBase
      class(MatrixBase), intent(in) :: self
    end subroutine print_iface

    subroutine to_skyline_iface(self, al, au, ad, jp)
      import :: MatrixBase
      class(MatrixBase), intent(in) :: self
      double precision, allocatable, intent(out) :: al(:), au(:), ad(:)
      integer, allocatable, intent(out) :: jp(:)
    end subroutine to_skyline_iface

  end interface

end module MatrixBase_mod