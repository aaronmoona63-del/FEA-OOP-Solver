module SparseMatrixBase_mod
    implicit none
    private
    public :: SparseMatrixBase

    !=========================================================
    !  抽象基类：稀疏矩阵基类（方阵）
    !
    !  责任非常简单：
    !    - 记录阶数 n
    !    - 提供统一接口：
    !        init(n)
    !        to_dense()
    !        spmv(x, y)
    !        get_nnz()
    !
    !  其它操作（add_entry、格式转换等）由各具体子类
    !  或者独立的转换模块 SparseConvert_mod 实现。
    !=========================================================
    type, abstract :: SparseMatrixBase
        integer :: n = 0      ! 方阵 n × n
    contains
        ! 初始化矩阵尺寸
        procedure(init_iface),    deferred :: init

        ! 转换成 dense（主要用来 debug / 验证）
        procedure(to_dense_iface), deferred :: to_dense

        ! 稀疏矩阵–向量乘法 y = A*x
        procedure(spmv_iface),    deferred :: spmv

        ! 返回非零元个数（统计/监控用）
        procedure(get_nnz_iface), deferred :: get_nnz

        procedure(print_iface), deferred :: print
    end type SparseMatrixBase


    !=========================================================
    !  抽象接口定义
    !=========================================================
    abstract interface

        !-------------------------------------
        ! init(self, n)
        !  初始化矩阵阶数 n×n
        !-------------------------------------
        subroutine init_iface(self, n)
            import :: SparseMatrixBase
            class(SparseMatrixBase), intent(inout) :: self
            integer, intent(in) :: n
        end subroutine init_iface


        !-------------------------------------
        ! to_dense(self) -> A(n,n)
        !-------------------------------------
        function to_dense_iface(self) result(A)
            import :: SparseMatrixBase
            class(SparseMatrixBase), intent(in) :: self
            real(8), allocatable :: A(:,:)
        end function to_dense_iface


        !-------------------------------------
        ! spmv(self, x, y)
        !   y = A * x
        !-------------------------------------
        subroutine spmv_iface(self, x, y)
            import :: SparseMatrixBase
            class(SparseMatrixBase), intent(in) :: self
            real(8), intent(in)  :: x(:)
            real(8), intent(out) :: y(:)
        end subroutine spmv_iface


        !-------------------------------------
        ! get_nnz(self) -> nz
        !-------------------------------------
        function get_nnz_iface(self) result(nz)
            import :: SparseMatrixBase
            class(SparseMatrixBase), intent(in) :: self
            integer :: nz
        end function get_nnz_iface

        !-------------------------------------
        ! print(self)
        !   统一调试接口（子类必须实现）
        !-------------------------------------
        subroutine print_iface(self)
            import :: SparseMatrixBase
            class(SparseMatrixBase), intent(in) :: self
        end subroutine print_iface

    end interface

end module SparseMatrixBase_mod
