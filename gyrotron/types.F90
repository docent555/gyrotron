module types
   use, intrinsic :: iso_c_binding

   type Input
      integer Ne
      real(c_double), allocatable :: ZAxis(:), TAxis(:)
      real(c_double) Delta, Ic, dz, dt
      integer(c_int) INTERVALT, INTERVALZ, OUTNz, Nz, Nt
      complex(c_double_complex), allocatable :: InitialField(:)
   end type Input

   type Output
      real(c_double), allocatable :: ZAxis(:), TAxis(:)
      complex(c_double_complex), allocatable :: OUTB(:, :), OUTCu(:, :)
   end type Output

   complex(c_double_complex) :: Im1 = (0.0d0, 1.0d0)
   !real(c_double), parameter :: pi = 2.0d0*dacos(0.0d0)
   real(c_double) :: pi = dacos(-1.0d0)

end module types
