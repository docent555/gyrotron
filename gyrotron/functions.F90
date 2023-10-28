module functions
   use, intrinsic :: iso_c_binding

   complex(c_double_complex), private :: C0, CR = 0, C2
   real(c_double), private :: SQR2 = dsqrt(2.0D0), SQR2M2 = 2.828427124746190, SQR2D2 = 0.707106781186548

   complex(c_double_complex), private :: IR, IROldPart, WNz, WNzm1, &
                                         WR_step_part, D_end_part, C0mSQRDZdDTm2, C0m2d3dDT, &
                                         DZmC0m2d3dDT, DZmC0d3dDT, C2mSQRDTm4d3
   real(c_double), private :: DeltaT, SQRDT, SQRDZ, SQRDTm4d3, SQRDTm2d3, &
                              DeltaZd3, DeltaZd6
   integer(c_int), private :: i, j, err_alloc = 1, step, jout = 1, Nzm1, inloopnum

   complex(c_double_complex), private :: Im1 = (0.0d0, 1.0d0)
   real(c_double), private :: pi = dacos(-1.0d0) !real(c_double), parameter :: pi = 2.0d0*dacos(0.0d0)

   complex(c_double_complex), private, allocatable ::  FNz(:), FNzm1(:), gam(:), WR(:)

   private :: Current, rhs, pendulumODE, u, calc_theta0, ltridag, rtridag

contains
   subroutine gyr(Nz, Nt, Ne, Delta, Ic, dz, dt, ZAxis, TAxis, &
                  InitialField, tol, INTT, INTZ, OUTNz, OUTNt, OUTF, OUTCu, OUTZAxis, OUTTAxis) bind(c, name='gyr')
      implicit none

      integer(c_int), intent(in) :: INTT, INTZ, OUTNz, OUTNt, Nz, Nt, Ne
      real(c_double), intent(in) :: Delta, Ic, dz, dt, tol
      real(c_double), intent(in) :: ZAxis(0:Nz), TAxis(0:Nt)
      complex(c_double_complex), intent(in) :: InitialField(0:Nz)
      real(c_double), intent(inout) :: OUTZAxis(0:OUTNz), OUTTAxis(0:OUTNt)
      complex(c_double_complex), intent(out) :: OUTF(0:OUTNz, 0:OUTNt), OUTCu(0:OUTNz, 0:OUTNt)

      real(c_double) maxdiff, maxfield, maxfield_p

      complex(c_double_complex), allocatable :: field(:), lfield(:), rfield(:), field_p(:), lfield_p(:), rfield_p(:), &
                                                A(:), B(:), C(:), D(:), CuNz(:), CuNzm1(:), &
                                                cu(:), cu_p(:), tmp(:), SigmaNz(:), SigmaNzm1(:), kpar2(:)
      real(c_double), allocatable :: theta(:, :), dthdz(:, :)
      integer(c_int), allocatable :: IZ(:), IT(:)

      allocate (field(0:Nz), lfield(0:Nz), rfield(0:Nz), field_p(0:Nz), lfield_p(0:Nz), rfield_p(0:Nz), &
                A(0:Nz), B(0:Nz - 1), C(1:Nz), D(0:Nz), gam(Nz + 1), WR(0:Nt), &
                FNz(0:Nt), FNzm1(0:Nt), theta(0:Nz, Ne), dthdz(0:Nz, Ne), &
                cu(0:Nz), cu_p(0:Nz), IZ(0:OUTNz), IT(0:OUTNt), CuNz(0:Nt), CuNzm1(0:Nt), &
                tmp(0:Nz), SigmaNz(0:Nt), SigmaNzm1(0:Nt), kpar2(0:Nz), stat=err_alloc)
      if (err_alloc /= 0) then
         print *, "allocation error"
         pause
         stop
      end if

      C0 = dcmplx(0, -1.0D0)
      C2 = 1.0D0/cdsqrt(Im1*pi)
      DeltaT = dt
      Nzm1 = Nz - 1
      SQRDT = dsqrt(dt)
      SQRDZ = dz*dz
#if __INTEL_COMPILER
      IZ(:) = (/0:Nz:INTZ/)
      IT(:) = (/0:Nt:INTT/)
#else
      do i = 0, OUTNz
         IZ(i) = i*INTZ
         OUTZAxis(i) = ZAxis(i*INTZ)
      end do
      do i = 0, OUTNt
         IT(i) = i*INTT
         OUTTAxis(i) = TAxis(i*INTT)
      end do
#endif
      do i = 0, OUTNz
         OUTZAxis(i) = ZAxis(i*INTZ)
      end do
      do i = 0, OUTNt
         OUTTAxis(i) = TAxis(i*INTT)
      end do

      WNz = -((2.0D0/3.0D0*C0*dz/dt + kpar2(Nz)*dz/3.0D0) - 1.0D0/dz)
      WNzm1 = -((C0/3.0D0*dz/dt + kpar2(Nzm1)*dz/6.0D0) + 1.0D0/dz); 
      A(0) = 1.0D0
      A(1:Nzm1) = -2.0D0*(1.0D0 - dz*dz/dt*C0 - dz*kpar2(1:Nzm1)*dz/2.0D0)
      A(Nz) = 1.0D0 + 4.0D0/3.0D0*C2*WNz*SQRDT
      B(0) = 0
      B(1:Nzm1) = 1.0D0
      C(1:Nzm1) = 1.0D0
      C(Nz) = 4.0D0/3.0D0*C2*WNzm1*SQRDT

      !SQRDTm4d3 = 4.0d0/3.0d0*SQRDT
      !SQRDTm2d3 = 2.0d0/3.0d0*SQRDT
      !C2mSQRDTm4d3 = C2*SQRDT*4.0d0/3.0d0
      !C0mSQRDZdDTm2 = 2.0d0*C0*SQRDZ/dt !C0mSQRDZdDTm2 = 2.0d0*C0*dz/dt*dz
      !C0m2d3dDT = C0*2.0d0/3.0d0/dt
      !DZmC0d3dDT = dz/dt*C0/3.0d0
      !DZmC0m2d3dDT = dz/dt*C0*2.0d0/3.0d0
      !DeltaZd3 = dz/3.0d0
      !DeltaZd6 = dz/6.0d0

      !Initial values
      kpar2(:) = 0
      field(:) = InitialField(:)
      OUTF(:, 0) = InitialField(:)
      SigmaNz(0) = 0
      SigmaNzm1(0) = 0
      call pendulumODE(theta, dthdz, field, Ne, Nz, dz)
      call calc_theta0(theta, dthdz, Ne, Delta)
      call Current(cu, theta, Ne, Nz, Ic)
      FNz(0) = field(Nz)
      FNzm1(0) = field(Nzm1)
      CuNz(0) = cu(Nz)
      CuNzm1(0) = cu(Nzm1)
      OUTCu(:, 0) = cu(IZ)

      WR(0) = dz*(1.0D0/6.0D0*(2.0D0*CuNz(0) + CuNzm1(0))); 
      write (*, '(/)')

      time_loop: do step = 1, Nt
         !vezde nije step == j

         WR(step) = dz*((C0*2.0D0/3.0D0/dt - kpar2(Nz)/3.0D0)*FNz(step - 1) &
                        + (C0/3.0D0/dt - kpar2(Nzm1)/6.0D0)*FNzm1(step - 1) &
                        + 1.0D0/6.0D0*(4.0D0*CuNz(step - 1) + 2.0D0*CuNzm1(step - 1)) &
                        - (2.0D0*SigmaNz(step - 1) + SigmaNzm1(step - 1)))
         if (step == 1) then
            IR = 0
         elseif (step == 2) then
            IR = 4.0D0/3.0D0*SQRDT*(u(0)*(1 - SQR2D2) + u(1)*(SQR2M2 - 2.5D0))
         else
            IR = 4.0D0/3.0D0*SQRDT*(u(0)*((step - 1.0D0)**(1.5) - (step - 1.5D0)*dsqrt(dble(step))) + u(step - 1)*(SQR2M2 - 2.5D0))
            do j = 1, step - 2
               IR = IR + 4.0D0/3.0D0*SQRDT*(u(j)*((step - j - 1.0D0)**(1.5) - 2.0D0*(step - j)**(1.5) + (step - j + 1.0D0)**(1.5)))
            end do
         end if

         D(1) = 0
         D(1:Nzm1) = SQRDZ*(2.0D0*cu(1:Nzm1)) &
                     + 2.0D0*(1.0D0 + C0*SQRDZ/dt - SQRDZ*kpar2(1:Nzm1)/2.0D0)*field(1:Nzm1) &
                     - (field(0:Nz - 2) + field(2:Nz))
         D(Nz) = -C2*(IR + 4.0D0/3.0D0*WR(step)*SQRDT + 2.0D0/3.0D0*SQRDT*(WNzm1*field(Nzm1) &
                                                                           + WNz*field(Nz) + WR(step - 1))*cdexp(CR*dt))

         call ltridag(C, A, B, D, lfield_p)
         call rtridag(C, A, B, D, rfield_p)
         field_p(:) = (lfield_p + rfield_p)/2.0D0 ! nesamosoglas. pole

         maxfield = maxval(cdabs(field))
         inloopnum = 0
         do
            !inloopnum = inloopnum + 1

            call pendulumODE(theta, dthdz, field_p, Ne, Nz, dz)
            call Current(cu_p, theta, Ne, Nz, Ic)

            WR(step) = dz*((C0*2.0D0/3.0D0/dt - kpar2(Nz)/3.0D0)*field(Nz) &
                           + (C0/3.0D0/dt - kpar2(Nzm1)/6.0D0)*field(Nzm1) &
                           + 1.0D0/6.0D0*(2.0D0*cu_p(Nz) + 2.0D0*cu(Nz) + cu_p(Nzm1) + cu(Nzm1)) &
                           - (2.0D0*SigmaNz(step - 1) + SigmaNzm1(step - 1)))

            D(1:Nzm1) = SQRDZ*(cu_p(1:Nzm1) + cu(1:Nzm1)) &
                        + 2.0D0*(1 + C0*SQRDZ/dt - SQRDZ*kpar2(1:Nzm1)/2.0D0) *field(1:Nzm1) &
                        - (field(0:Nz - 2) + field(2:Nz))
            D(Nz) = -C2*(IR + 4.0D0/3.0D0*WR(step)*SQRDT + 2.0D0/3.0D0*SQRDT*(WNzm1*field(Nzm1) &
                                                                              + WNz*field(Nz) + WR(step - 1))*cdexp(CR*dt))

            call ltridag(C, A, B, D, lfield_p)
            call rtridag(C, A, B, D, rfield_p)
            field_p(:) = (lfield_p + rfield_p)/2.0D0 ! samosoglas. pole

            maxfield_p = maxval(cdabs(field_p))
            maxdiff = dabs(maxfield - maxfield_p)/maxfield
            !print *, 'innerloop # ', inloopnum, 'maxdiff =', maxdiff
            if (maxdiff .le. tol) exit
            maxfield = maxfield_p
         end do

         field(:) = field_p(:)
         call pendulumODE(theta, dthdz, field, Ne, Nz, dz)
         call Current(cu, theta, Ne, Nz, Ic)

         FNz(step) =  field(Nz)
         FNzm1(step) = field(Nzm1)
         CuNz(step) = cu(Nz)
         CuNzm1(step) = cu(Nzm1)
    
         SigmaNz(step) = -(kpar2(Nz)/6.0D0 + C0/3.0D0/dt) * FNz(step) &
            + (C0/3.0D0/dt - kpar2(Nz)/6.0D0) * FNz(step-1) &
            + 1.0D0/6.0D0*(CuNz(step) + CuNz(step-1)) - SigmaNz(step-1)
         SigmaNzm1(step) = -(kpar2(Nzm1)/6.0D0 + C0/3.0D0/dt) * FNzm1(step) &
            + (C0/3.0D0/dt - kpar2(Nzm1)/6.0D0) * FNzm1(step-1) &
            + 1.0D0/6.0D0*(CuNzm1(step) + CuNzm1(step-1)) - SigmaNzm1(step-1)

         if (mod(step, INTT) .eq. 0) then
            OUTF(:, jout) = field(IZ)
            OUTCu(:, jout) = cu(IZ)
            jout = jout + 1
         end if

#if __INTEL_COMPILER
         write (*, '(a,i10,a,f10.5,a,1pe17.8,a,1pe17.8,a,1pe17.8,\,a)') 'Step = ', step, '   Time = ', dt*step, &
            '   MAX|B| = ', maxval(cdabs(field)), '   MAX|J| = ', maxval(cdabs(cu_p)), '   MAX DIFF = ', maxdiff, char(13)
#else
         write (*, '(a,i10,a,f10.5,a,1pe17.8,a,1pe17.8,a,1pe17.8,a)', advance="no") 'Step = ', step, '   Time = ', dt*step, &
            '   MAX|B| = ', maxval(cdabs(field)), '   MAX|J| = ', maxval(cdabs(cu_p)), '   MAX DIFF = ', maxdiff, char(13)
#endif

      end do time_loop

      OUTF(:, OUTNt) = field(IZ) ! sohranyaem posledniy result

      write (*, '(/)')

      return
101   print *, "error of file open"; pause; stop
102   print *, 'error of file read'; pause; stop
103   print *, 'error of file write'; pause; stop
   end subroutine

   subroutine Current(cu, theta, Ne, Nz, Ic)

      implicit none

      real(c_double), intent(in) :: theta(0:, :), Ic
      complex(c_double_complex), intent(inout) :: cu(0:)
      integer(c_int) Ne, Nz, i

      do i = 0, Nz
         cu(i) = Ic*2.0d0/Ne*sum(cdexp(-Im1*theta(i, :)))
      end do

   end subroutine Current

   subroutine rtridag(c, a, b, d, u)
      use util

      implicit none
      complex(c_double_complex), dimension(:), intent(in) :: c, a, b, d
      complex(c_double_complex), dimension(:), intent(out) :: u
      !complex(c_double_complex), dimension(size(a)) :: gam
      integer(c_int) :: n, j
      complex(c_double_complex) :: bet

      n = assert_eq((/size(c) + 1, size(a), size(b) + 1, size(d), size(u)/), 'tridag_ser')
      bet = a(1)
      if (bet == 0.0) call error('tridag_ser: error at code stage 1')
      u(1) = d(1)/bet
      do j = 2, n
         gam(j) = b(j - 1)/bet
         bet = a(j) - c(j - 1)*gam(j)
         if (bet == 0.0) &
            call error('tridag_ser: error at code stage 2')
         u(j) = (d(j) - c(j - 1)*u(j - 1))/bet
      end do
      do j = n - 1, 1, -1
         u(j) = u(j) - gam(j + 1)*u(j + 1)
      end do
   end subroutine rtridag

   subroutine ltridag(c, a, b, d, u)
      use util

      implicit none
      complex(c_double_complex), dimension(:), intent(in) :: c, a, b, d
      complex(c_double_complex), dimension(:), intent(out) :: u
      !complex(c_double_complex), dimension(size(a)) :: gam
      integer(c_int) :: n, j
      complex(c_double_complex) :: bet

      n = assert_eq((/size(c) + 1, size(a), size(b) + 1, size(d), size(u)/), 'tridag_ser')
      bet = a(n)
      if (bet == 0.0) call error('tridag_ser: error at code stage 1')
      u(n) = d(n)/bet
      do j = n - 1, 1, -1
         gam(j) = c(j)/bet
         bet = a(j) - b(j)*gam(j)
         if (bet == 0.0) &
            call error('tridag_ser: error at code stage 2')
         u(j) = (d(j) - b(j)*u(j + 1))/bet
      end do
      do j = 2, n
         u(j) = u(j) - gam(j - 1)*u(j - 1)
      end do
   end subroutine ltridag

   subroutine pendulumODE(theta, dthdz, F, Ne, Nz, h)
      implicit none

      real(c_double), intent(inout) :: theta(0:, :), dthdz(0:, :)
      complex(c_double_complex), intent(in) :: F(0:)
      integer(c_int) i, Ne, Nz
      real(c_double) rhs0(Ne), rhs1(Ne), h

      do i = 0, Nz - 1
         rhs0 = rhs(F(i), theta(i, :)); 
         theta(i + 1, :) = theta(i, :) + dthdz(i, :)*h + h/2.0D0*rhs0*h; 
         rhs1 = rhs(F(i + 1), theta(i + 1, :)); 
         theta(i + 1, :) = theta(i, :) + dthdz(i, :)*h + h/6.0D0*rhs0*h + h/3.0D0*rhs1*h; 
         dthdz(i + 1, :) = dthdz(i, :) + h/2.0D0*(rhs0 + rhs1); 
      end do
   end subroutine pendulumODE

   function rhs(F, th)
      implicit none

      real(c_double) th(:), rhs(size(th, 1))
      complex(c_double_complex) F

      rhs(:) = dimag(F*cdexp(Im1*th))
   end function rhs

   subroutine calc_theta0(th, dthdz, Ne, delta)
      implicit none

      real(c_double), intent(inout) :: th(:, :), dthdz(:, :)
      real(c_double), intent(in) ::  delta
      real(c_double) h
      integer(c_int), dimension(size(th, 2)) :: i
      integer(c_int) Ne, j

      h = 2.0d0*pi/Ne

#if __INTEL_COMPILER
      i = (/0:Ne - 1/)
#else
      do j = 1, size(th, 2)
         i(j) = j - 1
      end do
#endif

      th(1, :) = h*i
      dthdz(1, :) = delta
   end subroutine calc_theta0

   function u(j)
      implicit none

      integer(c_int) j
      complex(c_double_complex) u

      u = (WNzm1*FNzm1(j) + WNz*FNz(j) + WR(j))*cdexp(CR*DeltaT*(step - j))

   end function u

end module functions
