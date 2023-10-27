module functions
   use, intrinsic :: iso_c_binding

   complex(c_double_complex), private :: C0 = (0, -1), CR = 0, C2
   real(c_double), private :: SQR2 = dsqrt(2.0D0), SQR2M2 = 2.828427124746190, SQR2D2 = 0.707106781186548

   complex(c_double_complex), private :: IR, IROldPart, WNz, WNzm1, &
                                         WR_step_part, D_end_part, C0mSQRDZdDTm2, C0m2d3dDT, &
                                         DZmC0m2d3dDT, DZmC0d3dDT, C2mSQRDTm4d3
   real(c_double), private :: DeltaT, SQRDT, SQRDZ, SQRDTm4d3, SQRDTm2d3, &
                              DeltaZd3, DeltaZd6
   integer(c_int), private :: i, j, err_alloc = 1, step, jout = 1, Nzm1, inloopnum

   complex(c_double_complex), private :: Im1 = (0.0d0, 1.0d0)
   real(c_double), private :: pi = dacos(-1.0d0) !real(c_double), parameter :: pi = 2.0d0*dacos(0.0d0)

   complex(c_double_complex), private, allocatable :: Field_sam(:), lField_sam(:), rField_sam(:), Field_nesam(:), lField_nesam(:), rField_nesam(:), &
                                                      A(:), B(:), C(:), D(:), D_1_Nzm1_part(:), gam(:), WR(:), &
                                                      FNz(:), FNzm1(:), OldCu(:), Cup(:), tmp(:), SigmaNz(:), SigmaNzm1(:)
   real(c_double), private, allocatable :: theta(:, :), dthdz(:, :)!, kpar2(:)
   integer(c_int), private, allocatable :: IZ(:), IT(:)

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

      real(c_double) maxdiff

      allocate (Field_sam(0:Nz), lField_sam(0:Nz), rField_sam(0:Nz), Field_nesam(0:Nz), lField_nesam(0:Nz), rField_nesam(0:Nz), &
                A(0:Nz), B(0:Nz - 1), C(1:Nz), D(0:Nz), D_1_Nzm1_part(1:Nz - 1), gam(Nz + 1), WR(0:Nt), &
                FNz(0:Nt), FNzm1(0:Nt), theta(0:Nz, Ne), dthdz(0:Nz, Ne), &
                OldCu(0:Nz), Cup(0:Nz), IZ(0:OUTNz), IT(0:OUTNt), &
                tmp(0:Nz), SigmaNz(0:Nt), SigmaNzm1(0:Nt), stat=err_alloc)
      if (err_alloc /= 0) then
         print *, "allocation error"
         pause
         stop
      end if

      C2 = 1.0d0/cdsqrt(Im1*pi)
      DeltaT = dt
      Nzm1 = Nz - 1
      SQRDT = dsqrt(dt)
      SQRDZ = dz*dz
      Field_sam(:) = InitialField(:)
      OUTF(:, 0) = InitialField(:)
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

      WNz = -(2.0d0/3.0d0*C0*dz/dt - 1.0d0/dz)
      WNzm1 = -(C0/3.0d0*dz/dt + 1.0d0/dz)
      A(0) = 1
      A(1:Nzm1) = -2.0d0*(1.0d0 - dz/dt*C0*dz)
      A(Nz) = 1.0d0 + 4.0d0/3.0d0*C2*WNz*SQRDT
      B(0) = 0
      B(1:Nzm1) = 1
      C(1:Nzm1) = 1
      C(Nz) = 4.0d0/3.0d0*C2*WNzm1*SQRDT

      FNz(:) = 0
      FNzm1(:) = 0
      FNz(0) = Field_sam(Nz)
      FNzm1(0) = Field_sam(Nzm1)
      WR(:) = 0

      SQRDTm4d3 = 4.0d0/3.0d0*SQRDT
      SQRDTm2d3 = 2.0d0/3.0d0*SQRDT
      C2mSQRDTm4d3 = C2*SQRDT*4.0d0/3.0d0
      C0mSQRDZdDTm2 = 2.0d0*C0*SQRDZ/dt !C0mSQRDZdDTm2 = 2.0d0*C0*dz/dt*dz
      C0m2d3dDT = C0*2.0d0/3.0d0/dt
      DZmC0d3dDT = dz/dt*C0/3.0d0
      DZmC0m2d3dDT = dz/dt*C0*2.0d0/3.0d0
      DeltaZd3 = dz/3.0d0
      DeltaZd6 = dz/6.0d0

      call calc_theta0(theta, dthdz, Ne, Delta)

      write (*, '(/)')

      time_loop: do step = 1, Nt
         !vezde nije step == j
         !Old -> (step - 1); OldOld -> (step - 2)

         call pendulumODE(theta, dthdz, Field_sam, Ne, Nz, dz)
         call Current(OldCu, theta, Ne, Nz, Ic)

         if ((step /= 1) .and. (mod(step - 1, INTT) == 0)) then
            OUTCu(:, jout) = OldCu(IZ)
         end if

         WR_step_part = DZmC0m2d3dDT*FNz(step - 1) &
                        + DZmC0d3dDT*FNzm1(step - 1) &
                        !- (DeltaZm2*OldSigmaNz + dz*OldSigmaNzm1)
                        - (2.0d0*SigmaNz(step - 1) + SigmaNzm1(step - 1))
         !WR_step_part = dz*((C0m2d3dDT - kpar2(Nz)/3.0d0)*Field_sam(Nz) &
         !                       + (C0d3dDeltaT - kpar2(Nzm1)/6.0d0)*Field_sam(Nzm1) &
         !                       - (2.0d0*OldSigmaNz + OldSigmaNzm1))
         WR(step) = WR_step_part + 2.0d0*DeltaZd3*OldCu(Nz) + 2.0d0*DeltaZd6*OldCu(Nzm1)

         if (step == 1) then
            IR = 0.0d0
         elseif (step == 2) then
            IR = SQRDTm4d3*(u(0)*(1.0d0 - SQR2D2) + u(1)*(SQR2M2 - 2.5d0))
         else
            IR = u(0)*((step - 1.0d0)**(1.5) - (step - 1.5d0)*dsqrt(dble(step))) + u(step - 1)*(SQR2M2 - 2.5d0)
            do j = 1, step - 2
               IR = SQRDTm4d3*(IR + u(j)*((step - j - 1)**(1.5) - 2*(step - j)**(1.5) + (step - j + 1)**(1.5)))
            end do
         end if

         D(0) = 0
         D_1_Nzm1_part = (2.0d0 + C0mSQRDZdDTm2)*Field_sam(1:Nzm1) - (Field_sam(0:Nz - 2) + Field_sam(2:Nz))
         D(1:Nzm1) = D_1_Nzm1_part + SQRDZ*2.0d0*OldCu(1:Nzm1)
         D_end_part = -C2*(IR + SQRDTm2d3*u(step - 1))
         D(Nz) = D_end_part - C2mSQRDTm4d3*WR(step)

         call ltridag(C, A, B, D, lField_nesam)
         call rtridag(C, A, B, D, rField_nesam)
         Field_nesam(:) = (lField_nesam + rField_nesam)/2.0d0 ! nesamosoglas. pole

         inloopnum = 0
         do
            !inloopnum = inloopnum + 1

            call pendulumODE(theta, dthdz, Field_nesam, Ne, Nz, dz)
            call Current(Cup, theta, Ne, Nz, Ic)

            WR(step) = WR_step_part + DeltaZd3*(Cup(Nz) + OldCu(Nz)) + DeltaZd6*(Cup(Nzm1) + OldCu(Nzm1)); 
            D(1:Nzm1) = D_1_Nzm1_part + SQRDZ*(Cup(1:Nzm1) + OldCu(1:Nzm1))
            D(Nz) = D_end_part - C2mSQRDTm4d3*WR(step)

            call ltridag(C, A, B, D, lField_sam)
            call rtridag(C, A, B, D, rField_sam)
            Field_sam(:) = (lField_sam + rField_sam)/2.0d0 ! samosoglas. pole

            maxdiff = maxval(cdabs(Field_sam - Field_nesam))/maxval(cdabs(Field_sam))
            !print *, 'innerloop # ', inloopnum, 'maxdiff =', maxdiff

            if (maxdiff .le. tol) exit

            Field_nesam = Field_sam
         end do

         FNz(step) = Field_sam(Nz)
         FNzm1(step) = Field_sam(Nz - 1)

         SigmaNz(step) = -DZmC0d3dDT*FNz(step) &
                         + DZmC0d3dDT*FNz(step - 1) &
                         + DeltaZd6*(Cup(Nz) + OldCu(Nz)) - dz*SigmaNz(step - 1)
         SigmaNzm1(step) = -DZmC0d3dDT*FNzm1(step) &
                           + DZmC0d3dDT*FNzm1(step - 1) &
                           + DeltaZd6*(Cup(Nzm1) + OldCu(Nzm1)) - dz*SigmaNzm1(step - 1)

         if (mod(step, INTT) .eq. 0) then
            OUTF(:, jout) = Field_sam(IZ)
            jout = jout + 1
         end if

#if __INTEL_COMPILER
         write (*, '(a,i10,a,f10.5,a,1pe17.8,a,1pe17.8,a,1pe17.8,\,a)') 'Step = ', step, '   Time = ', dt*step, &
            '   MAX|B| = ', maxval(cdabs(Field_sam)), '   MAX|J| = ', maxval(cdabs(Cup)), '   MAX DIFF = ', maxdiff, char(13)
#else
         write (*, '(a,i10,a,f10.5,a,1pe17.8,a,1pe17.8,a,1pe17.8,a)', advance="no") 'Step = ', step, '   Time = ', dt*step, &
            '   MAX|B| = ', maxval(cdabs(Field_sam)), '   MAX|J| = ', maxval(cdabs(Cup)), '   MAX DIFF = ', maxdiff, char(13)
#endif

      end do time_loop

      OUTF(:, OUTNt) = Field_sam(IZ) ! sohranyaem posledniy result

      write (*, '(/)')

      return
101   print *, "error of file open"; pause; stop
102   print *, 'error of file read'; pause; stop
103   print *, 'error of file write'; pause; stop
   end subroutine

   subroutine Current(Cu, theta, Ne, Nz, Ic)

      implicit none

      real(c_double), intent(in) :: theta(0:, :), Ic
      complex(c_double_complex), intent(inout) :: Cu(0:)
      integer(c_int) Ne, Nz, i

      do i = 0, Nz
         Cu(i) = Ic*2.0d0/Ne*sum(cdexp(-Im1*theta(i, :)))
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
         theta(i + 1, :) = theta(i, :) + dthdz(i, :)*h + h/2.0*rhs0*h; 
         rhs1 = rhs(F(i + 1), theta(i + 1, :)); 
         theta(i + 1, :) = theta(i, :) + dthdz(i, :)*h + h/6.0*rhs0*h + h/3.0*rhs1*h; 
         dthdz(i + 1, :) = dthdz(i, :) + h/2.0*(rhs0 + rhs1); 
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
