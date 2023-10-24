module functions
   use, intrinsic :: iso_c_binding
   use types

   complex(c_double_complex) :: C0 = (0, 1), CR = 0, C2
   real(c_double) :: SQR2 = dsqrt(2.0D0), SQR2M2 = 2.828427124746190, &
                     SQR2D2 = 0.707106781186548, ONEminusSQR2D2, SQR2M2minus2p5

   complex(c_double_complex) :: IR, IROldPart, WNz, WNzm1, OldOldSigmaNz, OldOldSigmaNzm1, OldSigmaNz, OldSigmaNzm1, &
                                WR_step_part, D_end_part, C0mSQRDZdDeltaTm2, C0m2d3dDeltaT, &
                                DeltaZmC0m2d3dDeltaT, DeltaZmC0d3dDeltaT
   real(c_double) :: DeltaZ, DeltaT, SQRDT, SQRDZ, SQRDTm4d3, SQRDTm2d3, C2mSQRDTm4d3, &
                     DeltaZd2, DeltaZd3, DeltaZd6
   integer(c_int) :: i, j, err_alloc = 1, step, jout = 1, Nzm1

   integer, pointer :: Ne, OUTNz, Nz, Nt
   real(c_double), pointer :: ZAxis(:), TAxis(:), Delta, Ic
   integer(c_int), pointer :: INTERVALT, INTERVALZ
   complex(c_double_complex), pointer :: InitialField(:), OUTB(:, :), OUTCu(:, :)

   complex(c_double_complex), allocatable :: Field(:), lField(:), rField(:), Field_p(:), lField_p(:), rField_p(:), &
                                             A(:), B(:), C(:), D(:), D_2_Nzm1_part(:), gam(:), WR(:), &
                                             FNz(:), FNzm1(:), Cu(:), Cup(:), OldOldCu(:), tmp(:)!, CuNz(:), CuNzm1(:)
   real(c_double), allocatable :: TAxisNew(:), theta(:, :), dthdz(:, :), kpar2(:)
   integer(c_int), allocatable :: IZ(:)

contains
   subroutine oro(INP, OUTP) bind(c, name='oro')
      implicit none

      type(Input), target, intent(in) :: INP
      type(OUtput), target, intent(inout) :: OUTP

      C2 = 1.0d0/cdsqrt(-Im1*pi)
      SQR2M2minus2p5 = SQR2M2 - 2.5d0
      ONEminusSQR2D2 = 1.0d0 - SQR2D2
      Ne => INP%Ne
      ZAxis(0:) => INP%ZAxis(:)
      TAxis(0:) => INP%TAxis(:)
      Delta => INP%Delta
      Ic => INP%Ic
      INTERVALT => INP%INTERVALT
      INTERVALZ => INP%INTERVALZ
      InitialField(0:) => INP%InitialField(:)
      OUTB(0:, 0:) => OUTP%OUTB(:, :)
      OUTCu(0:, 0:) => OUTP%OUTCu(:, :)
      OUTNz => INP%OUTNz
      Nz => INP%Nz
      Nt => INP%Nt
      DeltaZ = INP%dz
      DeltaT = INP%dt

      allocate (Field(0:Nz), lField(0:Nz), rField(0:Nz), Field_p(0:Nz), lField_p(0:Nz), rField_p(0:Nz), &
                A(0:Nz), B(0:Nz - 1), C(1:Nz), D(0:Nz), D_2_Nzm1_part(1:Nz - 1), gam(Nz + 1), WR(-1:Nt), &
                FNz(-1:Nt), FNzm1(-1:Nt), OldOldCu(0:Nz), theta(0:Nz, Ne), dthdz(0:Nz, Ne), &
                Cu(0:Nz), Cup(0:Nz), IZ(0:OUTNz), kpar2(0:Nz), tmp(0:Nz), stat=err_alloc)
      if (err_alloc /= 0) then
         print *, "allocation error"
         pause
         stop
      end if

      Nzm1 = Nz - 1
      SQRDT = dsqrt(DeltaT)
      SQRDZ = DeltaZ*DeltaZ
      Field(:) = InitialField(:)
      IZ(:) = (/0:Nz:INTERVALZ/)
      OldOldSigmaNz = dcmplx(0)
      OldOldSigmaNzm1 = dcmplx(0)

      kpar2(:) = 0
      WNz = -(2.0d0/3.0d0*C0*DeltaZ/DeltaT + kpar2(Nz)*DeltaZ/3.0d0 - 1.0d0/DeltaZ)
      WNzm1 = -(C0/3.0d0*DeltaZ/DeltaT + kpar2(Nzm1)*DeltaZ/6.0d0 + 1.0d0/DeltaZ)
      A(0) = 1
      A(1:Nzm1) = -2.0d0*(1.0d0 - DeltaZ/DeltaT*C0*DeltaZ - DeltaZ*kpar2(1:Nzm1)*DeltaZ/2.0d0)
      A(Nz) = 1.0d0 + 4.0d0/3.0d0*C2*WNz*SQRDT
      B(0) = 0
      B(1:Nzm1) = 1
      C(1:Nzm1) = 1
      C(Nz) = 4.0d0/3.0d0*C2*WNzm1*SQRDT

      FNz(:) = 0
      FNzm1(:) = 0
      WR(:) = 0
      OldOldCu(:) = 0

      SQRDTm4d3 = 4.0d0/3.0d0*SQRDT
      SQRDTm2d3 = 2.0d0/3.0d0*SQRDT
      C2mSQRDTm4d3 = C2*SQRDT*4.0d0/3.0d0
      C0mSQRDZdDeltaTm2 = 2.0d0*C0*DeltaZ/DeltaT*DeltaZ
      C0m2d3dDeltaT = C0*2.0d0/3.0d0/DeltaT
      DeltaZmC0d3dDeltaT = DeltaZ/DeltaT*C0/3.0d0
      DeltaZmC0m2d3dDeltaT = DeltaZ/DeltaT*C0*2.0d0/3.0d0
      DeltaZd2 = DeltaZ/2.0d0
      DeltaZd3 = DeltaZ/3.0d0
      DeltaZd6 = DeltaZ/6.0d0

      call calc_theta0(theta, dthdz, Ne, Delta)

      write (*, '(/)')

      time_loop: do step = 1, Nt
         !vezde nije step == j
         !Old -> (step - 1); OldOld -> (step - 2)

         call pendulumODE(theta, dthdz, Field, Ne, Nz, DeltaZ)
         call Current(Cu, theta, Nz, Ic)

         if ((step /= 1) .and. (mod(step - 1, INTERVALT) == 0)) then
            OUTCu(:, jout) = Cu(IZ)
         end if

         FNz(step - 1) = Field(Nz)
         FNzm1(step - 1) = Field(Nzm1)

         !OldSigmaNz = -(kpar2(Nz)/6.0d0 + C0d3dDeltaT)*FNz(step - 1) &
         !             + (C0d3dDeltaT - kpar2(Nz)/6.0d0)*FNz(step - 2) &
         !             + 0.166666666666667*(Cu(Nz) + OldOldCu(Nz)) - OldOldSigmaNz
         !OldSigmaNzm1 = -(kpar2(Nzm1)/6.0d0 + C0d3dDeltaT)*FNzm1(step - 1) &
         !               + (C0d3dDeltaT - kpar2(Nzm1)/6.0d0)*FNzm1(step - 2) &
         !               + 0.166666666666667*(Cu(Nzm1) + OldOldCu(Nzm1)) - OldOldSigmaNzm1
         OldSigmaNz = -(DeltaZd6*kpar2(Nz) + DeltaZmC0d3dDeltaT)*FNz(step - 1) &
                      + (DeltaZmC0d3dDeltaT - DeltaZd6*kpar2(Nz))*FNz(step - 2) &
                      + DeltaZd6*(Cu(Nz) + OldOldCu(Nz)) - DeltaZ*OldOldSigmaNz
         OldSigmaNzm1 = -(DeltaZd6*kpar2(Nzm1) + DeltaZmC0d3dDeltaT)*FNzm1(step - 1) &
                        + (DeltaZmC0d3dDeltaT - DeltaZd6*kpar2(Nzm1))*FNzm1(step - 2) &
                        + DeltaZd6*(Cu(Nzm1) + OldOldCu(Nzm1)) - DeltaZ*OldOldSigmaNzm1

         OldOldSigmaNz = OldSigmaNz
         OldOldSigmaNzm1 = OldSigmaNzm1

         WR_step_part = (DeltaZmC0m2d3dDeltaT - DeltaZd3*kpar2(Nz))*Field(Nz) &
                        + (DeltaZmC0d3dDeltaT - DeltaZd6*kpar2(Nzm1))*Field(Nzm1) &
                        !- (DeltaZm2*OldSigmaNz + DeltaZ*OldSigmaNzm1)
                        - (2.0d0*OldSigmaNz + OldSigmaNzm1)
         !WR_step_part = DeltaZ*((C0m2d3dDeltaT - kpar2(Nz)/3.0d0)*Field(Nz) &
         !                       + (C0d3dDeltaT - kpar2(Nzm1)/6.0d0)*Field(Nzm1) &
         !                       - (2.0d0*OldSigmaNz + OldSigmaNzm1))
         WR(step) = WR_step_part + DeltaZ*Cu(Nz) - DeltaZd3*OldOldCu(Nz) + DeltaZd2*Cu(Nzm1) - DeltaZd6*OldOldCu(Nzm1)
         !WR(step) = WR_step_part + DeltaZ*(0.166666666666667*(6.0d0*Cu(Nz) - 2.0d0*OldOldCu(Nz) + 3.0d0*Cu(Nzm1) - OldOldCu(Nzm1)))

         if (step == 1) then
            IR = 0.0d0
         elseif (step == 2) then
            IR = SQRDTm4d3*(u(0)*(ONEminusSQR2D2) + u(1)*(SQR2M2minus2p5))
         else
            IR = u(0)*((step - 1.0d0)**(1.5) - (step - 1.5d0)*dsqrt(dble(step))) + u(step - 1)*(SQR2M2minus2p5)
            do j = 1, step - 2
               IR = SQRDTm4d3*(IR + u(j)*((step - j - 1)**(1.5) - 2*(step - j)**(1.5) + (step - j + 1)**(1.5)))
            end do
         end if

         D(0) = 0
         D_2_Nzm1_part = (2.0d0 + C0mSQRDZdDeltaTm2 - SQRDZ*kpar2(1:Nzm1))*Field(1:Nzm1) - (Field(0:Nz - 2) + Field(2:Nz))
         D(1:Nzm1) = D_2_Nzm1_part + SQRDZ*(2.0d0*Cu(1:Nzm1) - OldOldCu(1:Nzm1))
         D_end_part = -C2*(IR + SQRDTm2d3*u(step - 1))
         D(Nz) = D_end_part - C2mSQRDTm4d3*WR(step)

         call ltridag(C, A, B, D, lField_p)
         call rtridag(C, A, B, D, rField_p)
         Field_p(:) = (lField_p + rField_p)/2.0d0

         call pendulumODE(theta, dthdz, Field_p, Ne, Nz, DeltaZ)
         call Current(Cup, theta, Nz, Ic)

         !WR(step) = WR_step_part + DeltaZ*(0.166666666666667*(2.0d0*(Cup(Nz) + Cu(Nz)) + Cup(Nzm1) + Cu(Nzm1)));
         WR(step) = WR_step_part + DeltaZd3*(Cup(Nz) + Cu(Nz)) + DeltaZd6*(Cup(Nzm1) + Cu(Nzm1)); 
         D(1:Nzm1) = D_2_Nzm1_part + SQRDZ*(Cup(1:Nzm1) + Cu(1:Nzm1))
         D(Nz) = D_end_part - C2mSQRDTm4d3*WR(step)

         call ltridag(C, A, B, D, lField)
         call rtridag(C, A, B, D, rField)
         Field(:) = (lField + rField)/2.0d0

         OldOldCu(:) = Cu(:)

         if (mod(step, INTERVALT) .eq. 0) then
            OUTB(:, jout) = Field(IZ)
            jout = jout + 1
         end if

         write (*, '(a,i10,a,f10.5,a,1pe17.8,a,1pe17.8,\,a)') 'Step = ', step, '   Time = ', DeltaT*step, &
            '   MAX|B| = ', maxval(cdabs(Field)), '   MAX|J| = ', maxval(cdabs(Cu)), char(13)

      end do time_loop

      write (*, '(/)')

      return
101   print *, "error of file open"; pause; stop
102   print *, 'error of file read'; pause; stop
103   print *, 'error of file write'; pause; stop
   end subroutine

   subroutine Current(Cu, theta, Nz, Ic)

      implicit none

      real(c_double), intent(in) :: theta(0:, :), Ic
      complex(c_double_complex), intent(inout) :: Cu(0:)
      integer(c_int) Nz, i

      do i = 0, Nz
         Cu(i) = Ic*2.0d0/Ne*sum(cdexp(-Im1*theta(i, :)))
      end do

   end subroutine Current

   !function Current(theta, Nz, Ic) result(Cu)
   !
   !   implicit none
   !
   !   real(c_double) theta(0:, :), Ic
   !   complex(c_double_complex) Cu(0:size(theta, 1) - 1)
   !   integer(c_int) Nz,i
   !
   !   do i=0,Nz
   !      Cu(i) = Ic*2.0d0/Ne*sum(cdexp(-Im1*theta(i,:)))
   !   enddo
   !
   !end function Current

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
      integer(c_double), dimension(size(th, 2)) :: i
      integer(c_int) Ne

      h = 2.0d0*pi/Ne

      i = (/0:Ne - 1/)

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

!open (1, file='test.dat', err=101)
!do i = 0, Nz
!   write (1, '(f17.8,a,\)') ZAxis(i), '   '
!   do j = 1, Ne
!      write (1, '(F17.8,a,\)', err=103) theta(i, j), '   '
!   end do
!   write(1,'(/,\)')
!end do
!close (1)
!stop

!open (1, file='test.dat')
!do i = 0, Nz
!   write (1, '(1p3e17.8)') ZAxis(i), Field_p(i)
!end do
!close (1)
!stop
