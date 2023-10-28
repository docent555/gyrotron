program gyrotron
   use, intrinsic :: iso_c_binding
   !use types
   use functions
#if __INTEL_COMPILER
   use ifport
#endif

   implicit none

   integer(c_int) hours, minutes, seconds
   real(c_double) start_time, stop_time, calc_time, tol
   !real(c_double), parameter :: pi = 2.0d0*dacos(0.0d0)
   real(c_double) :: pi = dacos(-1.0d0)

   real(c_double) Lz, Tend, Delta, dz, dt, Ic
   real(c_double), allocatable :: ZAxis(:), TAxis(:), ZBEG, ZEND
   integer(c_int) i, j, Ne, Nz, Nt, err_alloc, INTT, INTZ, OUTNz, OUTNt
   real(c_double), allocatable :: OUTZAxis(:), OUTTAxis(:)
   complex(c_double_complex), allocatable :: InitialField(:), OUTF(:, :), OUTCu(:, :)
   logical(c_bool) :: IFR = .false.

   namelist /param/ Ne, Lz, Tend, Delta, Ic, dz, dt, tol

   open (unit=1, file='input_fortran.in', status='old', err=101)
   read (unit=1, nml=param, err=102)
   close (unit=1)

   write (*, nml=param)

   write (*, '(/)')
#if __INTEL_COMPILER
   start_time = dclock()
#endif

   Nz = Lz/dz ! Razmer vektora ZAxis (vvod)
   Nt = Tend/dt ! Razmer vektora TAxis (vvod)
   allocate (ZAxis(0:Nz), TAxis(0:Nt), InitialField(0:Nz), stat=err_alloc)
   if (err_alloc /= 0) then
      print *, "allocation error"
      pause
      stop
   end if

   ! Vychislit' razmery massivov dlya vyvoda
   INTT = Nt/500
   INTZ = Nz/500
   !INTT = 1
   !INTZ = 1

   if (INTT < 1) then
      print *, 'Too small "Tend"'
      pause
      stop
   end if

   if (INTZ < 1) then
      print *, 'Too small "Lz"'
      pause
      stop
   end if

   if (INTT > 1 .and. INTZ > 1) then
      OUTNt = Nt/INTT
      OUTNz = Nz/INTZ
   elseif (INTT == 1 .and. INTZ > 1) then
      OUTNt = Nt
      OUTNz = Nz/INTZ
   elseif (INTT > 1 .and. INTZ == 1) then
      OUTNt = Nt/INTT
      OUTNz = Nz
   else ! (INTT == 1) .and. (INTZ == 1) then
      OUTNt = Nt
      OUTNz = Nz
   end if

   ! Vydelit' pamyat' pod massivy vyvoda
   allocate (OUTF(0:OUTNz, 0:OUTNt), OUTCu(0:OUTNz, 0:OUTNt), OUTZAxis(0:OUTNz), OUTTAxis(0:OUTNt), stat=err_alloc)
   if (err_alloc /= 0) then
      print *, "allocation error"
      pause
      stop
   end if

   do i = 0, Nz
      ZAxis(i) = i*dz
   end do

   do i = 0, Nt
      TAxis(i) = i*dt
   end do

   if (IFR .eqv. .true.) then
      open (1, file='init_field.in', status='old')
      do i = 0, Nz
         read (1, '(1p2e17.7)') InitialField(i)
      end do
      close (1)
   else
      ZBEG = 0
      ZEND = 0.5
      !where ((ZAxis .GT. ZBEG) .AND. (ZAxis .LT. ZEND)) InitialField = 0.001*sin(pi*(ZAxis - ZBEG)/(ZEND - ZBEG))**2
      where ((ZAxis .GT. ZBEG) .AND. (ZAxis .LT. ZEND)) InitialField = sin(pi*(ZAxis - ZBEG)/(ZEND - ZBEG))**2
      !InitialField = dcmplx(10,10)
   end if

   !open (1, file='test.dat')
   !do i = 0, Nz
   !   write (1, '(1p2e17.7)') InitialField(i)
   !end do
   !close (1)
   !stop

   !call oro(INP, OUTP)
   !call oro(Nz, Nt, Ne, INTT, INTZ, OUTNz, Delta, Ic, dz, dt, &
   !         ZAxis, TAxis, InitialField, tol, OUTB, OUTCu)
   call gyr(Nz, Nt, Ne, Delta, Ic, dz, dt, ZAxis, TAxis, &
            InitialField, tol, INTT, INTZ, OUTNz, OUTNt, OUTF, OUTCu, OUTZAxis, OUTTAxis)

   !open(1, file = 'test.dat', err = 101)
   !write(1, '(4I)', err = 103)  size(TAxis,1), INTT, size(ZAxis,1), INTZ
   !do i=1,Nz
   !write(1, '(3F17.8)', err = 103) ZAxis(i), real(InitialField(i)), imag(InitialField(i))
   !enddo
   !close(1)

#if __INTEL_COMPILER
   stop_time = dclock()
   calc_time = stop_time - start_time
   hours = calc_time/3600
   minutes = (calc_time - hours*3600)/60
   seconds = calc_time - hours*3600 - minutes*60

   write (*, '(/)')
   print *, 'Calcualting took:', hours, 'h :', minutes, 'm :', seconds, 's'
#endif

   open (1, file='fvsz.dat')
   do i = 0, OUTNz
      write (1, fmt='(1x,f14.6)', advance="no") OUTZAxis(i)
      do j = 0, OUTNt
         write (1, fmt='(1x,1p2e17.7,a)', advance="no") OUTF(i, j)
      end do
      write (1, *) ! Assumes default "advance='yes'".
   end do
   close (1)

   open (1, file='ivsz.dat')
   do i = 0, OUTNz
      write (1, fmt='(1x,f14.6)', advance="no") OUTZAxis(i)
      do j = 0, OUTNt
         write (1, fmt='(1x,1p2e17.7,a)', advance="no") OUTCu(i, j)
      end do
      write (1, *) ! Assumes default "advance='yes'".
   end do
   close (1)

!   open (1, file='br.dat')
!   do i = 0, OUTNz
!      write (1, fmt='(1x,f14.6)', advance="no") OUTZAxis(i)
!      do j = 0, OUTNt
!!         write (1, '(1pe17.7,a,\)') dreal(OUTP%OUTB(i, j)), ' '
!         write (1, fmt='(1x,1pe17.7,a)', advance="no") dreal(OUTF(i, j))
!      end do
!      write (1, *) ! Assumes default "advance='yes'".
!!      write (1, '(/,\)')
!   end do
!   close (1)
!
!   open (1, file='bi.dat')
!   do i = 0, OUTNz
!      write (1, fmt='(1x,f14.6)', advance="no") OUTZAxis(i)
!      do j = 0, OUTNt
!!         write (1, '(1pe17.7,a,\)') dimag(OUTP%OUTB(i, j)), ' '
!         write (1, fmt='(1x,1pe17.7,a)', advance="no") dimag(OUTF(i, j))
!      end do
!      write (1, *) ! Assumes default "advance='yes'".
!!      write (1, '(/,\)')
!   end do
!   close (1)
!
!   open (1, file='ir.dat')
!   do i = 0, OUTNz
!      write (1, fmt='(1x,f14.6)', advance="no") OUTZAxis(i)
!      do j = 0, OUTNt
!!         write (1, '(1pe17.7,a,\)') dreal(OUTP%OUTCu(i, j)), ' '
!         write (1, fmt='(1x,1pe17.7,a)', advance="no") dreal(OUTCu(i, j))
!      end do
!      write (1, *) ! Assumes default "advance='yes'".
!!      write (1, '(/,\)')
!   end do
!   close (1)
!
!   open (1, file='ii.dat')
!   do i = 0, OUTNz
!      write (1, fmt='(1x,f14.6)', advance="no") OUTZAxis(i)
!      do j = 0, OUTNt
!!         write (1, '(1pe17.7,a,\)') dimag(OUTP%OUTCu(i, j)), ' '
!         write (1, fmt='(1x,1pe17.7,a)', advance="no") dimag(OUTCu(i, j))
!      end do
!      write (1, *) ! Assumes default "advance='yes'".
!!      write (1, '(/,\)')
!   end do
!   close (1)

   pause
   stop
101 print *, "error of file open"; pause; stop
102 print *, 'error of file read'; pause; stop
103 print *, 'error of file write'; pause; stop
end
