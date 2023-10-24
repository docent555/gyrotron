program orotron
   use, intrinsic :: iso_c_binding
   !use types
#if __INTEL_COMPILER
   use ifport
#endif

   implicit none

   !interface
   !   subroutine oro(INP, OUTP) bind(c, name='oro')
   !      use types
   !      type(Input), target, intent(in) :: INP
   !      type(Output), target, intent(inout) :: OUTP
   !   end subroutine
   !end interface
   interface
      subroutine oro(Nz, Nt, Ne, INTERVALT, INTERVALZ, OUTNz, Delta, Ic, dz, dt, &
                     ZAxis, TAxis, InitialField, OUTB, OUTCu) bind(c, name='oro')
         use, intrinsic :: iso_c_binding
         integer(c_int), intent(in) :: INTERVALT, INTERVALZ, OUTNz, Nz, Nt, Ne
         real(c_double), intent(in) :: Delta, Ic, dz, dt
         real(c_double), intent(in) :: ZAxis(:), TAxis(:)
         complex(c_double_complex), intent(in) :: InitialField(:)
         complex(c_double_complex), intent(out) :: OUTB(:, :), OUTCu(:, :)
      end subroutine
   end interface

   integer(c_int) hours, minutes, seconds
   real(c_double) start_time, stop_time, calc_time
   !real(c_double), parameter :: pi = 2.0d0*dacos(0.0d0)
   real(c_double) :: pi = dacos(-1.0d0)

   !type(Input) INP
   !type(Output) OUTP

   real(c_double) Lz, Tend, Delta, dz, dt, Ic
   real(c_double), allocatable :: ZAxis(:), TAxis(:), ZBEG, ZEND
   integer(c_int) i, j, Ne, Nz, Nt, err_alloc, INTERVALT, INTERVALZ, OUTNz, OUTNt
   complex(c_double_complex), allocatable :: InitialField(:), OUTB(:, :), OUTCu(:, :)
   logical(c_bool) :: IFR = .false.

   namelist /param/ Ne, Lz, Tend, Delta, Ic, dz, dt

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
   INTERVALT = Nt/500
   INTERVALZ = Nz/500
   !INTERVALT = 1
   !INTERVALZ = 1

   if (INTERVALT < 1) then
      print *, 'Too small "Tend"'
      pause
      stop
   end if

   if (INTERVALZ < 1) then
      print *, 'Too small "Lz"'
      pause
      stop
   end if

   if (INTERVALT > 1 .and. INTERVALZ > 1) then
      OUTNt = Nt/INTERVALT
      OUTNz = Nz/INTERVALZ
   elseif (INTERVALT == 1 .and. INTERVALZ > 1) then
      OUTNt = Nt
      OUTNz = Nz/INTERVALZ
   elseif (INTERVALT > 1 .and. INTERVALZ == 1) then
      OUTNt = Nt/INTERVALT
      OUTNz = Nz
   else ! (INTERVALT == 1) .and. (INTERVALZ == 1) then
      OUTNt = Nt
      OUTNz = Nz
   end if

   ! Vydelit' pamyat' pod massivy vyvoda
   allocate (OUTB(0:OUTNz, 0:OUTNt), OUTCu(0:OUTNz, 0:OUTNt), stat=err_alloc)
   !allocate (OUTP%OUTB(0:OUTNz, 0:OUTNt), OUTP%OUTCu(0:OUTNz, 0:OUTNt), stat=err_alloc)
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
      where ((ZAxis .GT. ZBEG) .AND. (ZAxis .LT. ZEND)) InitialField = 0.001*sin(pi*(ZAxis - ZBEG)/(ZEND - ZBEG))**2; 
   end if

   !INP%Ne = Ne
   !INP%ZAxis(:) = ZAxis(:)
   !INP%TAxis(:) = TAxis(:)
   !INP%Delta = Delta
   !INP%Ic = Ic
   !INP%INTERVALT = INTERVALT
   !INP%INTERVALZ = INTERVALZ
   !INP%InitialField(:) = InitialField
   !INP%OUTNz = OUTNz
   !INP%Nz = Nz
   !INP%Nt = Nt
   !INP%dz = dz
   !INP%dt = dt

   !open (1, file='test.dat')
   !do i = 0, Nz
   !   write (1, '(1p2e17.7)') InitialField(i)
   !end do
   !close (1)
   !stop

   !call oro(INP, OUTP)
   call oro(Nz, Nt, Ne, INTERVALT, INTERVALZ, OUTNz, Delta, Ic, dz, dt, &
            ZAxis, TAxis, InitialField, OUTB, OUTCu)

   !open(1, file = 'test.dat', err = 101)
   !write(1, '(4I)', err = 103)  size(TAxis,1), INTERVALT, size(ZAxis,1), INTERVALZ
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

   open (1, file='br.dat')
   do i = 0, OUTNz
      do j = 0, OUTNt
!         write (1, '(1pe17.7,a,\)') dreal(OUTP%OUTB(i, j)), ' '
         write (1, fmt='(1x,1pe17.7,a)', advance="no") dreal(OUTB(i, j))
      end do
      write (1, *) ! Assumes default "advance='yes'".
!      write (1, '(/,\)')
   end do
   close (1)

   open (1, file='bi.dat')
   do i = 0, OUTNz
      do j = 0, OUTNt
!         write (1, '(1pe17.7,a,\)') dimag(OUTP%OUTB(i, j)), ' '
         write (1, fmt='(1x,1pe17.7,a)', advance="no") dimag(OUTB(i, j))
      end do
      write (1, *) ! Assumes default "advance='yes'".
!      write (1, '(/,\)')
   end do
   close (1)

   open (1, file='ir.dat')
   do i = 0, OUTNz
      do j = 0, OUTNt
!         write (1, '(1pe17.7,a,\)') dreal(OUTP%OUTCu(i, j)), ' '
         write (1, fmt='(1x,1pe17.7,a)', advance="no") dreal(OUTCu(i, j))
      end do
      write (1, *) ! Assumes default "advance='yes'".
!      write (1, '(/,\)')
   end do
   close (1)

   open (1, file='ii.dat')
   do i = 0, OUTNz
      do j = 0, OUTNt
!         write (1, '(1pe17.7,a,\)') dimag(OUTP%OUTCu(i, j)), ' '
         write (1, fmt='(1x,1pe17.7,a)', advance="no") dimag(OUTCu(i, j))
      end do
      write (1, *) ! Assumes default "advance='yes'".
!      write (1, '(/,\)')
   end do
   close (1)

   pause
   stop
101 print *, "error of file open"; pause; stop
102 print *, 'error of file read'; pause; stop
103 print *, 'error of file write'; pause; stop
end
