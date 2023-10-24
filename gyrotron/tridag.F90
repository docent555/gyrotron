        subroutine tridag(c, a, b, d, u)
           use, intrinsic :: iso_c_binding
           use util
           implicit none
           real(c_double), dimension(:), intent(in) :: c, a, b, d
           real(c_double), dimension(:), intent(out) :: u
           real(c_double), dimension(size(a)) :: gam
           integer(c_int) :: n, j
           real(c_double) :: bet
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
        end subroutine tridag
