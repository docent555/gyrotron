module util

   interface assert_eq
      module procedure assert_eq2, assert_eq3, assert_eq4, assert_eqn
   end interface

contains

!bl
   function assert_eq2(n1, n2, string)
      character(len=*), intent(in) :: string
      integer, intent(in) :: n1, n2
      integer :: assert_eq2
      if (n1 == n2) then
         assert_eq2 = n1
      else
         write (*, *) 'error: an assert_eq failed with this tag:', &
            string
         stop 'program terminated by assert_eq2'
      end if
   end function assert_eq2
!bl
   function assert_eq3(n1, n2, n3, string)
      character(len=*), intent(in) :: string
      integer, intent(in) :: n1, n2, n3
      integer :: assert_eq3
      if (n1 == n2 .and. n2 == n3) then
         assert_eq3 = n1
      else
         write (*, *) 'error: an assert_eq failed with this tag:', &
            string
         pause
         stop 'program terminated by assert_eq3'
      end if
   end function assert_eq3
!bl
   function assert_eq4(n1, n2, n3, n4, string)
      character(len=*), intent(in) :: string
      integer, intent(in) :: n1, n2, n3, n4
      integer :: assert_eq4
      if (n1 == n2 .and. n2 == n3 .and. n3 == n4) then
         assert_eq4 = n1
      else
         write (*, *) 'error: an assert_eq failed with this tag:', &
            string
         pause
         stop 'program terminated by assert_eq4'
      end if
   end function assert_eq4
!bl
   function assert_eqn(nn, string)
      character(len=*), intent(in) :: string
      integer, dimension(:), intent(in) :: nn
      integer :: assert_eqn
      if (all(nn(2:) == nn(1))) then
         assert_eqn = nn(1)
      else
         write (*, *) 'error: an assert_eq failed with this tag:', &
            string
         pause
         stop 'program terminated by assert_eqn'
      end if
   end function assert_eqn
   subroutine error(string)
      character(len=*), intent(in) :: string
      write (*, *) 'error: ', string
      pause
      stop 'program terminated by error'
   end subroutine error
end module util
