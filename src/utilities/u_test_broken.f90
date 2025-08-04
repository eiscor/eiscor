#include "eiscor.h"

! Same as u_test_failed.f90 but returns a zero exit status

subroutine u_test_broken(LINENUM)

  implicit none

  ! input variables
  integer, intent(in) :: LINENUM

  ! print failure
  write(STDERR,'(a,I4)') 'BROKEN on line: ',LINENUM
  stop

end subroutine u_test_broken
