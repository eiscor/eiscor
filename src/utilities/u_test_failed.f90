#include "eiscor.h"
subroutine u_test_failed(LINENUM) 

  implicit none

  ! input variables
  integer, intent(in) :: LINENUM

  ! print failure
  if (VERBOSE) then
     write(STDERR,*) ''
  end if
  write(STDERR,'(a,I4)') 'FAILED on line: ',LINENUM
  stop
  
end subroutine u_test_failed
