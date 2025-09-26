#include "eiscor.h"
subroutine u_test_skipped()

  implicit none

  ! print skipped message
  write(STDERR,'(a)') 'SKIPPED due do PRNG changes'
  stop

end subroutine u_test_skipped
