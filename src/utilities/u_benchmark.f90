#include "eiscor.h"
subroutine u_benchmark(TIME,ACC) 

  implicit none

  ! input variables
  real(8), intent(in) :: TIME, ACC

  ! print failure
  write(STDERR,'(a,E10.5,a,E10.5,a)') 'TIME ',TIME,' secs, ACC ',ACC, ' points'
  
end subroutine u_benchmark