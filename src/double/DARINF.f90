#include "eiscor.h"
!
! DARINF (Double Auxiliary Routine INF check)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks a real number to see if it is a INF.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  NUM             REAL(8) 
!                    real number to be checked
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                    INFO equal to 1 implies NUM is an INF.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DARINF(NUM,INFO)

  implicit none
  
  ! input variables
  real(8), intent(in) :: NUM
  integer, intent(inout) :: INFO
  
  ! parameter
  real(8) :: temp
  real(8), parameter :: infdef = huge(temp)
  
  ! initialize INFO
  INFO = 0
  
  ! check magnitude 
  if (abs(NUM)>=infdef) then
    INFO = 1
  end if

end subroutine DARINF
