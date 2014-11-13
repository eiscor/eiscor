#include "eiscor.h"
!
! IARINF (Integer Auxiliary Routine INF check)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks an integer number to see if it is a INF.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  NUM             INTEGER
!                    integer number to be checked
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                    INFO equal to 1 implies NUM is an INF.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine IARINF(NUM,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: NUM
  integer, intent(inout) :: INFO
  
  ! parameter
  integer :: temp
  integer, parameter :: infdef = huge(temp)
  
  ! initialize INFO
  INFO = 0
  
  ! check magnitude of real and imaginary parts
  if (abs(NUM)>infdef) then
    INFO = 1
  end if

end subroutine IARINF
