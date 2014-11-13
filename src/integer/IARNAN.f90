#include "eiscor.h"
!
! IARNAN (Integer Auxiliary Routine NAN check)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks an integer number to see if it is a NAN.
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
!                    INFO equal to 1 implies NUM is a NAN.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine IARNAN(NUM,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: NUM
  integer, intent(inout) :: INFO
  
  ! initialize INFO
  INFO = 0
  
  ! check self equality
  if (NUM.NE.NUM) then
    INFO = 1
  end if

end subroutine IARNAN
