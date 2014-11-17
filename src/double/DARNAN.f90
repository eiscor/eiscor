#include "eiscor.h"
!
! DARNAN (Double Auxiliary Routine NAN check)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks a real number to see if it is a NAN.
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
!                    INFO equal to 1 implies NUM is a NAN.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DARNAN(NUM,INFO)

  implicit none
  
  ! input variables
  real(8), intent(in) :: NUM
  integer, intent(inout) :: INFO
  
  ! initialize INFO
  INFO = 0
  
  ! check self equality
  if (NUM.NE.NUM) then
    INFO = 1
  end if

end subroutine DARNAN
