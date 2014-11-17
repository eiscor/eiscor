#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZARNAN (Zomplex Auxiliary Routine NAN check)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks a complex number to see if it is a NAN.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  NUM             COMPLEX(8) 
!                    complex number to be checked
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                    INFO equal to 1 implies NUM is a NAN.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZARNAN(NUM,INFO)

  implicit none
  
  ! input variables
  complex(8), intent(in) :: NUM
  integer, intent(inout) :: INFO
  
  ! initialize INFO
  INFO = 0
  
  ! check self equality
  if (NUM.NE.NUM) then
    INFO = 1
  end if

end subroutine ZARNAN
