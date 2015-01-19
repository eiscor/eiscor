#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_scalar_nancheck
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
!                    INFO = -1 implies NUM is a NAN.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_scalar_nancheck(NUM,INFO)

  implicit none
  
  ! input variables
  complex(8), intent(in) :: NUM
  integer, intent(inout) :: INFO
  
  ! initialize INFO
  INFO = 0
  
  ! check self equality
  if (NUM.NE.NUM) then
    INFO = -1
  end if

end subroutine z_scalar_nancheck
