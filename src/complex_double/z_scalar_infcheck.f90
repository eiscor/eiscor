#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_scalar_infcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks a complex number to see if it is a INF.
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
!                    INFO = -1 implies NUM is an INF.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_scalar_infcheck(NUM,INFO)

  implicit none
  
  ! input variables
  complex(8), intent(in) :: NUM
  integer, intent(inout) :: INFO
  
  ! parameter
  real(8) :: temp
  real(8), parameter :: infdef = huge(temp)
  
  ! initialize INFO
  INFO = 0
  
  ! check magnitude of real and imaginary parts
  if ((abs(dble(NUM))>infdef).OR.(abs(aimag(NUM))>infdef)) then
    INFO = -1
  end if

end subroutine z_scalar_infcheck
