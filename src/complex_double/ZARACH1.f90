#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZARACH1 (Zomplex Auxiliary Routine Array CHeck 1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks a one dimensional complex array for INFs and NANs.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    positive integer, dimension of A
!
!  A               COMPLEX(8) array of dimension (N)
!                    complex array to be checked
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                    INFO equal to 1 implies A contains an INF or NAN.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZARACH1(N,A,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  complex(8), intent(in) :: A(N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii
  
  ! initialize INFO
  INFO = 0
  
  ! check array
  do ii = 1,N
  
    ! check for NAN
    call ZARNAN(A(ii),INFO)
    if (INFO.NE.0) then
      return
    end if
    
    ! check for INF
    call ZARINF(A(ii),INFO)
    if (INFO.NE.0) then
      return
    end if
    
  end do

end subroutine ZARACH1
