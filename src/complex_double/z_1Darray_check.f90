#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_1Darray_check 
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
subroutine z_1Darray_check(N,A,INFO)

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
    call z_scalar_nancheck(A(ii),INFO)
    if (INFO.NE.0) then
      return
    end if
    
    ! check for INF
    call z_scalar_infcheck(A(ii),INFO)
    if (INFO.NE.0) then
      return
    end if
    
  end do

end subroutine z_1Darray_check
