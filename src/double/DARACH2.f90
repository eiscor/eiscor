#include "eiscor.h"
!
! DARACH2 (Double Auxiliary Routine Array CHeck 2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks a one dimensional complex array for INFs and NANs.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  M, N            INTEGER
!                    positive integers, dimensions of A
!
!  A               REAL(8) array of dimension (M,N)
!                    real array to be checked
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                    INFO equal to -1 implies A contains an INF or NAN.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DARACH2(M,N,A,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: M, N
  real(8), intent(in) :: A(M,N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii, jj
  
  ! initialize INFO
  INFO = 0
  
  ! check array
  do ii = 1,M
    do jj = 1,N
  
      ! check for NAN
      call DARNAN(A(ii,jj),INFO)
      if (INFO.NE.0) then
        return
      end if
    
      ! check for INF
      call DARINF(A(ii,jj),INFO)
      if (INFO.NE.0) then
        return
      end if
    
    end do
  end do

end subroutine DARACH2
