#include "eiscor.h"
!
! d_2Darray_check
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
!                    INFO = 0 implies A contains no INF nor NAN.
!                    INFO = -1 implies N is invalid.
!                    INFO = -2 implies M is invalid.
!                    INFO = -3 implies A contains an INF or NAN.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_2Darray_check(M,N,A,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: M, N
  real(8), intent(in) :: A(M,N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii, jj
  
  ! initialize INFO
  INFO = 0
  
  if (DEBUG) then
    ! check M
    if (M<1) then
      INFO = -1
      call u_infocode_check(__FILE__,__LINE__,"M must be at least 1",INFO,INFO)
      return
    end if

    ! check N
    if (N<1) then
      INFO = -2
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 1",INFO,INFO)
      return
    end if
  end if

  ! check array
  do ii = 1,M
    do jj = 1,N
  
      ! check for NAN
      call d_scalar_nancheck(A(ii,jj),INFO)
      if (INFO.NE.0) then
        INFO = -3
        return
      end if
    
      ! check for INF
      call d_scalar_infcheck(A(ii,jj),INFO)      
      if (INFO.NE.0) then
        INFO = -3
        return
      end if
    
    end do
  end do

end subroutine d_2Darray_check
