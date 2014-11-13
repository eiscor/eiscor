#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DARGTD (Double Auxiliary Routine Givens Through Diagonal)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine passes the generator for a Givens rotation represented
! by 2 real numbers: a strictly real cosine and a scrictly real sine 
! through a 2x2 complex diagonal matrix represented by 4 real numbers.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  JOB             CHARACTER
!                    'L': pass rotation from Left of diagonal to right
!                    'R': pass rotation from Right of diagonal to left
!
!  D               REAL(8) array of dimension (4)
!                    array of generators for complex diagonal matrix
!                    the entries of D must be 1 or -1
!
!  B               REAL(8) array of dimension (2)
!                    generator for a Givens rotation
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies JOB is invalid
!                    INFO = -2 implies D is invalid
!                    INFO = -3 implies B is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DARGTD(JOB,D,B,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: JOB
  integer, intent(inout) :: INFO
  real(8), intent(inout) :: D(4), B(2)
  
  ! initialize INFO
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then 
  
    ! check JOB
    if ((JOB.NE.'L').AND.(JOB.NE.'R')) then
      INFO = -1
      call UARERR(__FILE__,__LINE__,"JOB is invalid",INFO,-1)
      return
    end if
  
    ! check D
    if ((abs(D(1)).NE.1d0).AND.(abs(D(2)).NE.0d0)) then
      INFO = -2
      call UARERR(__FILE__,__LINE__,"D must be diagonal with entries 1 or -1",INFO,-2)
      return
    end if
    if ((abs(D(3)).NE.1d0).AND.(abs(D(4)).NE.0d0)) then
      INFO = -2
      call UARERR(__FILE__,__LINE__,"D must be diagonal with entries 1 or -1",INFO,-2)
      return
    end if
      
    ! check B
    call DARACH1(2,B,INFO)
    call UARERR(__FILE__,__LINE__,"B is invalid",INFO,-3)
    if (INFO.NE.0) then 
      return 
    end if 
  
  end if
  
  ! return if scalar
  if (D(1).EQ.D(3)) then
    return
  end if 
  
  !if JOB == L
  if (JOB.EQ.'L')then
    
    ! update B
    B(2) = -B(2)
     
  ! if JOB == R
  else
    
    ! update B
    B(2) = -B(2)
     
  end if
  
end subroutine DARGTD
