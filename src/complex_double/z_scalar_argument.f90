#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_scalar_argument
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the argument of a complex number A + BI.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  A,B            REAL(8)
!                    real and imaginary part of complex number
!
! OUTPUT VARIABLES:
!
!  ARG            REAL(8)
!                    arg of A + Bi in the interval [0,2pi)
!
!  INFO           INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies A is invalid
!                    INFO = -2 implies B is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_scalar_argument(A,B,ARG,INFO)

  ! input variables
  real(8), intent(in) :: A, B
  real(8), intent(inout) :: ARG
  integer, intent(inout) :: INFO
  
  ! compute variables
  real(8), parameter :: PI = 3.14159265358979323846264338327950d0
  
  ! initialize INFO
  INFO = 0  
  
  ! check error in debug mode
  if (DEBUG) then
  
    ! check A
    call d_scalar_check(A,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"A is invalid",INFO,-1)
      return
    end if

    ! check B
    call d_scalar_check(B,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"B is invalid",INFO,-2)
      return
    end if

  end if
    
  ! compute arg 1
  if ((abs(A).EQ.0d0).AND.(abs(B).EQ.0d0)) then
    ARG = 0d0
   
  ! abs(A) > abs(B)
  else if (abs(A) > abs(B)) then
  
    ! compute ARG between [-pi/2,pi/2]
    ARG = atan(abs(B/A))
  
  ! abs(A) < abs(B)
  else
  
    ! compute ARG between [-pi/2,pi/2]
    ARG = atan(abs(A/B))
    ARG = PI/2d0 - ARG
    
  end if
  
  ! correct for [0,2pi)
  ! second quadrant
  if ((B >= 0).AND.(A < 0)) then
    ARG = PI-ARG
     
  ! third quadrant
  else if ((B < 0).AND.(A < 0)) then
    ARG = PI+ARG
      
  ! fourth quadrant
  else if ((B < 0).AND.(A > 0)) then
    ARG = 2d0*PI-ARG
  end if
      
end subroutine z_scalar_argument
