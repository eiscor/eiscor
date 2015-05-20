#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! l_upr1fact_upperhess
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This function returns the new position flag for a single shift 
! iteration on a twisted upr1 pencil based on the position flags from
! the current factorization. The returned flag is always .FALSE. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function l_upr1fact_upperhess(N,P)

  implicit none
  
  ! return type
  logical :: l_upr1fact_upperhess
  
  ! input variables
  integer, intent(in) :: N
  logical, intent(in) :: P(N-2)
  
  ! set output
  l_upr1fact_upperhess = .FALSE.

end function l_upr1fact_upperhess
