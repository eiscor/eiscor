#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot4_rot2fuse
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the product of two complex Givens rotations and 
! and stores the output in one of the input arrays.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  DIR             LOGICAL
!                    .TRUE.: A = B*A
!                    .FALSE.: A = A*B
!
!  A               REAL(8) array of dimension (4)
!                    generators for rotation
!
!  B               REAL(8) array of dimension (2)
!                    generators for rotation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot4_rot2fuse(DIR,A,B)

  implicit none

  ! input variables
  logical, intent(in) :: DIR
  real(8), intent(inout) :: A(4), B(2)

  ! compute variables
  real(8) :: temp(4)

  ! pad with a zero
  temp(1) = B(1)
  temp(2) = 0d0
  temp(3) = B(2)
  temp(4) = 0d0

  ! call z_rot4_rot4fuse
  call z_rot4_rot4fuse(DIR,A,temp)

end subroutine z_rot4_rot2fuse
