#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot4_rot4fuse
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
!  B               REAL(8) array of dimension (4)
!                    generators for rotation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot4_rot4fuse(DIR,A,B)

  implicit none

  ! input variables
  logical, intent(in) :: DIR
  real(8), intent(inout) :: A(4), B(4)

  ! compute variables
  real(8) :: cr, ci, sr, si, nrm

  ! fuse from left
  if(DIR)then

    ! compute new generators
    cr = B(1)*A(1) - B(2)*A(2) - B(3)*A(3) - B(4)*A(4)
    ci = B(1)*A(2) + B(2)*A(1) - B(3)*A(4) + B(4)*A(3) 
    sr = B(3)*A(1) - B(4)*A(2) + B(1)*A(3) + B(2)*A(4)
    si = B(3)*A(2) + B(4)*A(1) + B(1)*A(4) - B(2)*A(3)
    
    ! store in A
    call z_rot4_vec4gen(cr,ci,sr,si,A(1),A(2),A(3),A(4),nrm)
    
  ! fuse from right
  else

    ! compute new generators
    cr = A(1)*B(1) - A(2)*B(2) - A(3)*B(3) - A(4)*B(4)
    ci = A(1)*B(2) + A(2)*B(1) - A(3)*B(4) + A(4)*B(3) 
    sr = A(3)*B(1) - A(4)*B(2) + A(1)*B(3) + A(2)*B(4)
    si = A(3)*B(2) + A(4)*B(1) + A(1)*B(4) - A(2)*B(3)
    
    ! store in A
    call z_rot4_vec4gen(cr,ci,sr,si,A(1),A(2),A(3),A(4),nrm)
    
  end if

end subroutine z_rot4_rot4fuse
