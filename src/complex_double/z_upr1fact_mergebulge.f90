#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_mergebulge 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! There are only 4 possibilities for this fusion:
!
! 1)    top, left  <=>      TOP.AND..NOT.P(1)
!
! 2)    top, right <=>      TOP.AND.P(1)
!
! 3) bottom, left  <=> .NOT.TOP.AND.P(1)
!
! 4) bottom, right <=> .NOT.TOP.AND..NOT.P(1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  TOP             LOGICAL
!                    .TRUE.: fuse at top
!                    .FALSE.: fuse at bottom
!
!  P               LOGICAL array of dimension (2)
!                    position flags for Q
!
!  Q               REAL(8) array of dimension (8)
!                    array of generators for givens rotations
!
!  G               REAL(8) array of dimension (3)
!                    generators for rotation that will be fused
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_mergebulge(TOP,P,Q,G)

  implicit none
  
  ! input variables
  logical, intent(in) :: TOP, P(2)
  real(8), intent(inout) :: Q(8)
  real(8), intent(in) :: G(3)
  
  ! compute variables
  real(8) :: B(4), nrm

! This is only here for debbugging purposes and should be removed when
! code is finalized.
! check to see if P(1) == P(2)
if (.NOT.(P(1).EQV.P(2))) then

    print*,""
    print*,""
    print*,"Inside z_upr1fact_mergebulge."
    print*,"  P(1) must equal P(2)!"
    print*,""
    print*,""

end if
  
  ! fusion at top from left
  if (TOP.AND..NOT.P(1)) then

    ! temporary storage
    B(1:2) = G(1:2)
    B(3) = G(3)*Q(1)
    B(4) = -G(3)*Q(2)

    ! fuse with Q 
    call z_rot4_rot4fuse(.TRUE.,Q(5:8),B) 

  ! fusion at top from right
  else if (TOP.AND.P(1)) then

    ! temporary storage
    B(1:2) = G(1:2)
    B(3) = G(3)*Q(1)
    B(4) = G(3)*Q(2)

    ! fuse with Q 
    call z_rot4_rot4fuse(.FALSE.,Q(5:8),B) 

  ! fusion at bottom from right
  else if (.NOT.TOP.AND..NOT.P(1)) then
  
    ! temporary storage
    B(1:2) = G(1:2)
    B(3) = G(3)*Q(5)
    B(4) = G(3)*Q(6)

    ! fuse with Q 
    call z_rot4_rot4fuse(.FALSE.,Q(1:4),B) 

  ! fusion at bottom from left
  else
  
    ! temporary storage
    B(1:2) = G(1:2)
    B(3) = G(3)*Q(5)
    B(4) = -G(3)*Q(6)

    ! fuse with Q 
    call z_rot4_rot4fuse(.TRUE.,Q(1:4),B) 

  end if

end subroutine z_upr1fact_mergebulge
