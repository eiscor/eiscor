#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_mergebulge 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine fuses a Givens rotation represented by 3 real numbers: 
! the real and imaginary parts of a complex cosine and a scrictly real 
! sine into the unitary extended hessenberg part of a upr1 pencil.
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
  
  ! fusion at top from left
  if (TOP.AND..NOT.P(1)) then
  
    ! fuse with Q 
    call z_rot4_rot3fuse(.TRUE.,Q(1:4),G) 

  ! fusion at top from right
  else if (TOP.AND.P(1)) then
  
    ! fuse with Q 
    call z_rot4_rot3fuse(.FALSE.,Q(1:4),G) 

  ! fusion at bottom from right
  else if (.NOT.TOP.AND..NOT.P(1)) then
    
    ! next Q is not in the way
    if (P(2)) then
  
      ! fuse with Q 
      call z_rot4_rot3fuse(.FALSE.,Q(1:4),G) 

    ! next Q is in the way
    else
  
      ! temporary storage
      B(1:2) = G(1:2)
      B(3) = G(3)*Q(5)
      B(4) = G(3)*Q(6)

      ! fuse with Q 
      call z_rot4_rot4fuse(.FALSE.,Q(1:4),B) 

    end if

  ! fusion at bottom from left
  else
    
    ! next Q is not in the way
    if (.NOT.P(2)) then
  
      ! fuse with Q 
      call z_rot4_rot3fuse(.TRUE.,Q(1:4),G) 

    ! next Q is in the way
    else
  
      ! temporary storage
      B(1:2) = G(1:2)
      B(3) = G(3)*Q(5)
      B(4) = -G(3)*Q(6)

      ! fuse with Q 
      call z_rot4_rot4fuse(.TRUE.,Q(1:4),B) 

    end if

  end if

end subroutine z_upr1fact_mergebulge
