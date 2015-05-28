#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_correctend 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine adjusts either the top 2 or bottom 2 Givens rotations
! of an extended-Hessenberg Q so that core chasing can be done with
! minimal flops.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  TOP             LOGICAL
!                    .TRUE.: adjust the top
!                    .FALSE.: adjust the bottom
!
!  P               LOGICAL array of dimension (2)
!                    position flags for Q
!
!  Q               REAL(8) array of dimension (8)
!                    array of generators for givens rotations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_correctend(TOP,P,Q)

  implicit none
  
  ! input variables
  logical, intent(in) :: TOP
  logical, intent(inout) :: P(2)
  real(8), intent(inout) :: Q(8)
  
  ! compute variables
  real(8) :: nrm
  complex(8) :: temp
  
  ! adjust top if necessary
  if (TOP) then
  
    ! P(2) == .FALSE. and P(1) == .TRUE.
    if ((.NOT.P(2)).AND.P(1)) then 

      temp = cmplx(Q(1),-Q(2),kind=8)*cmplx(Q(7),Q(8),kind=8)
      call z_rot4_vec4gen(Q(5),Q(6),dble(temp),aimag(temp) &
      ,Q(5),Q(6),Q(7),Q(8),nrm)
      P(1) = P(2)

    ! P(2) == .TRUE. and P(1) == .FALSE.
    else if (P(2).AND.(.NOT.P(1))) then 

      temp = cmplx(Q(1),Q(2),kind=8)*cmplx(Q(7),Q(8),kind=8)
      call z_rot4_vec4gen(Q(5),Q(6),dble(temp),aimag(temp) &
      ,Q(5),Q(6),Q(7),Q(8),nrm)
      P(1) = P(2)

    end if

  ! adjust bottom if necessary
  else
  
    ! P(1) == .FALSE. and P(2) == .TRUE.
    if ((.NOT.P(1)).AND.P(2)) then 

      temp = cmplx(Q(1),Q(2),kind=8)*cmplx(Q(7),Q(8),kind=8)
      call z_rot4_vec4gen(Q(5),Q(6),dble(temp),aimag(temp) &
      ,Q(5),Q(6),Q(7),Q(8),nrm)
      P(2) = P(1)

    ! P(1) == .TRUE. and P(2) == .FALSE.
    else if (P(1).AND.(.NOT.P(2))) then 

      temp = cmplx(Q(1),-Q(2),kind=8)*cmplx(Q(7),Q(8),kind=8)
      call z_rot4_vec4gen(Q(5),Q(6),dble(temp),aimag(temp) &
      ,Q(5),Q(6),Q(7),Q(8),nrm)
      P(2) = P(1)

    end if

  end if

end subroutine z_upr1fact_correctend
