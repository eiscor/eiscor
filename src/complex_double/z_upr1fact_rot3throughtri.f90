#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_rot3throughtri
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine passes a rotation through an upper-triangular matrix 
! stored as the product of 2 sequences of Givens rotations.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  DIR             LOGICAL
!                    .TRUE.: pass rotation from left to right
!                    .FALSE.: pass rotation from right to left
!
!  C               REAL(8) array of dimension (6)
!                    array of generators for upper-triangular parts
!
!  B               REAL(8) array of dimension (6)
!                    array of generators for upper-triangular parts
!
!  G               REAL(8) array of dimension (3)
!                    generator for rotation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_rot3throughtri(DIR,C,B,G)

  implicit none
  
  ! input variables
  logical, intent(in) :: DIR
  real(8), intent(inout) :: C(6), B(6), G(3)
 
  ! compute variables
  logical :: SYM
  real(8) :: T(3)

  ! initialize SYM
  if ((C(1).EQ.B(1)).AND.(C(2).EQ.-B(2)).AND.(C(3).EQ.-B(3)).AND. &
             (C(4).EQ.B(4)).AND.(C(5).EQ.-B(5)).AND.(C(6).EQ.-B(6))) then
    SYM = .TRUE.
  else
    SYM = .FALSE.
  end if 

  ! L2R
  if (DIR) then
  
    ! through C and B, symmetric case
    if (SYM) then
    
      ! store G
      T = G

      ! through C
      call z_rot3_turnover(C(1:3),C(4:6),T)

      ! update B
      B(1) = C(1)
      B(2) = -C(2)
      B(3) = -C(3)
      B(4) = C(4)
      B(5) = -C(5)
      B(6) = -C(6)

    ! general case
    else

      ! though C
      call z_rot3_turnover(C(1:3),C(4:6),G)
      
      ! through B
      call z_rot3_turnover(B(4:6),B(1:3),G)
    
    end if
  
  ! R2L
  else
  
    ! through C and B, symmetric case
    if (SYM) then
    
      ! store G
      T = G

      ! through B
      call z_rot3_turnover(B(1:3),B(4:6),T)

      ! update B
      C(1) = B(1)
      C(2) = -B(2)
      C(3) = -B(3)
      C(4) = B(4)
      C(5) = -B(5)
      C(6) = -B(6)

    ! general case
    else
  
      ! through B
      call z_rot3_turnover(B(1:3),B(4:6),G)
      
      ! through C
      call z_rot3_turnover(C(4:6),C(1:3),G)
      
    end if
    
  end if

end subroutine z_upr1fact_rot3throughtri
