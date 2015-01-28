#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_rot3throughtri
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine passes a rotation through an upper-triangular matrix 
! stored as the product of a diagonal matrix and 2 sequences of 
! rotations.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  DIR             LOGICAL
!                    .TRUE.: pass rotation from left to right
!                    .FALSE.: pass rotation from right to left
!
!  FLAG            LOGICAL
!                    .TRUE.: pass rotation through C and B
!                    .FALSE.: ignore C and B
!
!  D               REAL(8) array of dimension (4)
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
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
subroutine z_upr1fact_rot3throughtri(DIR,FLAG,D,C,B,G)

  implicit none
  
  ! input variables
  logical, intent(in) :: DIR, FLAG
  real(8), intent(inout) :: D(4), C(6), B(6), G(3)
  
  ! L2R
  if (DIR) then
  
    ! through D
    call z_rot3_swapdiag(DIR,D,G)
    
    ! through C and B
    if (FLAG) then
    
      ! though C
      call z_rot3_turnover(C(1:3),C(4:6),G)
      
      ! through B
      call z_rot3_turnover(B(4:6),B(1:3),G)
    
    end if
  
  ! R2L
  else
  
    ! through C and B
    if (FLAG) then
  
      ! through B
      call z_rot3_turnover(B(1:3),B(4:6),G)
      
      ! through C
      call z_rot3_turnover(C(4:6),C(1:3),G)
      
    end if
    
    ! through D
    call z_rot3_swapdiag(DIR,D,G)
  
  end if

end subroutine z_upr1fact_rot3throughtri
