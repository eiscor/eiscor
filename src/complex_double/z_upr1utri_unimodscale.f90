#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1utri_unimodscale
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine scales either a row or column of a unitary plus rank
! one upper-triangular matrix by a complex unimodular number.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  ROW             LOGICAL
!                    .TRUE.: scale row
!                    .FALSE.: scale column
!
!  D               REAL(8) array of dimension (2)
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C               REAL(8) array of dimension (3)
!                    array of generators for upper-triangular parts
!
!  B               REAL(8) array of dimension (3)
!                    array of generators for upper-triangular parts
!
!  SCL             COMPLEX(8) 
!                    unimodular scalar
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1utri_unimodscale(ROW,D,C,B,SCL)

  implicit none
  
  ! input variables
  logical, intent(in) :: ROW
  real(8), intent(inout) :: D(2), C(3), B(3)
  complex(8), intent(in) :: SCL
 
  ! compute variables
  logical :: SYM
  real(8) :: nrm
  complex(8) :: temp

  ! check to see if C and B are symmetric. if so we can avoid updating them
  if ((C(1).EQ.B(1)).AND.(C(2).EQ.-B(2)).AND.(C(3).EQ.-B(3))) then
    SYM = .TRUE.
  else
    SYM = .FALSE.
  end if 

  ! update D regardless of row or column
  temp = SCL*cmplx(D(1),D(2),kind=8)
  call d_rot2_vec2gen(dble(temp),aimag(temp),D(1),D(2),nrm)
  
  ! update B and C if SYM == .FALSE.
  if ((.NOT.SYM).AND.(.NOT.ROW)) then

      ! update B
    temp = SCL*cmplx(B(1),B(2),kind=8)
    call d_rot2_vec2gen(dble(temp),aimag(temp),B(1),B(2),nrm)
    B(1:2) = B(1:2)*nrm
        
    ! update C
    temp = conjg(SCL)*cmplx(C(1),C(2),kind=8)
    call d_rot2_vec2gen(dble(temp),aimag(temp),C(1),C(2),nrm)
    C(1:2) = C(1:2)*nrm
  
  end if

end subroutine z_upr1utri_unimodscale
