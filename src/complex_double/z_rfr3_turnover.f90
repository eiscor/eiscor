#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rfr3_turnover 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  W               COMPLEX(8) 
!                    array of generators for givens rotations
!
!  CC, SS          REAL(8) 
!                    array of generators for givens rotations
!
!  U               COMPLEX(8) 
!                    array of generators for givens rotations
!
!  VV              REAL(8) 
!                    array of generators for givens rotations
!
!  RHO             COMPLEX(8) 
!                    array of generators for givens rotations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rfr3_turnover(W,CC,SS,U,VV,RHO)

  implicit none
  
  ! input variables
  real(8), intent(inout) :: CC, SS, VV
  complex(8), intent(inout) :: W, U
  complex(8), intent(in) :: RHO
  
  ! compute variables
  real(8) :: nn, zz, xx
  complex(8) :: z
  complex(8) :: t, uold, D1(3), D2(3)

  ! store old SS and U
  xx = SS
  z = U
  uold = U
  D1(1) = -conjg(RHO)*(CC*W - SS*U)
  D1(2) = -RHO*(CC*conjg(W)*U-SS)
  D1(3) = conjg(U) 

  ! new U
  U = -conjg(RHO)*(CC*W - SS*U)

  ! new W, CC and SS
  z = z + W + 10d0*(EISCOR_DBL_EPS)
  zz = dble(z)**2 + aimag(z)**2
  nn = VV + CC*zz
  W = -RHO*conjg(W)*((z*z)/zz)
  CC = CC*zz/nn
  SS = VV/nn
  
  ! new VV
  VV = xx*nn

  ! ensure normality
  xx = dble(U)**2 + aimag(U)**2 + VV
  U = U*(3d0 - xx)/2
  VV = VV/xx
  xx = dble(W)**2 + aimag(W)**2
  W = W*(3d0 - xx)/2
  
  D2(1) = U
  D2(2) = RHO*SS - conjg(RHO)*CC*W*conjg(D1(1)) 
  D2(3) = conjg(RHO)*SS*conjg(D1(1)) - RHO*CC*conjg(W) 
if ((zz < 1d-30)) then
print*,abs(D1)
print*,abs(D2)
end if

end subroutine z_rfr3_turnover
