#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rfr3_turnover 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a square root free turnover. 
!
! | cg        -s |                | c       -s | | -conj(rho)     0 |
! |  s  conj(cg) | | u       -v | | s  conj(c) | |          0  -rho |
!                  | v  conj(u) |
!
!
!                                = 
!
!                                     | u       -v |       
! | c       -s | | -conj(rho)     0 | | v  conj(u) | | cg        -s |
! | s  conj(c) | |          0  -rho |                |  s  conj(cg) |
!
! Due to the symmetry it is not necessary to store the values c, s
! or g explicitly. The following variables take their place:
!
! CC = |c|^2
! SS = s^2
!  W = gc^2/|c|^2
! VV = v^2
!
! The input W, CC, SS, U, VV and RHO must satisfy the following:
!
!    CC + SS = 1
!        |W| = 1
! |U|^2 + VV = 1
!      |RHO| = 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  W               COMPLEX(8) 
!                    unimodular complex number
!
!  CC, SS          REAL(8) 
!                    moduli squared of cosine and sine
!
!  U               COMPLEX(8) 
!                    complex component of Givens rotation
!
!  VV              REAL(8) 
!                    real component squared of Givens rotation
!
!  RHO             COMPLEX(8) 
!                    unimodular shift
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
  complex(8) :: z, uold

  ! store old SS and U
  xx = SS
  uold = U

  ! new U
  U = (CC*W - SS*U)

  ! new W, CC and SS
  z = uold + W
  zz = dble(z)**2 + aimag(z)**2
  CC = CC*zz
  nn = VV + CC
  if ( CC.EQ.0d0 ) then
    CC = 0d0
    SS = 1d0
    W = RHO*W
  else
    CC = CC/nn
    SS = VV/nn
    W = -RHO*(SS*U + uold)/CC
  end if

  ! new U
  U = -conjg(RHO)*U

  ! new VV
  VV = xx*nn

  ! ensure normality
  xx = dble(U)**2 + aimag(U)**2 + VV
  U = U*(3d0 - xx)/2
  VV = VV/xx
  xx = dble(W)**2 + aimag(W)**2
  W = W*(3d0 - xx)/2

end subroutine z_rfr3_turnover
