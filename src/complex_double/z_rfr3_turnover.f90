#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rfr3_turnover 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rfr3_turnover(W,CC,SS,U,VV,RHO)

  implicit none
  
  ! input variables
  real(8), intent(inout) :: CC, SS, VV
  complex(8), intent(inout) :: W, U
  complex(8), intent(in) :: RHO
  
  ! compute variables
  complex(8) :: z, Uh, Wh
  real(8) :: nn, zinf, zz, xx, CCh, SSh, VVh

  ! z and zz
  z = U + W
  zinf = max(abs(dble(z)),abs(aimag(z)))
  zz = dble(z)**2 + aimag(z)**2

  ! new nn
  nn = CC*zz + VV

  ! nn > 0
  if ( nn > 0d0 ) then
    CCh = CC*zz/nn
    SSh = VV/nn
    Uh  = conjg(RHO)*(SS*U - CC*W)
    VVh = nn*SS
    ! zz > 0
    if ( zz > 0d0 ) then
      Wh = -RHO*conjg(W)*z**2/zz
    else
      Wh = cmplx(1d0,0d0,kind=8)
    end if
  else
    CCh = 1d0
    SSh = 0d0
    Uh  = conjg(RHO)*U
    VVh = 0d0
    Wh  = -RHO*U
  end if

  ! ensure normality
  xx = dble(Uh)**2 + aimag(Uh)**2 + VVh - 1d0
  U  = Uh*(1d0 - xx/2d0)
  VV = VVh/(1d0+xx)
  CC = CCh
  SS = SSh
  xx = dble(W)**2 + aimag(W)**2 - 1d0
  W  = Wh*(1d0 - xx/2d0)

end subroutine z_rfr3_turnover
