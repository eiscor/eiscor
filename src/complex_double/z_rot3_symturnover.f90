#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot3_symturnover 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a square root free turnover. 
!
! | c       -s |                | cg        -s | | -conj(rho)     0 |
! | s  conj(c) | | u       -v | |  s  conj(cg) | |          0  -rho |
!                | v  conj(u) |
!
!
!                                = 
!
!                                       | u       -v |       
! | cg        -s | | -conj(rho)     0 | | v  conj(u) | | c       -s |
! |  s  conj(cg) | |          0  -rho |                | s  conj(c) |
!
! The input W, C, S, U, V and RHO must satisfy the following:
!
!         |W| = 1
! |C|^2 + S^2 = 1
! |U|^2 + V^2 = 1
!       |RHO| = 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  W               COMPLEX(8) 
!                    unimodular complex number
!
!  C               COMPLEX(8) 
!                    real component of Givens rotation
!
!  S               REAL(8) 
!                    moduli squared of cosine and sine
!
!  U               COMPLEX(8) 
!                    complex component of Givens rotation
!
!  V               REAL(8) 
!                    real component of Givens rotation
!
!  RHO             COMPLEX(8) 
!                    unimodular shift
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3_symturnover(W,C,S,U,V,RHO)

  implicit none
  
  ! input variables
  real(8), intent(inout) :: C, S, V
  complex(8), intent(inout) :: W, U
  complex(8), intent(in) :: RHO
  
  ! compute variables
  complex(8) :: z, Wh, Uh
  real(8) :: zz, zabs, n, xx, Ch, Sh, Vh

  ! z and zz
  z = U + W
  zz = dble(z)**2 + aimag(z)**2

  ! new n
  n = sqrt(C**2*zz + V**2)
  zabs = sqrt(zz)

  ! n > 0
  if ( n > 0d0 ) then
    Ch = C*zabs/n
    Sh = V/n
    Uh  = conjg(RHO)*(S**2*U - C**2*W)
    Vh = n*S
    ! zz > 0
    if ( zz > 0d0 ) then
!      Wh = -RHO*conjg(W)*z**2/zz
      Wh = -RHO*(U + z*(V/zabs)**2)
    else
      Wh = cmplx(1d0,0d0,kind=8)
    end if
  else
    Ch = 1d0
    Sh = 0d0
    Uh  = conjg(RHO)*U
    Vh = 0d0
    Wh  = -RHO*U
  end if

  ! ensure normality
  xx = dble(Uh)**2 + aimag(Uh)**2 + Vh**2 - 1d0
  U  = Uh*(1d0 - xx/2d0)
  V  = Vh*(1d0 - xx/2d0)
  C  = Ch
  S  = Sh
  xx = dble(Wh)**2 + aimag(Wh)**2 - 1d0
  W  = Wh*(1d0 - xx/2d0)

end subroutine z_rot3_symturnover
