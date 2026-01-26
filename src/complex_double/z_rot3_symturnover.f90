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
! The input G, C, S, U, V and RHO must satisfy the following:
!
!         |G| = 1
! |C|^2 + S^2 = 1
! |U|^2 + V^2 = 1
!       |RHO| = 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  G               COMPLEX(8) 
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
subroutine z_rot3_symturnover(G,C,SIG,S,U,V,RHO)

  implicit none
  
  ! input variables
  real(8), intent(inout) :: C, S, V
  complex(8), intent(inout) :: G, SIG, U
  complex(8), intent(in) :: RHO
  
  ! compute variables
  real(8) :: n, xx, pp
  complex(8) :: p, w, uold
real(8) :: modz
complex(8) :: z

  ! W
  w = G*SIG*SIG

  ! store old V 
  xx = V
uold = U

  ! p, pp and n
  z = w + U
  modz = abs(z)
  n = sqrt(C**2*modz**2 + V**2)

  ! new U and V
  U = -conjg(RHO)*(C**2*z - U)
  V = S*n

  ! new C and S
  if ( C*modz.EQ.0d0 ) then
    C = 1d0
    S = 0d0
    SIG = cmplx(1d0,0d0,kind=8)
!  else if (xx.EQ.0d0) then
!    C = 1d0
!    S = 0d0
!    SIG = -RHO*sqrt(conjg(G)*uold)
  else
    C = C*modz/n
    S = xx/n
w = conjg(G)*(uold + xx**2/modz**2*z)
if (dble(w) > 0) then
w = 1d0+w
SIG = -RHO*w/abs(w)
else
w = cmplx(0d0,1d0,kind=8)*(w-1d0)
SIG = -RHO*w/abs(w)
end if
  end if

  ! new G
  G = -conjg(RHO)*G

  ! ensure normality of U and V
  xx = dble(U)**2 + aimag(U)**2 + V**2
  U = 5d-1*U*(3d0 - xx)
  V = 5d-1*V*(3d0 - xx)

  xx = dble(G)**2 + aimag(G)**2
  G = 5d-1*G*(3d0 - xx)

  xx = dble(SIG)**2 + aimag(SIG)**2
  SIG = 5d-1*SIG*(3d0 - xx)

  xx = C**2 + S**2
  C = 5d-1*C*(3d0 - xx)
  S = 5d-1*S*(3d0 - xx)

end subroutine z_rot3_symturnover
