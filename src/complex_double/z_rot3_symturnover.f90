#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot3_symturnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the standard symmetric turnover.
!
! NORMALIZE controls the final renormalization:
!
!       NORMALIZE(1) = .TRUE.  : renormalize U and V
!       NORMALIZE(2) = .TRUE.  : renormalize SIGMA
!       NORMALIZE(3) = .TRUE.  : renormalize GAMMA
!
! The three flags are independent.
!
! The input variables GAMMA, SIGMA, C, S, U, V and RHO determine the
! two phase parameters in the left/right rotations through
!
!       GAMMA*SIGMA**2.
!
! The scalar Z used below is the scaled quantity
!
!       Z = 1 + T,     T = conjg(GAMMA*SIGMA**2)*U,
!
! rather than the unscaled sum GAMMA*SIGMA**2 + U.  Since
! GAMMA*SIGMA**2 is unimodular, both quantities have the same modulus.
!
! If real(T) < 0, Z is computed with the cancellation-avoiding formula
!
!       Z = ( V**2 + 2*i*aimag(T) )/( 1 - conjg(T) ).
!
! The input GAMMA, SIGMA, C, S, U, V and RHO must satisfy the following to
! machine precision.
!
!       |GAMMA| = 1
!       |SIGMA| = 1
!     C^2 + S^2 = 1
!    |U|^2+ V^2 = 1
!         |RHO| = 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  NORMALIZE       LOGICAL array of dimension 3
!                    NORMALIZE(1): renormalize U and V
!                    NORMALIZE(2): renormalize SIGMA
!                    NORMALIZE(3): renormalize GAMMA
!
! INPUT/OUTPUT VARIABLES:
!
!  GAMMA           COMPLEX(8)
!                    unimodular phase parameter
!
!  SIGMA           COMPLEX(8)
!                    unimodular phase parameter for the left/right rotations
!
!  C               REAL(8)
!                    nonnegative cosine modulus
!
!  S               REAL(8)
!                    nonnegative sine
!
!  U               COMPLEX(8)
!                    complex component of the middle core transformation
!
!  V               REAL(8)
!                    nonnegative sine of the middle core transformation
!
!  RHO             COMPLEX(8)
!                    unimodular shift
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3_symturnover(NORMALIZE,GAMMA,SIGMA,C,S,U,V,RHO)

  implicit none

  ! input variables
  logical, intent(in) :: NORMALIZE(3)

  ! input/output variables
  real(8), intent(inout) :: C, S, V
  complex(8), intent(inout) :: GAMMA, SIGMA, U
  complex(8), intent(in) :: RHO

  ! compute variables
  complex(8) :: gs2, gammah, sigmah, z, t, Uh
  real(8) :: zz, zabs, n, xx, Ch, Sh, Vh

  ! compute t and z
  gs2 = GAMMA*SIGMA*SIGMA
  t = conjg(gs2)*U

  if ( dble(t) >= 0d0 ) then
    z = 1d0 + t
  else
    z = cmplx(V**2,2d0*aimag(t),kind=8)/(1d0 - conjg(t))
  end if

  zz = dble(z)**2 + aimag(z)**2

  ! compute n
  n = sqrt(C**2*zz + V**2)

  ! n > 0
  if ( n > 0d0 ) then

    zabs = sqrt(zz)

    Uh = conjg(RHO)*(S**2*U - C**2*gs2)
    Vh = n*S
    Ch = C*zabs/n
    Sh = V/n

    ! z /= 0
    if ( zz > 0d0 ) then
      sigmah = -RHO*SIGMA*z/zabs
    else
      sigmah = cmplx(1d0,0d0,kind=8)
    end if

  else

    Ch = 0d0
    Sh = 1d0
    Uh = conjg(RHO)*U
    Vh = 0d0
    sigmah = cmplx(1d0,0d0,kind=8)

  end if

  gammah = -conjg(RHO)*GAMMA

  ! store middle core variables
  U = Uh
  V = Vh

  ! optionally renormalize the middle core transformation
  if (NORMALIZE(1)) then
    xx = dble(U)**2 + aimag(U)**2 + V**2
    xx = 5d-1*(3d0-xx)
    U = U*xx
    V = V*xx
  end if

  ! store the left/right core variables
  C = Ch
  S = Sh

  ! store phase variables
  SIGMA = sigmah
  GAMMA = gammah

  ! optionally renormalize sigma
  if (NORMALIZE(2)) then
    xx = dble(SIGMA)**2 + aimag(SIGMA)**2
    xx = 5d-1*(3d0-xx)
    SIGMA = SIGMA*xx
  end if

  ! optionally renormalize gamma
  if (NORMALIZE(3)) then
    xx = dble(GAMMA)**2 + aimag(GAMMA)**2
    xx = 5d-1*(3d0-xx)
    GAMMA = GAMMA*xx
  end if

end subroutine z_rot3_symturnover
