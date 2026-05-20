#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rfr3_symturnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the square-root-free symmetric turnover corresponding
! to z_rot3_symturnover.
!
! The standard variables
!
!     gamma, sigma, c, s, u, v, rho
!
! are replaced by the root-free variables
!
!     omega = gamma*sigma^2,
!     cc    = c^2,
!     ss    = s^2,
!     vv    = v^2.
!
! The scalar Z used below is the scaled quantity
!
!     Z = 1 + T,     T = conjg(OMEGA)*U,
!
! rather than the unscaled sum OMEGA + U. Since OMEGA is unimodular,
! both quantities have the same modulus.
!
! If real(T) < 0, Z is computed with the cancellation-avoiding formula
!
!     Z = ( VV + 2*i*aimag(T) )/( 1 - conjg(T) ).
!
! NORMALIZE controls the final renormalization:
!
!     NORMALIZE(1) = .TRUE.  : renormalize U and VV
!     NORMALIZE(2) = .TRUE.  : renormalize OMEGA
!
! The two flags are independent.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  NORMALIZE       LOGICAL array of dimension 2
!                    NORMALIZE(1): renormalize U and VV
!                    NORMALIZE(2): renormalize OMEGA
!
! INPUT/OUTPUT VARIABLES:
!
!  OMEGA           COMPLEX(8)
!                    root-free phase variable gamma*sigma^2
!
!  CC, SS          REAL(8)
!                    squared cosine and sine variables, c^2 and s^2
!
!  U               COMPLEX(8)
!                    complex component of the middle core transformation
!
!  VV              REAL(8)
!                    squared sine variable, v^2
!
!  RHO             COMPLEX(8)
!                    unimodular shift
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rfr3_symturnover(NORMALIZE,OMEGA,CC,SS,U,VV,RHO)

  implicit none

  ! input/output variables
  logical, intent(in) :: NORMALIZE(2)
  real(8), intent(inout) :: CC, SS, VV
  complex(8), intent(inout) :: OMEGA, U
  complex(8), intent(in) :: RHO

  ! compute variables
  complex(8) :: z, t, Uh, OMEGAh
  real(8) :: zz, m, xx, CCh, SSh, VVh

  ! compute t and scaled z
  t = conjg(OMEGA)*U

  if ( dble(t) >= 0d0 ) then
    z = 1d0 + t
  else
    z = cmplx(VV,2d0*aimag(t),kind=8)/(1d0 - conjg(t))
  end if

  zz = dble(z)**2 + aimag(z)**2

  ! compute m = n^2
  m = CC*zz + VV

  ! m > 0
  if ( m > 0d0 ) then

    CCh = CC*zz/m
    SSh = VV/m
    Uh  = conjg(RHO)*(SS*U - CC*OMEGA)
    VVh = SS*m

    ! z /= 0
    if ( zz > 0d0 ) then
      OMEGAh = -RHO*OMEGA*z**2/zz
    else
      OMEGAh = cmplx(1d0,0d0,kind=8)
    end if

  else

    CCh = 0d0
    SSh = 1d0
    Uh  = conjg(RHO)*U
    VVh = 0d0
    OMEGAh = cmplx(1d0,0d0,kind=8)

  end if

  ! store squared left/right core variables
  CC = CCh
  SS = SSh

  ! store middle core variables
  U  = Uh
  VV = VVh

  ! store phase variable
  OMEGA = OMEGAh

  ! optionally renormalize U and VV
  !
  ! If |U|^2 + VV = N, use
  !
  !     U  <- U*(3-N)/2,
  !     VV <- VV*(2-N).
  if (NORMALIZE(1)) then
    xx = dble(U)**2 + aimag(U)**2 + VV
    U  = U*(5d-1*(3d0-xx))
    VV = VV*(2d0-xx)
  end if

  ! optionally renormalize OMEGA
  !
  ! If |OMEGA|^2 = N, use
  !
  !     OMEGA <- OMEGA*(3-N)/2.
  if (NORMALIZE(2)) then
    xx = dble(OMEGA)**2 + aimag(OMEGA)**2
    OMEGA = OMEGA*(5d-1*(3d0-xx))
  end if

end subroutine z_rfr3_symturnover
