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
! The input OMEGA, CC, SS, U, VV and RHO must satisfy the following to
! machine precision:
!
!        |OMEGA| = 1
!        CC + SS = 1
!     |U|^2 + VV = 1
!           |RHO| = 1
!
! On output the same variables contain the post-turnover root-free data.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
subroutine z_rfr3_symturnover(OMEGA,CC,SS,U,VV,RHO)

  implicit none

  ! input/output variables
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

  ! ensure normality of the middle core transformation
  !
  ! If |Uh|^2 + VVh = N, then use
  !
  !     Uh <- Uh*(3-N)/2,
  !     VVh <- VVh*(2-N).
  xx = dble(Uh)**2 + aimag(Uh)**2 + VVh
  U  = Uh*(5d-1*(3d0-xx))
  VV = VVh*(2d0-xx)

  ! store squared left/right core variables
  CC = CCh
  SS = SSh

  ! ensure unimodularity of omega
  xx = dble(OMEGAh)**2 + aimag(OMEGAh)**2
  OMEGA = OMEGAh*(5d-1*(3d0-xx))
!  OMEGA = OMEGAh

end subroutine z_rfr3_symturnover
