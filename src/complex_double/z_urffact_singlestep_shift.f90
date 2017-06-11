#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_urffact_singlestep_shift 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift algorithm on a
! unitary upper hessenberg matrix that is stored as a product of givens
! rotations. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER 
!                    dimension of matrix, must be >= 2
!
!  U               REAL(8) array of dimension (3*N)
!                    array of generators for givens rotations
!
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_urffact_singlestep_shift(N,U,VV,ITCNT,SHIFT)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: ITCNT
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: VV(N)
  complex(8), intent(in) :: SHIFT
  
  ! compute variables
  integer :: ii
  real(8) :: nn, zz, cc, ss, xx
  complex(8) :: z, w, rho
  complex(8) :: block(2,2), t1(2,2), t2(2,2)

  rho = SHIFT

  ! initialize
  w = -rho
  z = U(1) + w
  zz = dble(z)**2 + aimag(z)**2
  nn = VV(1) + zz
  cc = 1d0

!print*,""
  ! main chasing loop
  do ii=1,(N-1)
    cc = cc*zz/nn
    ss = VV(ii)/nn
    w = -rho*conjg(w)*((z*z)/zz)
!print*,w,cc,ss,cc+ss
    xx = dble(w)**2 + aimag(w)**2
    w = w*(3d0 - xx)/2
    z = U(ii+1) + w
    zz = dble(z)**2 + aimag(z)**2
if (zz < 1d-3) then
print*,ii,w,U(ii+1),z,zz,VV(ii+1)
end if
    nn = VV(ii+1) + cc*zz
    !U(ii) = conjg(rho)*(ss*U(ii+1) - cc*w)
    U(ii) = -conjg(rho)*(w - ss*z)
    VV(ii) = ss*nn
    xx = dble(U(ii))**2 + aimag(U(ii))**2 + VV(ii)
    U(ii) = U(ii)*(3d0 - xx)/2
    VV(ii) = VV(ii)/xx
    xx = dble(U(ii))**2 + aimag(U(ii))**2 + VV(ii)
    U(ii) = U(ii)*(3d0 - xx)/2
    VV(ii) = VV(ii)/xx
  end do
  
end subroutine z_urffact_singlestep_shift
