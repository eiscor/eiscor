#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_usymfact_singlestep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' single-shift algorithm on a
! unitary upper Hessenberg matrix stored as a product of Givens rotations.
!
! NORMALIZE controls the final renormalization inside z_rot3_symturnover:
!
!     NORMALIZE(1) = .TRUE.  : renormalize U and V
!     NORMALIZE(2) = .TRUE.  : renormalize SIGMA
!     NORMALIZE(3) = .TRUE.  : renormalize GAMMA
!
! The three flags are independent.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_usymfact_singlestep(NORMALIZE,VEC,N,U,V,NU,M,Z,ITCNT)

  implicit none

  ! input/output variables
  logical, intent(in) :: NORMALIZE(3)
  logical, intent(in) :: VEC
  integer, intent(in) :: N, M
  integer, intent(inout) :: ITCNT
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: V(N)
  complex(8), intent(in) :: NU
  complex(8), intent(inout) :: Z(M,N)

  ! compute variables
  integer :: ii, jj
  real(8) :: vt, c, s, xx
  complex(8) :: ut, rho
  complex(8) :: gamma, sigma
  complex(8) :: z1, z2
  complex(8) :: block(2,2), t1(2,2), t2(2,2)

  ! get trailing 2 x 2 block used for the projected Wilkinson shift
  if (N.EQ.2) then
    block(1,1) =  NU*U(N-1)
    block(2,2) =  conjg(NU*U(N-1))
  else
    block(1,1) =  U(N-1)
    block(2,2) =  conjg(U(N-1))
  end if

  block(1,2) = -V(N-1)
  block(2,1) =  V(N-1)
  block(:,2) =  block(:,2)*U(N)

  if (N > 2) then
    xx = abs(U(N-2))
    if (xx.GT.0d0) then
      block(1,:) = conjg(U(N-2))*block(1,:)/xx
    end if
  end if

  ! compute eigenvalues of the trailing 2 x 2 block
  t1 = block
  call z_2x2array_eig(.FALSE.,t1,t1,t2,t2)

  ! choose Wilkinson shift
  if (abs(block(2,2)-t1(1,1)) < abs(block(2,2)-t1(2,2))) then
    rho = t1(1,1)
  else
    rho = t1(2,2)
  end if

  ! project shift onto the unit circle
  xx = abs(rho)

  if (xx == 0d0) then
    call random_number(xx)
    rho = cmplx(cos(xx),sin(xx),kind=8)
  else
    rho = rho/xx
  end if

  ! initialize symmetric turnover parameters
  !
  ! Older implementations used a single phase W = -rho.
  ! In the current standard representation,
  !
  !     W = gamma*sigma^2.
  !
  ! We choose sigma = 1 and gamma = -rho.
  gamma = -rho
  sigma = cmplx(1d0,0d0,kind=8)
  c = 1d0
  s = 0d0

  ! main chasing loop
  do ii = 0,(N-1)

    ! set ut and vt
    ut = NU*U(ii+1)
    vt = V(ii+1)

    ! symmetric turnover
    call z_rot3_symturnover(NORMALIZE,gamma,sigma,c,s,ut,vt,rho)

    ! update eigenvectors / similarity accumulator
    !
    ! The computed similarity core is
    !
    !     Q =
    !       [ c*sigma       -s              ]
    !       [ s              c*conjg(sigma)]
    !
    ! acting on columns ii+1 and ii+2.
    if (VEC .AND. ii < N-1) then

      do jj = 1,M

        z1 = Z(jj,ii+1)
        z2 = Z(jj,ii+2)

        Z(jj,ii+1) = z1*(c*sigma) + z2*s
        Z(jj,ii+2) = -z1*s + z2*(c*conjg(sigma))

      end do

    else if (VEC) then

      ! Final diagonal similarity factor Q_N.
      do jj = 1,M
        Z(jj,N) = Z(jj,N)*sigma
      end do

    end if

    ! store ut and vt
    if (ii > 0) then
      U(ii) = conjg(NU)*ut
      V(ii) = vt
    end if

  end do

  ! increment iteration counter
  ITCNT = ITCNT + 1

end subroutine z_usymfact_singlestep
