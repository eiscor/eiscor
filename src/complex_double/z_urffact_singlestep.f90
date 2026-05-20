#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_urffact_singlestep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' single-shift algorithm on a
! unitary upper Hessenberg matrix stored in square-root-free symmetric factored
! form.
!
! The matrix is represented by
!
!   | nu  0 | | u1       -v1 |
!   |  0  1 | | v1  conj(u1) | | u2       -v2 |
!                              | v2  conj(u2) | | u3       -v3 | | 1   0 |
!                                               | v3  conj(u3) | | 0  u4 |
!
! but the square-root-free algorithm stores only vi^2.  Thus the arrays U and
! VV contain
!
!     U(i)  = ui,
!     VV(i) = vi^2.
!
! The input must satisfy
!
!     |U(i)|^2 + VV(i) = 1,   i = 1,...,N-1,
!                  VV(N) = 0,
!                  |U(N)| = 1,
!                    |NU| = 1.
!
! The turnover variables are
!
!     omega = gamma*sigma^2,
!     cc    = c^2,
!     ss    = s^2.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix, must be >= 2
!
!  U               COMPLEX(8) array of dimension N
!                    array of complex generators for core transformations
!
!  VV              REAL(8) array of dimension N
!                    array of squared real generators vi^2
!
!  NU              COMPLEX(8)
!                    unimodular leading phase
!
!  ITCNT           INTEGER
!                    contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_urffact_singlestep(N,U,VV,NU,ITCNT)

  implicit none

  ! input/output variables
  integer, intent(in) :: N
  integer, intent(inout) :: ITCNT
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: VV(N)
  complex(8), intent(in) :: NU

  ! compute variables
  integer :: ii
  real(8) :: vvt, cc, ss, xx
  complex(8) :: ut, omega, rho
  complex(8) :: block(2,2), t1(2,2), t2(2,2)

  ! get trailing 2 x 2 block used for the projected Wilkinson shift
  if (N.EQ.2) then
    block(1,1) =  NU*U(N-1)
    block(2,2) =  conjg(NU*U(N-1))
  else
    block(1,1) =  U(N-1)
    block(2,2) =  conjg(U(N-1))
  end if

  block(1,2) = -sqrt(VV(N-1))
  block(2,1) =  sqrt(VV(N-1))
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

  ! initialize root-free symmetric turnover parameters
  !
  ! In the standard symmetric turnover, the initial phase product is
  !
  !     gamma*sigma^2 = -rho.
  !
  ! Therefore, in the square-root-free representation,
  !
  !     omega = -rho,
  !     cc    = 1,
  !     ss    = 0.
  omega = -rho
  cc = 1d0
  ss = 0d0

  ! main chasing loop
  do ii = 0,(N-1)

    ! set ut and vvt
    ut = NU*U(ii+1)
    vvt = VV(ii+1)

    ! root-free symmetric turnover
    call z_rfr3_symturnover(omega,cc,ss,ut,vvt,rho)

    ! store ut and vvt
    if (ii > 0) then
      U(ii) = conjg(NU)*ut
      VV(ii) = vvt
    end if

  end do

end subroutine z_urffact_singlestep
