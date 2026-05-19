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
! The matrix is represented by
!
!   | nu  0 | | u1       -v1 |
!   |  0  1 | | v1  conj(u1) | | u2       -v2 |
!                              | v2  conj(u2) | | u3       -v3 | | 1   0 |
!                                               | v3  conj(u3) | | 0  u4 |
!
! The arrays U and V contain
!
!     U(i) = ui,
!     V(i) = vi.
!
! The input must satisfy
!
!     |U(i)|^2 + V(i)^2 = 1,   i = 1,...,N-1,
!                  V(N) = 0,
!                  |NU| = 1.
!
! If VEC = .TRUE., the eigenvector/similarity matrix Z is updated by the
! same similarity transformations used in the QR step.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: update eigenvectors/similarity matrix
!                    .FALSE.: do not update Z
!
!  N               INTEGER
!                    dimension of matrix, must be >= 2
!
!  U               COMPLEX(8) array of dimension N
!                    complex generators of the core transformations
!
!  V               REAL(8) array of dimension N
!                    real generators of the core transformations
!
!  NU              COMPLEX(8)
!                    leading unimodular phase
!
!  M               INTEGER
!                    leading dimension / number of rows of Z
!
!  Z               COMPLEX(8) array of dimension (M,N)
!                    if VEC = .TRUE., updated by right multiplication
!                    if VEC = .FALSE., unused
!
!  ITCNT           INTEGER
!                    iteration counter since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_usymfact_singlestep(VEC,N,U,V,NU,M,Z,ITCNT)

  implicit none

  ! input/output variables
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
    call z_rot3_symturnover(gamma,sigma,c,s,ut,vt,rho)

    ! update eigenvectors / similarity accumulator
    !
    ! The computed similarity core is
    !
    !     Q =
    !       [ c*sigma       -s              ]
    !       [ s              c*conjg(sigma)]
    !
    ! acting on columns ii+1 and ii+2.
    !
    ! The final pass ii = N-1 updates the trailing diagonal core and does not
    ! correspond to a 2-column update inside the N x N eigenvector matrix.
    if (VEC .AND. ii < N-1) then
    
      do jj = 1,M
    
        z1 = Z(jj,ii+1)
        z2 = Z(jj,ii+2)
    
        Z(jj,ii+1) = z1*(c*sigma) + z2*s
        Z(jj,ii+2) = -z1*s + z2*(c*conjg(sigma))
    
      end do
    
    else if (VEC) then
    
      ! Final diagonal similarity factor Q_N.
      !
      ! At ii = N-1, the turnover returns the final one-dimensional
      ! core transformation.  It acts on the last column of Z by sigma.
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
