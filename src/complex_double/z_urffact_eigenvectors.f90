#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_urffact_eigenvectors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine diagonalizes a unitary upper hessenberg matrix that is stored as
! a product of N Givens rotations, without computing square roots.
!
! | u1       -v1 |
! | v1  conj(u1) | | u2       -v2 | 
!                  | v2  conj(u2) | | u3       -v3 | | 1   0 |
!                                   | v3  conj(u3) | | 0  u4 | 
!                                                                     
! The square root free algorithm only requires the storage of the vi^2,
! so the arrays U and VV contain the following:
!
!  U(i) = ui
! VV(i) = vi^2
!
! The input must satisfy the following:
!
!  |U(i)|^2 + VV(i) = 1
!             VV(N) = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  U               COMPLEX(8) array of dimension N
!                    array of complex generators for Givens rotations
!                    on output contains eigenvalues
!
!  VV              REAL(8) array of dimension N
!                    array of real generators (squared) for Givens rotations
!
! OUTPUT VARIABLES:
!
!  ITS             INTEGER array of dimension N-1
!                    contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO = 1 implies no convergence
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies N is invalid
!                    INFO = -2 implies U or VV is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_urffact_eigenvectors(N,U,VV,E,M,Z,WORK)

  implicit none
  
  ! input variables
  integer, intent(in) :: N,M
  complex(8), intent(in) :: U(N), E(N)
  real(8), intent(in) :: VV(N)
  complex(8), intent(inout) :: Z(M,N), WORK(N)
  
  ! compute variables
  integer :: ii, sgn 
  
  ! initialize Z
  Z = cmplx(0d0,0d0,kind=8)
  Z(1,:) = cmplx(1d0,0d0,kind=8)
  
  ! initialize WORK
  WORK = cmplx(1d0,0d0,kind=8)
 
  ! initialize sgn
  sgn = -1 
 
  ! fill Z
  do ii = 1,M-1

    ! row ii+1 of Z
    Z(ii+1,:) = (E*Z(ii,:) + dble(sgn)*U(ii)*WORK)/sqrt(VV(ii))
    
    ! update WORK
    WORK = (dble(sgn)*conjg(U(ii))*E*Z(ii,:) + WORK)/sqrt(VV(ii))
    
    ! update sgn
    sgn = -sgn 

  end do

  ! compute column norms
  WORK = abs(Z(1,:))**2
  do ii = 1,M-1

    ! row ii+1 of Z
    WORK = WORK + abs(Z(ii+1,:))**2

  end do
  WORK = sqrt(WORK)

  ! scale columns
  do ii = 1,n
    Z(:,ii) = conjg(Z(:,ii))/WORK(ii)
  end do

end subroutine z_urffact_eigenvectors
