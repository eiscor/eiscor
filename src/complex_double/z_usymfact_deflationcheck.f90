#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_usymfact_deflationcheck 
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
! so the arrays U and V contain the following:
!
!  U(i) = ui
!  V(i) = vi
!
! The input must satisfy the following:
!
!  |U(i)|^2 + V(i)^2 = 1
!               V(N) = 0
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
!
!  V               REAL(8) array of dimension N
!                    array of real generators for Givens rotations
!
! OUTPUT VARIABLES:
!
!  ZERO            INTEGER
!                    largest index such that V(i) < tol
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_usymfact_deflationcheck(N,U,V,ZERO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: ZERO
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: V(N)

  ! compute variables
  integer :: ii
  real(8), parameter :: tol = EISCOR_DBL_EPS
  real(8) :: xx

  ! intialize ZERO
  ZERO = 0
  
  ! check for deflation
  do ii=1,N
  
    ! deflate if subdiagonal is small enough
    if (V(N+1-ii) < tol) then
        
      ! set ZERO
      ZERO = N+1-ii

      ! set rotation to diagonal
      V(ZERO) = 0d0
        
      ! renormalize U
      xx = dble(U(ZERO))**2 + aimag(U(ZERO))**2
      U(ZERO) = 5d-1*U(ZERO)*(3d0-xx)

    end if

  end do
  
end subroutine z_usymfact_deflationcheck
