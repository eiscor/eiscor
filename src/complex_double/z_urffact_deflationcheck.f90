#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_urffact_deflationcheck 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks for deflations in a unitary upper hessenberg matrix
! that is stored as a product of givens rotations. When a deflation occurs the
! corresponding rotation is set to the identity matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  U               complex(8) array of dimension (N)
!                    array of generators for givens rotations
!
!  VV              REAL(8) array of dimension (N)
!                    array of generators for complex diagonal matrix
!
! OUTPUT VARIABLES:
!
!  ZERO            INTEGER
!                     on output contains index of newest deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_urffact_deflationcheck(N,U,VV,ZERO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: ZERO
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: VV(N)

  ! compute variables
  integer :: ii
  real(8), parameter :: tol = (EISCOR_DBL_EPS)**2
  real(8) :: xx

  ! intialize ZERO
  ZERO = 0
  
  ! check for deflation
  do ii=1,N
  
    ! deflate if subdiagonal is small enough
    if (VV(N+1-ii) < tol) then
        
      ! set ZERO
      ZERO = N+1-ii

      ! set rotation to diagonal
      VV(ZERO) = 0d0
        
      ! renormalize U
      xx = dble(U(ZERO))**2 + aimag(U(ZERO))**2
      U(ZERO) = 5d-1*U(ZERO)*(3d0-xx)

    end if

  end do
  
end subroutine z_urffact_deflationcheck
