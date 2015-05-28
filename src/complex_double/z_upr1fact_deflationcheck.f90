#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_deflationcheck 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks for deflations in a unitary 
! matrix that is stored as a product of Givens rotations.
! When a deflation occurs the corresponding rotation is set diagonal.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (4*N)
!                    array of generators for givens rotations
!                    generators must be orthogonal to working precision
!
! OUTPUT VARIABLES:
!
!  ZERO            INTEGER
!                     index of the last known deflation
!                     on output contains index of newest deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_deflationcheck(N,Q,ZERO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(inout) :: Q(4*(N-1))
  integer, intent(inout) :: ZERO

  ! compute variables
  integer :: ii
  real(8), parameter :: tol = EISCOR_DBL_EPS
  real(8) :: nrm
  
  ! check for deflation
  do ii=1,(N-1)
  
    ! deflate if subdiagonal is small enough
    nrm = abs(cmplx(Q(4*(N-ii)-1),Q(4*(N-ii)),kind=8))
    if(nrm < tol)then

      ! set ZERO
      ZERO = max(0,N-ii) ! why 0?
      
      ! update Q
      call z_rot4_vec4gen(Q(4*(N-ii)-3),Q(4*(N-ii)-2),0d0,0d0 &
      ,Q(4*(N-ii)-3),Q(4*(N-ii)-2),Q(4*(N-ii)-1),Q(4*(N-ii)),nrm)

      ! exit loop 
      exit
        
    end if
    
  end do
  
end subroutine z_upr1fact_deflationcheck
