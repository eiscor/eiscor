#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DOFTDB (Double Orthogonal hessenberg Factored Two by two Diagonal Block)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a two by two diagonal block of an orthogonal
! upper hessenberg matrix that is stored as a product of givens 
! rotations and a diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  K               INTEGER
!                    index of the diagonal block
!
!  Q               REAL(8) array of dimension (2*(N-1))
!                    array of generators for givens rotations
!                    generators must be orthogonal to working precision
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
! OUTPUT VARIABLES:
!
!  H               COMPLEX(8) array of dimension (2,2)
!                    on exit contains the desired 2x2 block
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies N is invalid
!                    INFO = -2 implies K is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DOFTDB(N,K,Q,D,H,INFO)
  
  implicit none
  
  ! input variables
  integer, intent(in) :: N, K
  integer, intent(inout) :: INFO
  real(8), intent(in) :: Q(2*(N-1)), D(2*N)
  real(8), intent(inout) :: H(2,2)
  
  ! compute variables
  integer :: strt
  
  ! initialize INFO
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
  
    ! check N
    if (N < 2) then
      INFO = -1
      call UARERR(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
      return
    end if
  
    ! check K
    if ((K < 1).OR.(K > N-1)) then
      INFO = -2
      call UARERR(__FILE__,__LINE__,"K must be 1 <= K <= N-1",INFO,INFO)
      return
    end if
    
  end if
  
  ! set index
  strt = 2*(K-1)
  
  ! initialize H
  H(1,1) = Q(strt+1)
  H(2,1) = Q(strt+2)
  H(1,2) = -H(2,1)
  H(2,2) = H(1,1)

  ! apply lower rotation
  if(K < (N-1))then
      H(:,2) = H(:,2)*Q(strt+3)
  end if
  
  ! apply upper rotation
  if (K > 1) then
    H(1,:) = H(1,:)*Q(strt-1)
  end if
    
  ! apply diagonal
  H(:,1) = H(:,1)*D(strt+1)
  H(:,2) = H(:,2)*D(strt+3)

end subroutine DOFTDB
