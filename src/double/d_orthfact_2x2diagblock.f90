#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthfact_2x2diagblock 
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
!                    INFO = -3 implies Q is invalid
!                    INFO = -4 implies D is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthfact_2x2diagblock(N,K,Q,D,H,INFO)
  
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
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
      return
    end if
  
    ! check K
    if ((K < 1).OR.(K > N-1)) then
      INFO = -2
      call u_infocode_check(__FILE__,__LINE__,"K must be 1 <= K <= N-1",INFO,INFO)
      return
    end if

    ! check Q and D
    call d_orthfact_factorcheck(N,Q,D,INFO)
    if (INFO.EQ.-1) then
      call u_infocode_check(__FILE__,__LINE__,"N is invalid",INFO,-1)
      return
    end if
    if (INFO.EQ.-2) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,-3)
      return
    end if
    if (INFO.EQ.-3) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,-4)
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

end subroutine d_orthfact_2x2diagblock
