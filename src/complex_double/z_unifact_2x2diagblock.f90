#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_2x2diagblock
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a two by two diagonal block of a unitary upper 
! hessenberg matrix that is stored as a product of givens rotations 
! and a complex diagonal matrix. 
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
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!                    generators must be orthogonal to working precision
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
! OUTPUT VARIABLES:
!
!  H              COMPLEX(8) array of dimension (2,2)
!                   on exit contains the desired 2x2 block
!
!  INFO           INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies N is invalid
!                   INFO = -2 implies K is invalid
!                   INFO = -3 implies Q is invalid
!                   INFO = -4 implies D is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_2x2diagblock(N,K,Q,D,H,INFO)
  
  implicit none
  
  ! input variables
  integer, intent(in) :: N, K
  integer, intent(inout) :: INFO
  real(8), intent(in) :: Q(3*(N-1)), D(2*N)
  complex(8), intent(inout) :: H(2,2)
  
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
      call u_infocode_check(__FILE__,__LINE__,"K must 1 <= K <= N-1",INFO,INFO)
      return
    end if 

    ! check Q and D
    call z_unifact_factorcheck(N,Q,D,INFO)
    if (INFO.EQ.-1) then
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
      return
    end if
    if (INFO.EQ.-2) then
      INFO = -3
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,INFO)
      return
    end if
    if (INFO.EQ.-3) then
      INFO = -4
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,INFO)
      return
    end if

  end if
  
  ! set index
  strt = 3*(K-1)
  
  ! initialize H
  H(1,1) = cmplx(Q(strt+1),Q(strt+2),kind=8)
  H(2,1) = cmplx(Q(strt+3),0d0,kind=8)
  H(1,2) = -H(2,1)
  H(2,2) = conjg(H(1,1))
    
  ! apply upper rotation
  if (K > 1) then
    H(1,:) = H(1,:)*cmplx(Q(strt-2),-Q(strt-1),kind=8)
    end if
    
  ! apply lower rotation
  if (K < (N-1)) then
    H(:,2) = H(:,2)*cmplx(Q(strt+4),Q(strt+5),kind=8)
  end if
    
  ! apply diagonal
  strt = 2*(k-1)
  H(:,1) = H(:,1)*cmplx(D(strt+1),D(strt+2),kind=8)
  H(:,2) = H(:,2)*cmplx(D(strt+3),D(strt+4),kind=8)

end subroutine z_unifact_2x2diagblock
