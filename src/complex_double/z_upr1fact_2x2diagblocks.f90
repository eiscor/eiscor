#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_2x2diagblocks
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a set of two by two diagonal blocks of a 
! unitary plus rank one matrix pencil stored in factored form.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  ALG             CHARACTER(2)
!                    'QR': second triangular factor is assumed to be identity
!                    'QZ': second triangular factor is assumed nonzero
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D               REAL(8) array of dimension (2,2*(N+1))
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!                    D1 = D(1,:)
!                    D2 = D(2,:)
!
!  R               REAL(8) array of dimension (4,3*N)
!                    array of generators for upper-triangular parts
!                    of the pencil
!                    C1 = R(1,:)
!                    B1 = R(2,:)
!                    C2 = R(3,:)
!                    B2 = R(4,:)
!
! OUTPUT VARIABLES:
!
!  A              COMPLEX(8) array of dimension (2,2)
!                   on exit contains the desired 2x2 block from
!                   the extended hessenberg matrix 
!
!  B              COMPLEX(8) array of dimension (2,2)
!                   on exit contains the desired 2x2 block from
!                   the upper-triangular matrix
!
!  INFO           INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies N is invalid
!                   INFO = -2 implies K is invalid
!                   INFO = -3 implies ALG is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_2x2diagblocks(N,K,ALG,P,Q,D,R,A,B,INFO)
  
  implicit none
  
  ! input variables
  character(2), intent(in) :: ALG
  integer, intent(in) :: N, K
  logical, intent(in) :: P(N-2)
  real(8), intent(in) :: Q(3*(N-1)), D(2,2*(N+1)), R(4,3*N)
  complex(8), intent(inout) :: A(2,2), B(2,2)
  
  ! compute variables
  integer :: ind
  
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
    
    ! check ALG
    if ((ALG.NE.'QR').AND.(ALG.NE.'QZ')) then
      INFO = -3
      call u_infocode_check(__FILE__,__LINE__,"ALG must be 'QR' or 'QZ'",INFO,INFO)
      return
    end if

  end if
  
  ! set index
  ind = 3*(K-1)
  
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

end subroutine z_upr1fact_2x2diagblocks
