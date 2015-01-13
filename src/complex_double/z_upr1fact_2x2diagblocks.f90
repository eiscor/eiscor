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
!                   INFO = -1 implies ALG, N, Q, D or R is invalid
!                   INFO = -2 implies K is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_2x2diagblocks(N,K,ALG,P,Q,D,R,A,B,INFO)
  
  implicit none
  
  ! input variables
  character(2), intent(in) :: ALG
  integer, intent(in) :: N, K
  integer, intent(inout) :: INFO
  logical, intent(in) :: P(N-2)
  real(8), intent(in) :: Q(3*(N-1)), D(2,2*(N+1)), R(4,3*N)
  complex(8), intent(inout) :: A(2,2), B(2,2)
  
  ! compute variables
  integer :: ind
  complex(8) :: H(3,3), T(3,2)
  
  ! initialize INFO
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
    
    ! check factorization
    !call z_upr1fact_factorcheck(ALG,N,P,Q,D,R,INFO) ! z_upr1fact_factorcheck does not know P
     call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"ALG, N, Q, D or R is invalid",INFO,INFO)
      INFO = -1
      return
    end if
    
    ! check K
    if ((K < 1).OR.(K > N-1)) then
      INFO = -2
      call u_infocode_check(__FILE__,__LINE__,"K must 1 <= K <= N-1",INFO,INFO)
      return
    end if 

  end if
  
  ! compute A
  ! set index
  ind = 3*(K-1)
  
  ! initialize 
  T = cmplx(0d0,0d0,kind=8)

  ! first column of T
  ind = 3*(K-1)
  T(2,1) = cmplx(-R(2,ind+3)/R(1,ind+3),0d0,kind=8)

  ! if not at top  
  if (K > 1) then
    T(1,1) = (cmplx(-R(2,ind-2),R(2,ind-1),kind=8)*cmplx(R(2,ind+1),R(2,ind+2),kind=8) &
      + T(2,1)*cmplx(R(1,ind-2),R(1,ind-1),kind=8)*cmplx(R(1,ind+1),-R(1,ind+2),kind=8))/cmplx(R(1,ind),0d0,kind=8)
  end if
      
  ! second column of T
  ind = 3*K
  T(3,2) = cmplx(-R(2,ind+3)/R(1,ind+3),0d0)
  T(2,2) = (cmplx(-R(2,ind-2),R(2,ind-1),kind=8)*cmplx(R(2,ind+1),R(2,ind+2),kind=8) &
      + T(3,2)*cmplx(R(1,ind-2),R(1,ind-1),kind=8)*cmplx(R(1,ind+1),-R(1,ind+2),kind=8))/cmplx(R(1,ind),0d0,kind=8)
  
  ! if not at top
  if (K > 1) then    
    T(1,2) = (cmplx(R(2,ind-5),-R(2,ind-4),kind=8)*cmplx(R(2,ind),0d0,kind=8)*cmplx(R(2,ind+1),R(2,ind+2),kind=8) - &
      cmplx(R(1,ind-5),R(1,ind-4),kind=8)/cmplx(R(1,ind),0d0,kind=8)* &
      (cmplx(R(1,ind-2),-R(1,ind-1),kind=8)*cmplx(R(2,ind-2),-R(2,ind-1),kind=8)*cmplx(R(2,ind+1),R(2,ind+2),kind=8) - &
      cmplx(R(1,ind+1),-R(1,ind+2),kind=8)*T(3,2)))/cmplx(R(1,ind-3),0d0,kind=8)
  end if
  
  ! apply diagonal
  ind = 2*(K-1)
  T(2,:) = cmplx(D(1,ind+1),D(1,ind+2),kind=8)*T(2,:)
  T(3,:) = cmplx(D(1,ind+3),D(1,ind+4),kind=8)*T(3,:)
  
  ! if not at top
  if (K > 1) then
    T(1,:) = cmplx(D(1,ind-1),D(1,ind),kind=8)*T(1,:)
  end if

  ! build local Q
  H = cmplx(0d0,0d0,kind=8)
  H(1,1) = cmplx(1d0,0d0,kind=8)
  
  ind = 3*(K-1)  
  H(2,2) = cmplx(Q(ind+1),Q(ind+2),kind=8)
  H(3,2) = cmplx(Q(ind+3),0d0,kind=8)
  H(2,3) = cmplx(-Q(ind+3),0d0,kind=8)
  H(3,3) = cmplx(Q(ind+1),-Q(ind+2),kind=8)
  
  ! if not at top
  if (K > 1) then
    A(1,1) = cmplx(Q(ind-2),Q(ind-1),kind=8)
    A(2,1) = cmplx(Q(ind),0d0,kind=8)
    A(1,2) = cmplx(-Q(ind),0d0,kind=8)
    A(2,2) = cmplx(Q(ind-2),-Q(ind-1),kind=8)
    
    if (P(K-1).EQV..FALSE.) then  
      H(1:2,:) = matmul(A,H(1:2,:))
    else
      H(:,1:2) = matmul(H(:,1:2),A)
    end if
  end if
  
  ! if not at bottom
  if (K < (N-1)) then
    if (P(K).EQV..FALSE.) then  
      H(:,3) = H(:,3)*cmplx(Q(ind-2),Q(ind-1),kind=8)
    else
      H(3,:) = H(3,:)*cmplx(Q(ind-2),Q(ind-1),kind=8)
    end if
  end if
  
  ! set output
  A = matmul(H(2:3,:),T)

  ! compute B
  ! set to I if ALG == QR
  if (ALG.EQ."QR") then
  
    B = cmplx(0d0,0d0,kind=8)
    B(1,1) = cmplx(1d0,0d0,kind=8)
    B(2,2) = cmplx(1d0,0d0,kind=8)
  
  ! compute upper-triangular part otherwise
  else
  
    ! set index
    ind = 3*(K-1)
    
    ! initialize 
    B = cmplx(0d0,0d0,kind=8)

    ! first column of T
    ind = 3*(K-1)
    B(1,1) = cmplx(-R(4,ind+3)/R(3,ind+3),0d0,kind=8)
        
    ! second column of T
    ind = 3*K
    B(2,2) = cmplx(-R(4,ind+3)/R(3,ind+3),0d0)
    B(1,2) = (cmplx(-R(4,ind-2),R(4,ind-1),kind=8)*cmplx(R(4,ind+1),R(4,ind+2),kind=8) &
        + B(2,2)*cmplx(R(3,ind-2),R(3,ind-1),kind=8)*cmplx(R(3,ind+1),-R(3,ind+2),kind=8))/cmplx(R(3,ind),0d0,kind=8)
    
    ! apply diagonal
    ind = 2*(K-1)
    B(1,:) = cmplx(D(2,ind+1),D(2,ind+2),kind=8)*B(1,:)
    B(2,:) = cmplx(D(2,ind+3),D(2,ind+4),kind=8)*B(2,:)  
    
  end if

end subroutine z_upr1fact_2x2diagblocks
