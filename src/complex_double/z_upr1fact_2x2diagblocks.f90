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
!  JOB             CHARACTER
!                    'T': only upper triangular parts are returned
!                    'H': the extended hessenberg part is included
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
!  D1, D2          REAL(8) array of dimension (2*(N+1))
!                    arrays of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1, B1, C2, B2  REAL(8) array of dimension (4,3*N)
!                    arrays of generators for upper-triangular parts
!                    of the pencil
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
subroutine z_upr1fact_2x2diagblocks(JOB,ALG,N,K,P,Q,D1,C1,B1,D2,C2,B2,A,B,INFO)
  
  implicit none
  
  ! input variables
  character, intent(in) :: JOB
  character(2), intent(in) :: ALG
  integer, intent(in) :: N, K
  integer, intent(inout) :: INFO
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*(N+1)), D2(2*(N+1))
  real(8), intent(inout) :: C1(3*N), B1(3*N), C2(3*N) ,B2(3*N)
  complex(8), intent(inout) :: A(2,2), B(2,2)
  
  ! compute variables
  integer :: ind
  complex(8) :: H(3,3), T(3,2)
  
  ! initialize INFO
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
    
    ! check factorization
    call z_upr1fact_factorcheck(ALG,N,Q,D1,C1,B1,D2,C2,B2,INFO)
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
  T(2,1) = cmplx(-B1(ind+3)/C1(ind+3),0d0,kind=8)

  ! if not at top  
  if (K > 1) then
    T(1,1) = (cmplx(-B1(ind-2),B1(ind-1),kind=8)*cmplx(B1(ind+1),B1(ind+2),kind=8) &
      + T(2,1)*cmplx(C1(ind-2),C1(ind-1),kind=8)*cmplx(C1(ind+1),-C1(ind+2),kind=8))/cmplx(C1(ind),0d0,kind=8)
  end if
      
  ! second column of T
  ind = 3*K
  T(3,2) = cmplx(-B1(ind+3)/C1(ind+3),0d0)
  T(2,2) = (cmplx(-B1(ind-2),B1(ind-1),kind=8)*cmplx(B1(ind+1),B1(ind+2),kind=8) &
      + T(3,2)*cmplx(C1(ind-2),C1(ind-1),kind=8)*cmplx(C1(ind+1),-C1(ind+2),kind=8))/cmplx(C1(ind),0d0,kind=8)
  
  ! if not at top
  if (K > 1) then    
    T(1,2) = (cmplx(B1(ind-5),-B1(ind-4),kind=8)*cmplx(B1(ind),0d0,kind=8)*cmplx(B1(ind+1),B1(ind+2),kind=8) - &
      cmplx(C1(ind-5),C1(ind-4),kind=8)/cmplx(C1(ind),0d0,kind=8)* &
      (cmplx(C1(ind-2),-C1(ind-1),kind=8)*cmplx(B1(ind-2),-B1(ind-1),kind=8)*cmplx(B1(ind+1),B1(ind+2),kind=8) - &
      cmplx(C1(ind+1),-C1(ind+2),kind=8)*T(3,2)))/cmplx(C1(ind-3),0d0,kind=8)
  end if
  
  ! apply diagonal
  ind = 2*(K-1)
  T(2,:) = cmplx(D1(ind+1),D1(ind+2),kind=8)*T(2,:)
  T(3,:) = cmplx(D1(ind+3),D1(ind+4),kind=8)*T(3,:)
  
  ! if not at top
  if (K > 1) then
    T(1,:) = cmplx(D1(ind-1),D1(ind),kind=8)*T(1,:)
  end if

  ! extended hessenberg part
  if (JOB.EQ.'H') then

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
        H(:,3) = H(:,3)*cmplx(Q(ind+4),Q(ind+5),kind=8)
      else
        H(3,:) = H(3,:)*cmplx(Q(ind+4),Q(ind+5),kind=8)
      end if
    end if
    
    ! set output
    A = matmul(H(2:3,:),T)
    
  ! upper triangular part only
  else
   
    ! set output
    A = T(2:3,1:2)
    
  end if

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
    B(1,1) = cmplx(-B2(ind+3)/C2(ind+3),0d0,kind=8)
        
    ! second column of T
    ind = 3*K
    B(2,2) = cmplx(-B2(ind+3)/C2(ind+3),0d0)
    B(1,2) = (cmplx(-B2(ind-2),B2(ind-1),kind=8)*cmplx(B2(ind+1),B2(ind+2),kind=8) &
        + B(2,2)*cmplx(C2(ind-2),C2(ind-1),kind=8)*cmplx(C2(ind+1),-C2(ind+2),kind=8))/cmplx(C2(ind),0d0,kind=8)
    
    ! apply diagonal
    ind = 2*(K-1)
    B(1,:) = cmplx(D2(ind+1),D2(ind+2),kind=8)*B(1,:)
    B(2,:) = cmplx(D2(ind+3),D2(ind+4),kind=8)*B(2,:)  
    
  end if

end subroutine z_upr1fact_2x2diagblocks
