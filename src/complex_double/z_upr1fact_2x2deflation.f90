#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_2x2deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine deflates a 2x2 block in a upr1 pencil
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  ALG             CHARACTER(2)
!                    'QR': second triangular factor is assumed to be identity
!                    'QZ': second triangular factor is assumed nonzero
!
!  COMPZ           CHARACTER
!                    'N': no schurvectors
!                    'I': schurvectors, initializing V and W to the identity
!                    'V': schurvectors, assume V and W already initialized
!
!  N               INTEGER
!                    dimension of matrix
!
!  K               INTEGER
!                    index of block to be deflated
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
!  V              COMPLEX(8) array of dimension (N,N)
!                   right schur vectors
!                   if COMPZ = 'N' unused
!                   if COMPZ = 'I' stores right schurvectors in V 
!                   if COMPZ = 'V' update V to store right schurvectors 
!
!  W              COMPLEX(8) array of dimension (N,N)
!                   left schur vectors
!                   if COMPZ = 'N' unused
!                   if COMPZ = 'I' stores left schurvectors in W 
!                   if COMPZ = 'V' update W to store left schurvectors 
!
!  INFO           INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies ALG, N, Q, D or R is invalid
!                   INFO = -2 implies K is invalid
!                   INFO = -3 implies COMPZ is invalid
!                   INFO = -8 implies V is invalid
!                   INFO = -9 implies W is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_2x2deflation(ALG,COMPZ,N,K,P,Q,D,R,V,W,INFO)

  implicit none
  
  ! input variables
  character(2), intent(in) :: ALG
  character, intent(in) :: COMPZ
  integer, intent(in) :: N, K
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D(2,2*(N+1)), R(4,3*N)
  complex(8), intent(inout) :: V(N,N), W(N,N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  complex(8) :: A(2,2), B(2,2), Vt(2,2), Wt(2,2)
  
  ! initialize info
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
  
    ! check factorization
    call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"ALG, N, Q, D or R is invalid",INFO,-1)
      return
    end if
    
    ! check K
    if ((K < 1).OR.(K > N-1)) then
      INFO = -2
      call u_infocode_check(__FILE__,__LINE__,"K must 1 <= K <= N-1",INFO,INFO)
      return
    end if 
  
    ! check COMPZ
    if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
      INFO = -3
      call u_infocode_check(__FILE__,__LINE__,"COMPZ must be 'N', 'I' or 'V'",INFO,INFO)
      return
    end if
    
    ! check V
    if (COMPZ.EQ.'V') then
      call z_2Darray_check(N,N,V,INFO)
      if (INFO.NE.0) then
        call u_infocode_check(__FILE__,__LINE__,"V is invalid",INFO,-8)
        return
      end if
    end if   
    
    ! check W
    if ((ALG.EQ.'QZ').AND.(COMPZ.EQ.'V')) then
      call z_2Darray_check(N,N,W,INFO)
      if (INFO.NE.0) then
        call u_infocode_check(__FILE__,__LINE__,"W is invalid",INFO,-9)
      end if
    end if
    
  end if
  
  ! compute 2x2 blocks
  call z_upr1fact_2x2diagblocks(N,K,ALG,P,Q,D,R,A,B,INFO)
    
  ! check INFO in debug mode
  if (DEBUG) then
    call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_2x2diagblocks failed",INFO,INFO)
    if (INFO.NE.0) then 
      return 
    end if 
  end if
      
  ! compute generalized Schur decomposition
  call z_2x2array_geneig(A,B,Wt,Vt,INFO)
    
  ! check INFO in debug mode
  if (DEBUG) then
    call u_infocode_check(__FILE__,__LINE__,"z_2x2array_geneig failed",INFO,INFO)
    if (INFO.NE.0) then 
      return 
    end if 
  end if



end subroutine z_upr1fact_2x2deflation
