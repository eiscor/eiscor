#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_buildbulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the first transformation for a sinlge shift
! iteration on a upr1 pencil.
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
!  K               INTEGER
!                    index for the top of the relevant submatrix
!                    1 <= K < N-1 (never bulge chase at bottom)
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
!  SHFT            COMPLEX(8) 
!                    contains the shift needed for the first transformation
!
! OUTPUT VARIABLES:
!
!  G               REAL(8) array of dimension (3)
!                    on exit contains the generators for the first
!                    transformation
!
!  INFO            INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies ALG is invalid
!                   INFO = -2 implies N is invalid
!                   INFO = -3 implies K is invalid
!                   INFO = -5 implies Q is invalid
!                   INFO = -6 implies D is invalid
!                   INFO = -7 implies R is invalid
!                   INFO = -8 implies SHFT is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_buildbulge(ALG,N,K,P,Q,D1,C1,B1,D2,C2,B2,SHFT,G,INFO)
  
  implicit none
  
  ! input variables
  character(2), intent(in) :: ALG
  integer, intent(in) :: N, K
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*(N+1)), D2(2*(N+1))
  real(8), intent(inout) :: C1(3*N), B1(3*N), C2(3*N) ,B2(3*N)
  complex(8), intent(in) :: SHFT
  real(8), intent(inout) :: G(3)
  integer, intent(inout) :: INFO
  
  ! compute variables
  real(8) :: nrm
  complex(8) :: A(2,2), B(2,2), vec1(2), vec2(2)
  
  ! initialize info
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
    
    ! check factorization
    call z_upr1fact_factorcheck(ALG,N,Q,D1,C1,B1,D2,C2,B2,INFO)
    if (INFO.EQ.-1) then
      call u_infocode_check(__FILE__,__LINE__,"ALG must be 'QR' or 'QZ'",INFO,-1)
      return
    end if
    if (INFO.EQ.-2) then
      call u_infocode_check(__FILE__,__LINE__,"N is invalid",INFO,-2)
      return
    end if
    if (INFO.EQ.-3) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,-5)
      return
    end if
    if (INFO.EQ.-4) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,-6)
      return
    end if
    if (INFO.EQ.-5) then
      call u_infocode_check(__FILE__,__LINE__,"R is invalid",INFO,-7)
      return
    end if
    
    ! check K
    if ((K < 1).OR.(K >= N-1)) then
      INFO = -3
      call u_infocode_check(__FILE__,__LINE__,"K must 1 <= K < N-1",INFO,INFO)
      return
    end if 
  
    ! check SHFT
    call z_scalar_nancheck(SHFT,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"SHFT is invalid",INFO,-8)
      return
    end if
    call z_scalar_infcheck(SHFT,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"SHFT is invalid",INFO,-8)
      return
    end if
  
  end if
  
  ! get 2x2 blocks
  call z_upr1fact_2x2diagblocks('T',ALG,N,K,P,Q,D1,C1,B1,D2,C2,B2,A,B,INFO)
      
  ! check INFO in debug mode
  if (DEBUG) then
    call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_2x2diagblocks failed",INFO,1)
    if (INFO.NE.0) then 
      return 
    end if 
  end if
  
  ! compute first columns
  ! P(K) == FALSE
  if (P(K).EQV..FALSE.) then
   
    ! first column of A
    vec1(1) = cmplx(Q(3*K-2),Q(3*K-1),kind=8)*A(1,1)
    vec1(2) = cmplx(Q(3*K),0d0,kind=8)*A(1,1)
    
    ! first column of B
    vec2(1) = B(1,1)
    vec2(2) = B(2,1)
  
  ! P(K) == TRUE
  else
    
    ! Q(K)*e1
    vec2(1) = cmplx(Q(3*K-2),-Q(3*K-1),kind=8)
    vec2(2) = cmplx(-Q(3*K),0d0,kind=8)
    
    ! back solve with A
    vec2(2) = vec2(2)/A(2,2)
    vec2(1) = (vec2(1) - A(1,2)*vec2(2))/A(1,1)
    
    ! B^-1 e1
    vec1(1) = 1d0/B(1,1)
    vec1(2) = 0d0
  
  end if
  
  ! insert shift
  vec1 = vec1 - SHFT*vec2
  
  ! compute eliminator
  call z_rot3_vec4gen(dble(vec1(1)),aimag(vec1(1)),dble(vec1(2)),aimag(vec1(2)),G(1),G(2),G(3),nrm,INFO)
      
  ! check INFO in debug mode
  if (DEBUG) then
    call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec4gen failed",INFO,1)
    if (INFO.NE.0) then 
      return 
    end if 
  end if

end subroutine z_upr1fact_buildbulge
