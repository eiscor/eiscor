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
!                   INFO = 1 implies subroutine failed
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies ALG is invalid
!                   INFO = -2 implies COMPZ is invalid
!                   INFO = -3 implies N is invalid
!                   INFO = -4 implies K is invalid
!                   INFO = -6 implies Q is invalid
!                   INFO = -7 implies D1 is invalid
!                   INFO = -8 implies C1 is invalid
!                   INFO = -9 implies B1 is invalid
!                   INFO = -10 implies D2 is invalid
!                   INFO = -11 implies C2 is invalid
!                   INFO = -12 implies B2 is invalid
!                   INFO = -13 implies V is invalid
!                   INFO = -14 implies W is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_2x2deflation(ALG,COMPZ,N,K,P,Q,D1,C1,B1,D2,C2,B2,V,W,INFO)

  implicit none
  
  ! input variables
  character(2), intent(in) :: ALG
  character, intent(in) :: COMPZ
  integer, intent(in) :: N, K
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*(N+1)), D2(2*(N+1))
  real(8), intent(inout) :: C1(3*N), B1(3*N), C2(3*N) ,B2(3*N)
  complex(8), intent(inout) :: V(N,N), W(N,N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  real(8) :: nrm, G1(3), G2(3), G3(3)
  complex(8) :: A(2,2), B(2,2), Vt(2,2), Wt(2,2)
  
  ! initialize info
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
  
    ! check factorization
    call z_upr1fact_factorcheck(ALG,N,Q,D1,C1,B1,D2,C2,B2,INFO)
    if (INFO.EQ.-1) then
      call u_infocode_check(__FILE__,__LINE__,"ALG must be 'QR' or 'QZ'",INFO,-1)
      return
    else if (INFO.EQ.-2) then
      call u_infocode_check(__FILE__,__LINE__,"N is invalid",INFO,-3)
      return
    else if (INFO.EQ.-3) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,-6)
      return
    else if (INFO.EQ.-4) then
      call u_infocode_check(__FILE__,__LINE__,"D1 is invalid",INFO,-7)
      return
    else if (INFO.EQ.-5) then
      call u_infocode_check(__FILE__,__LINE__,"C1 is invalid",INFO,-8)
      return
    else if (INFO.EQ.-6) then
      call u_infocode_check(__FILE__,__LINE__,"B1 is invalid",INFO,-9)
      return
    else if (INFO.EQ.-7) then
      call u_infocode_check(__FILE__,__LINE__,"D2 is invalid",INFO,-10)
      return
    else if (INFO.EQ.-8) then
      call u_infocode_check(__FILE__,__LINE__,"C2 is invalid",INFO,-11)
      return
    else if (INFO.EQ.-9) then
      call u_infocode_check(__FILE__,__LINE__,"B2 is invalid",INFO,-12)
      return
    else if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_factorcheck failed",INFO,1)
      return
    end if   
    
    ! check K
    if ((K < 1).OR.(K > N-1)) then
      INFO = -4
      call u_infocode_check(__FILE__,__LINE__,"K must 1 <= K <= N-1",INFO,INFO)
      return
    end if 
  
    ! check COMPZ
    if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
      INFO = -2
      call u_infocode_check(__FILE__,__LINE__,"COMPZ must be 'N', 'I' or 'V'",INFO,INFO)
      return
    end if
    
    ! check V
    if (COMPZ.EQ.'V') then
      call z_2Darray_check(N,N,V,INFO)
      if (INFO.NE.0) then
        call u_infocode_check(__FILE__,__LINE__,"V is invalid",INFO,-13)
        return
      else if (INFO.NE.0) then
        call u_infocode_check(__FILE__,__LINE__,"z_2Darray_check failed",INFO,1)
        return
      end if 
    end if   
    
    ! check W
    if ((ALG.EQ.'QZ').AND.(COMPZ.EQ.'V')) then
      call z_2Darray_check(N,N,W,INFO)
      if (INFO.NE.0) then
        call u_infocode_check(__FILE__,__LINE__,"W is invalid",INFO,-14)
        return
      else if (INFO.NE.0) then
        call u_infocode_check(__FILE__,__LINE__,"z_2Darray_check failed",INFO,1)
        return
      end if 
    end if
    
  end if
  
  ! compute 2x2 blocks
  call z_upr1fact_2x2diagblocks('H',ALG,N,K,P,Q,D1,C1,B1,D2,C2,B2,A,B,INFO)
    
  ! check INFO in debug mode
  if ((DEBUG).AND.(INFO.NE.0)) then
    call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_2x2diagblocks failed",INFO,1)
    return 
  end if
      
  ! compute standard Schur decomposition
  if (ALG.EQ.'QR') then
  
    ! standard schur decomposition
    call z_2x2array_geneig('S',A,B,Wt,Vt,INFO)
      
    ! check INFO in debug mode
    if ((DEBUG).AND.(INFO.NE.0)) then
      call u_infocode_check(__FILE__,__LINE__,"z_2x2array_geneig failed",INFO,1)
      return 
    end if
    
    ! replace Vt with rotation G1
    call z_rot3_vec4gen(dble(Vt(1,1)),aimag(Vt(1,1)),dble(Vt(2,1)),aimag(Vt(2,1)),G1(1),G1(2),G1(3),nrm,INFO)
    
    ! check INFO in debug mode
    if ((DEBUG).AND.(INFO.NE.0)) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec4gen",INFO,1)
      return 
    end if
    
    ! pass G1 through triangular factor
    G3 = G1
    call z_upr1fact_rot3throughtri('R2L',N,K,D1,C1,B1,G3,INFO)
    
    ! check INFO in debug mode
    if ((DEBUG).AND.(INFO.NE.0)) then
      call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_rot3throughtri",INFO,1)
      return 
    end if
    
    ! similarity transform of Q
    A(1,1) = cmplx(Q(3*K-2),Q(3*K-1),kind=8)
    A(2,1) = cmplx(Q(3*K),0d0,kind=8)
    A(1,2) = -A(2,1)
    A(2,2) = conjg(A(1,1))
    
    B(1,1) = cmplx(G3(1),G3(2),kind=8)
    B(2,1) = cmplx(G3(3),0d0,kind=8)
    B(1,2) = -B(2,1)
    B(2,2) = conjg(B(1,1))
    
    A = matmul(A,B)
    
    ! update V
    if (COMPZ.NE.'N') then
    
      B(1,1) = cmplx(G1(1),G1(2),kind=8)
      B(2,1) = cmplx(G1(3),0d0,kind=8)
      B(1,2) = -B(2,1)
      B(2,2) = conjg(B(1,1))
      
      V(:,K:(K+1)) = matmul(V(:,K:(K+1)),B)
    
    end if

    B(1,1) = cmplx(G1(1),-G1(2),kind=8)
    B(2,1) = cmplx(G1(3),0d0,kind=8)
    B(1,2) = -B(2,1)
    B(2,2) = conjg(B(1,1))    

    A = matmul(B,A)
    
    Q(3*K-2) = 1d0
    Q(3*K-1) = 0d0
    Q(3*K) = 0d0
        
    ! deflate into diagonal D
    nrm = dble(A(1,1))*D1(2*K-1) - aimag(A(1,1))*D1(2*K)
    D1(2*K) = dble(A(1,1))*D1(2*K) + aimag(A(1,1))*D1(2*K-1)
    D1(2*K-1) = nrm
    nrm = sqrt(D1(2*K-1)**2 + D1(2*K)**2)
    D1(2*K-1) = D1(2*K-1)/nrm
    D1(2*K) = D1(2*K)/nrm
    
    nrm = dble(A(2,2))*D1(2*K+1) - aimag(A(2,2))*D1(2*K+2)
    D1(2*K+2) = dble(A(2,2))*D1(2*K+2) + aimag(A(2,2))*D1(2*K+1)
    D1(2*K+1) = nrm
    nrm = sqrt(D1(2*K+1)**2 + D1(2*K+2)**2)
    D1(2*K+1) = D1(2*K+1)/nrm
    D1(2*K+2) = D1(2*K+2)/nrm
    
  ! compute generalized Schur decomposition
  else
  
    ! generalized schur decomposition
    call z_2x2array_geneig('G',A,B,Wt,Vt,INFO)
      
    ! check INFO in debug mode
    if ((DEBUG).AND.(INFO.NE.0)) then
      call u_infocode_check(__FILE__,__LINE__,"z_2x2array_geneig",INFO,1)
      return 
    end if
    
    ! replace Vt with rotation G1
    call z_rot3_vec4gen(dble(Vt(1,1)),aimag(Vt(1,1)),dble(Vt(2,1)),aimag(Vt(2,1)),G1(1),G1(2),G1(3),nrm,INFO)
    
    ! check INFO in debug mode
    if ((DEBUG).AND.(INFO.NE.0)) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec4gen",INFO,1)
      return 
    end if
    
    ! replace Wt with rotation G2
    call z_rot3_vec4gen(dble(Wt(1,1)),aimag(Wt(1,1)),dble(Wt(2,1)),aimag(Wt(2,1)),G2(1),G2(2),G2(3),nrm,INFO)
    
    ! check INFO in debug mode
    if ((DEBUG).AND.(INFO.NE.0)) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec4gen",INFO,1)
      return 
    end if

    ! pass G1 through right triangular factor
    G3 = G1
    call z_upr1fact_rot3throughtri('R2L',N,K,D2,C2,B2,G3,INFO)
    
    ! check INFO in debug mode
    if ((DEBUG).AND.(INFO.NE.0)) then
      call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_rot3throughtri",INFO,1)
      return 
    end if
    
    ! equivalence transform
    A(1,1) = cmplx(G2(1),G2(2),kind=8)
    A(2,1) = cmplx(G2(3),0d0,kind=8)
    A(1,2) = -A(2,1)
    A(2,2) = conjg(A(1,1))
    
    B(1,1) = cmplx(G3(1),G3(2),kind=8)
    B(2,1) = cmplx(G3(3),0d0,kind=8)
    B(1,2) = -B(2,1)
    B(2,2) = conjg(B(1,1))
    
    A = matmul(transpose(conjg(A)),B)
    
    ! deflate into right diagonal D
    nrm = dble(A(1,1))*D2(2*K-1) - aimag(A(1,1))*D2(2*K)
    D2(2*K) = dble(A(1,1))*D2(2*K) + aimag(A(1,1))*D2(2*K-1)
    D2(2*K-1) = nrm
    nrm = sqrt(D2(2*K-1)**2 + D2(2*K)**2)
    D2(2*K-1) = D2(2*K-1)/nrm
    D2(2*K) = D2(2*K)/nrm
    
    nrm = dble(A(2,2))*D2(2*K+1) - aimag(A(2,2))*D2(2*K+2)
    D2(2*K+2) = dble(A(2,2))*D2(2*K+2) + aimag(A(2,2))*D2(2*K+1)
    D2(2*K+1) = nrm
    nrm = sqrt(D2(2*K+1)**2 + D2(2*K+2)**2)
    D2(2*K+1) = D2(2*K+1)/nrm
    D2(2*K+2) = D2(2*K+2)/nrm
    
    ! pass G1 through left triangular factor
    G3 = G1
    call z_upr1fact_rot3throughtri('R2L',N,K,D1,C1,B1,G3,INFO)
    
    ! check INFO in debug mode
    if ((DEBUG).AND.(INFO.NE.0)) then
      call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_rot3throughtri",INFO,1)
      return 
    end if
    
    ! equivalence transform of Q
    A(1,1) = cmplx(Q(3*K-2),Q(3*K-1),kind=8)
    A(2,1) = cmplx(Q(3*K),0d0,kind=8)
    A(1,2) = -A(2,1)
    A(2,2) = conjg(A(1,1))
    
    B(1,1) = cmplx(G3(1),G3(2),kind=8)
    B(2,1) = cmplx(G3(3),0d0,kind=8)
    B(1,2) = -B(2,1)
    B(2,2) = conjg(B(1,1))
    
    A = matmul(A,B)
    
    B(1,1) = cmplx(G2(1),-G2(2),kind=8)
    B(2,1) = cmplx(-G2(3),0d0,kind=8)
    B(1,2) = -B(2,1)
    B(2,2) = conjg(B(1,1))
    
    A = matmul(B,A)
    
    Q(3*K-2) = 1d0
    Q(3*K-1) = 0d0
    Q(3*K) = 0d0
    
    ! deflate into left diagonal D
    nrm = dble(A(1,1))*D1(2*K-1) - aimag(A(1,1))*D1(2*K)
    D1(2*K) = dble(A(1,1))*D1(2*K) + aimag(A(1,1))*D1(2*K-1)
    D1(2*K-1) = nrm
    nrm = sqrt(D1(2*K-1)**2 + D1(2*K)**2)
    D1(2*K-1) = D1(2*K-1)/nrm
    D1(2*K) = D1(2*K)/nrm
    
    nrm = dble(A(2,2))*D1(2*K+1) - aimag(A(2,2))*D1(2*K+2)
    D1(2*K+2) = dble(A(2,2))*D1(2*K+2) + aimag(A(2,2))*D1(2*K+1)
    D1(2*K+1) = nrm
    nrm = sqrt(D1(2*K+1)**2 + D1(2*K+2)**2)
    D1(2*K+1) = D1(2*K+1)/nrm
    D1(2*K+2) = D1(2*K+2)/nrm
    
    ! update vecs
    if (COMPZ.NE.'N') then
    
      ! update V
      B(1,1) = cmplx(G1(1),G1(2),kind=8)
      B(2,1) = cmplx(G1(3),0d0,kind=8)
      B(1,2) = -B(2,1)
      B(2,2) = conjg(B(1,1))
      
      V(:,K:(K+1)) = matmul(V(:,K:(K+1)),B) 
      
      ! update W
      B(1,1) = cmplx(G2(1),G2(2),kind=8)
      B(2,1) = cmplx(G2(3),0d0,kind=8)
      B(1,2) = -B(2,1)
      B(2,2) = conjg(B(1,1))
      
      W(:,K:(K+1)) = matmul(W(:,K:(K+1)),B)
    
    end if
  
  end if

end subroutine z_upr1fact_2x2deflation
