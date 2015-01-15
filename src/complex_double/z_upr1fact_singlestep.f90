#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_singlestep 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift 
! algorithm on a upr1 pencil. 
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
!                    'I' or 'V': update schurvectors
!
!  N               INTEGER
!                    dimension of matrix
!
!  STR             INTEGER
!                    index of the top most givens rotation where 
!                    the iteration begins
!
!  STP             INTEGER
!                    index of the bottom most givens rotation where 
!                    the iteration ends
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  FUN             LOGICAL FUNCTION FUN(N,P)
!                    takes integer N and logical array P of 
!                    dimension N-1 and outputs a logical 
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*(N+1))
!                    array of generators for complex diagonal matrices
!
!  R               REAL(8) array of dimension (4,3*N)
!                    array of generators for upper-triangular parts
!                    of the pencil
!
!  V              COMPLEX(8) array of dimension (N,N)
!                   right schur vectors
!                   if COMPZ = 'N' unused
!                   if COMPZ = 'I' or 'V' update V to store right schurvectors 
!
!  W              COMPLEX(8) array of dimension (N,N)
!                   left schur vectors
!                   if COMPZ = 'N' unused
!                   if COMPZ = 'I' or 'V' update W to store left schurvectors
!                   if ALG = 'QR' unused
!
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies ALG is invalid
!                   INFO = -2 implies COMPZ is invalid
!                   INFO = -3 implies N is invalid
!                   INFO = -4 implies STR is invalid
!                   INFO = -5 implies STP is invalid
!                   INFO = -8 implies Q is invalid
!                   INFO = -9 implies D is invalid
!                   INFO = -10 implies R is invalid
!                   INFO = -11 implies V is invalid
!                   INFO = -12 implies W is invalid
!                   INFO = -13 implies ITCNT is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_singlestep(ALG,COMPZ,N,STR,STP,P,FUN,Q,D,R,V,W,ITCNT,INFO)

  implicit none
  
  ! input variables
  character(2), intent(in) :: ALG
  character, intent(in) :: COMPZ
  integer, intent(in) :: N, STR, STP
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D(2,2*(N+1)), R(4,3*N)
  complex(8), intent(inout) :: V(N,N), W(N,N)
  integer, intent(inout) :: INFO, ITCNT
  interface
    logical function FUN(N,P)
      integer, intent(in) :: N
      logical, intent(in) :: P(N-2)
    end function FUN
  end interface
  
  ! compute variables
  integer :: ii
  logical :: final_flag
  real(8) :: G1(3), G2(3), G3(3)
  complex(8) :: shift, rho, A(2,2), B(2,2), Vt(2,2), Wt(2,2)
  
  ! initialize INFO
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
  
    ! check factorization
    call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    if (INFO.EQ.-1) then
      call u_infocode_check(__FILE__,__LINE__,"ALG must be 'QR' or 'QZ'",INFO,-1)
      return
    end if
    if (INFO.EQ.-2) then
      call u_infocode_check(__FILE__,__LINE__,"N is invalid",INFO,-3)
      return
    end if
    if (INFO.EQ.-3) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,-8)
      return
    end if
    if (INFO.EQ.-4) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,-9)
      return
    end if
    if (INFO.EQ.-5) then
      call u_infocode_check(__FILE__,__LINE__,"R is invalid",INFO,-10)
      return
    end if
  
    ! check COMPZ
    if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
      INFO = -2
      call u_infocode_check(__FILE__,__LINE__,"COMPZ must be 'N', 'I' or 'V'",INFO,INFO)
      return
    end if
    
    ! check STR
    if ((STR < 1).OR.(STR > N-1)) then
      INFO = -4
      call u_infocode_check(__FILE__,__LINE__,"STR must 1 <= STR <= N-1",INFO,INFO)
      return
    end if 
    
    ! check STP
    if ((STP < STR).OR.(STP > N-1)) then
      INFO = -5
      call u_infocode_check(__FILE__,__LINE__,"STP must STR <= STP <= N-1",INFO,INFO)
      return
    end if  
    
    ! check V
    if ((COMPZ.EQ.'I').OR.(COMPZ.EQ.'V')) then
      call z_2Darray_check(N,N,V,INFO)
      if (INFO.NE.0) then
        call u_infocode_check(__FILE__,__LINE__,"V is invalid",INFO,-11)
        return
      end if 
    end if 
    
    ! check W
    if (ALG.EQ.'QZ') then
      if ((COMPZ.EQ.'I').OR.(COMPZ.EQ.'V')) then
        call z_2Darray_check(N,N,W,INFO)
        if (INFO.NE.0) then
          call u_infocode_check(__FILE__,__LINE__,"W is invalid",INFO,-12)
          return
        end if 
      end if
    end if 
    
    ! check ITCNT
    if (ITCNT < 0) then
      INFO = -13
      call u_infocode_check(__FILE__,__LINE__,"ITCNT must be non-negative",INFO,INFO)
      return
    end if 

  end if
  
  ! compute final_flag
  final_flag = FUN(N,P)
  
  ! compute shift
  ! random shift
  if(mod(ITCNT+1,11) == 0)then
    call random_number(G1(1))
    call random_number(G1(2))
    shift = cmplx(G1(1),G1(2),kind=8)
          
  ! wilkinson shift
  else
  
    ! get 2x2 blocks
    call z_upr1fact_2x2diagblocks(N,STP,ALG,P,Q,D,R,A,B,INFO)
      
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_2x2diagblocks failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! store bottom right entries
    shift = A(2,2)
    rho = B(2,2)
        
    ! compute eigenvalues and eigenvectors
    call z_2x2array_geneig('G',A,B,Wt,Vt,INFO)
      
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_2x2array_geneig failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
          
    ! choose wikinson shift
    ! complex abs does not matter here
    if(abs(A(2,2)-shift)+abs(B(2,2)-rho) < abs(A(1,1)-shift)+abs(B(1,1)-rho))then
      shift = A(2,2)
      rho = B(2,2)
    else
      shift = A(1,1)
      rho = B(1,1)
    end if
    
    ! avoid zero division
    if ((dble(rho).EQ.0).AND.(aimag(rho).EQ.0)) then
      shift = 1d16 ! not sure if this is a good idea?
    else
      shift = shift/rho
    end if

  end if

  ! build bulge
  call z_upr1fact_buildbulge(ALG,N,STR,P,Q,D,R,shift,G2,INFO)
        
  ! check INFO in debug mode
  if (DEBUG) then
    call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_buildbulge failed",INFO,INFO)
    if (INFO.NE.0) then 
      return 
    end if 
  end if
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! iteration for QR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (ALG.EQ.'QR') then
  
    ! update V
    if (COMPZ.NE.'N') then
      
      A(1,1) = cmplx(G1(1),G1(2),kind=8)
      A(2,1) = cmplx(G1(3),0d0,kind=8)
      A(1,2) = -A(2,1)
      A(2,2) = conjg(A(1,1))
      
      V(:,STR:(STR+1)) = matmul(V(:,STR:(STR+1)),A)
      
    end if
    
    ! set G1 as G2^-1 for turnover
    G1(1) = G2(1)
    G1(2) = -G2(2)
    G1(3) = -G2(3)
    
    ! merge with Q if necessary
    if (P(STR).EQV..FALSE.) then
    
      ! merge from left
      call z_upr1fact_mergebulge('L',N,STR,STP,STR,P,Q,D,G1,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_mergebulge failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
      
      ! set G1 for turnover
      G1(1) = Q(3*STR-2)
      G1(2) = Q(3*STR-1)
      G1(3) = Q(3*STR)

    end if
    
    ! pass G2 through triangular part
    call z_upr1fact_rot3throughtri('R2L',N,STR,D(1,:),R(1,:),R(2,:),G2,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_rot3throughtri failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! set G3 for turnover
    G3 = G2
    
    ! merge with Q if necessary
    if (P(STR).EQV..TRUE.) then
    
      ! merge from left
      call z_upr1fact_mergebulge('R',N,STR,STP,STR,P,Q,D,G2,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_mergebulge failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
      
      ! set G3 for turnover
      G3(1) = Q(3*STR-2)
      G3(2) = Q(3*STR-1)
      G3(3) = Q(3*STR)

    end if
    
    ! set G2 for turnover
    G2(1) = Q(3*STR+1)
    G2(2) = Q(3*STR+2)
    G2(3) = Q(3*STR+3)    
  
    ! chase bulge
    do ii=STR,(STP-2)
    
      ! execute turnover of G1G2G3
      call z_rot3_turnover(G1,G2,G3,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_rot3_turnover failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
      
      ! set Q(ii)
      Q(3*ii-2) = G1(1)
      Q(3*ii-1) = G1(2)
      Q(3*ii) = G1(3)
      
      ! prepare for next turnover based on P(ii+1)
      ! hess
      if (P(ii+1).EQV..FALSE.) then
      
        ! set P(ii)
        P(ii) = P(ii+1)
        
        ! set Q(ii+1)
        Q(3*ii+1) = G2(1)
        Q(3*ii+2) = G2(2)
        Q(3*ii+3) = G2(3)        
 
        ! set G1 for turnover
        G1 = G2     
        
        ! set G2 for turnover
        G2(1) = Q(3*ii+4)
        G2(2) = Q(3*ii+5)
        G2(3) = Q(3*ii+6)
        
        ! update V
        if (COMPZ.NE.'N') then
          
          A(1,1) = cmplx(G3(1),G3(2),kind=8)
          A(2,1) = cmplx(G3(3),0d0,kind=8)
          A(1,2) = -A(2,1)
          A(2,2) = conjg(A(1,1))
          
          V(:,(ii+1):(ii+2)) = matmul(V(:,(ii+1):(ii+2)),A)
          
        end if
        
        ! pass G3 through upper triangular part
        call z_upr1fact_rot3throughtri('R2L',N,ii+1,D(1,:),R(1,:),R(2,:),G3,INFO)
        
        ! check INFO in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_rot3throughtri failed",INFO,INFO)
          if (INFO.NE.0) then 
            return 
          end if 
        end if        
      
      ! inverse hess
      else
      
        ! set P(ii)
        P(ii) = P(ii+1)
        
        ! set Q(ii+1)
        Q(3*ii+1) = G3(1)
        Q(3*ii+2) = G3(2)
        Q(3*ii+3) = G3(3)  
        
        ! pass G2 through upper triangular part
        call z_upr1fact_rot3throughtri('L2R',N,ii+1,D(1,:),R(1,:),R(2,:),G2,INFO)
        
        ! check INFO in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_rot3throughtri failed",INFO,INFO)
          if (INFO.NE.0) then 
            return 
          end if 
        end if  
        
        ! update V
        if (COMPZ.NE.'N') then
          
          A(1,1) = cmplx(G2(1),-G2(2),kind=8)
          A(2,1) = cmplx(-G2(3),0d0,kind=8)
          A(1,2) = -A(2,1)
          A(2,2) = conjg(A(1,1))
          
          V(:,(ii+1):(ii+2)) = matmul(V(:,(ii+1):(ii+2)),A)
          
        end if

        ! set G1 for turnover
        G1 = G2    
        
        ! set G2 for turnover
        G2(1) = Q(3*ii+4)
        G2(2) = Q(3*ii+5)
        G2(3) = Q(3*ii+6)
        
      end if

    end do
    
    ! final turnover
    call z_rot3_turnover(G1,G2,G3,INFO)
      
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_turnover failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! set P(STP-1)
    P(STP-1) = final_flag
    
    ! finish transformation based on P(STP-1)
    ! hess
    if (P(STP-1).EQV..FALSE.) then
    
      ! set Q(STP-1)
      Q(3*(STP-1)-2) = G1(1)
      Q(3*(STP-1)-1) = G1(2)
      Q(3*(STP-1)) = G1(3)   
      
      ! set Q(STP)
      Q(3*STP-2) = G2(1)
      Q(3*STP-1) = G2(2)
      Q(3*STP) = G2(3)  
      
      ! update V
      if (COMPZ.NE.'N') then
          
        A(1,1) = cmplx(G3(1),G3(2),kind=8)
        A(2,1) = cmplx(G3(3),0d0,kind=8)
        A(1,2) = -A(2,1)
        A(2,2) = conjg(A(1,1))
          
        V(:,(STP):(STP+1)) = matmul(V(:,(STP):(STP+1)),A)
         
      end if    
    
      ! pass G3 through upper triangular part
      call z_upr1fact_rot3throughtri('R2L',N,STP,D(1,:),R(1,:),R(2,:),G3,INFO)
        
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_rot3throughtri failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if  
      
      ! merge bulge 
      call z_upr1fact_mergebulge('R',N,STR,STP,STP,P,Q,D,G3,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_mergebulge failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if     
    
    ! inverse hess
    else
    
      ! set Q(STP-1)
      Q(3*(STP-1)-2) = G1(1)
      Q(3*(STP-1)-1) = G1(2)
      Q(3*(STP-1)) = G1(3)   
      
      ! set Q(STP)
      Q(3*STP-2) = G3(1)
      Q(3*STP-1) = G3(2)
      Q(3*STP) = G3(3)  
      
      ! pass G2 through upper triangular part
      call z_upr1fact_rot3throughtri('L2R',N,STP,D(1,:),R(1,:),R(2,:),G2,INFO)
        
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_rot3throughtri failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if 
      
      ! update V
      if (COMPZ.NE.'N') then
          
        A(1,1) = cmplx(G2(1),-G2(2),kind=8)
        A(2,1) = cmplx(-G2(3),0d0,kind=8)
        A(1,2) = -A(2,1)
        A(2,2) = conjg(A(1,1))
          
        V(:,(STP):(STP+1)) = matmul(V(:,(STP):(STP+1)),A)
         
      end if  
      
      ! merge bulge 
      call z_upr1fact_mergebulge('L',N,STR,STP,STP,P,Q,D,G2,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_mergebulge failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if 
    
    end if
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! iteration for QZ
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  else
  
    ! chase bulge
    do ii=STR,(STP-1)
    

    end do  
  
  end if
  
  ! update ITCNT
  ITCNT = ITCNT + 1

end subroutine z_upr1fact_singlestep
