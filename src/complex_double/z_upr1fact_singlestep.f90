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
  complex(8) :: shift
  
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
  
  ! compute shift
  ! random shift
  if(mod(ITCNT+1,11) == 0)then
    call random_number(s1)
    call random_number(s2)
    shift = cmplx(s1,s2,kind=8)
          
  ! wilkinson shift
  else
    ! get 2x2 block
    call z_upr1fact_2x2diagblocks(N,STP,Q,D,block,INFO) 
      
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_2x2diagblocks failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
        
    ! compute eigenvalues and eigenvectors
    call z_2x2array_geneig(block,eigs,temp,INFO)
      
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_2x2array_geneig failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
          
    ! choose wikinson shift
    ! complex abs does not matter here
    if(abs(block(2,2)-eigs(1)) < abs(block(2,2)-eigs(2)))then
      shift = eigs(1)
    else
      shift = eigs(2)
    end if

  end if

  ! build bulge
  call z_upr1fact_buildbulge(N,STR,Q,D,shift,bulge,INFO)
        
  ! check INFO in debug mode
  if (DEBUG) then
    call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_buildbulge failed",INFO,INFO)
    if (INFO.NE.0) then 
      return 
    end if 
  end if
  
  ! main chasing loop
  do ii=STR,(STP-1)
  

  end do
  
  ! update ITCNT
  ITCNT = ITCNT + 1

end subroutine z_upr1fact_singlestep
