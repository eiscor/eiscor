#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_singlestep 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift 
! algorithm on a unitary upper hessenberg matrix that is stored as a 
! product of givens rotations and a complex diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  COMPZ           CHARACTER
!                    'N': no eigenvectors
!                    'I': eigenvectors, initializing Z to the identity
!                    'V': eigenvectors, assume Z already initialized
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
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  Z               COMPLEX(8) array of dimension (N,N)
!                    if COMPZ = 'N' unused
!                    if COMPZ = 'I' stores eigenvectors in Z 
!                    if COMPZ = 'V' update Z to store eigenvectors 
!
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies COMPZ is invalid
!                   INFO = -2 implies N is invalid
!                   INFO = -3 implies STR is invalid
!                   INFO = -4 implies STP is invalid
!                   INFO = -7 implies Z is invalid
!                   INFO = -8 implies ITCNT is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_singlestep(COMPZ,N,STR,STP,Q,D,Z,ITCNT,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: COMPZ
  integer, intent(in) :: N, STR, STP
  integer, intent(inout) :: ITCNT, INFO
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N)
  complex(8), intent(inout) :: Z(N,N)
  
  ! compute variables
  integer :: ii, ind1, ind2
  real(8) :: s1, s2
  real(8) :: bulge(3),binv(3)
  complex(8) :: shift
  complex(8) :: block(2,2), temp(2,2), eigs(2)
  
  ! initialize INFO
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
  
    ! check COMPZ
    if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
      INFO = -1
      call u_infocode_check(__FILE__,__LINE__,"COMPZ must be 'N', 'I' or 'V'",INFO,INFO)
      return
    end if
    
    ! check N
    call i_scalar_nancheck(N,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"N is invalid",INFO,-2)
      return
    end if
    call i_scalar_infcheck(N,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"N is invalid",INFO,-2)
      return
    end if
    if (N < 2) then
      INFO = -2
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
      return
    end if
    
    ! check STR
    if ((STR < 1).OR.(STR > N-1)) then
      INFO = -3
      call u_infocode_check(__FILE__,__LINE__,"STR must 1 <= STR <= N-1",INFO,INFO)
      return
    end if 
    
    ! check STP
    if ((STP < STR).OR.(STP > N-1)) then
      INFO = -4
      call u_infocode_check(__FILE__,__LINE__,"STP must STR <= STP <= N-1",INFO,INFO)
      return
    end if  
    
    ! check Z
    if (COMPZ.EQ.'V') then
      call z_2Darray_check(N,N,Z,INFO)
      if (INFO.NE.0) then
        call u_infocode_check(__FILE__,__LINE__,"Z is invalid",INFO,-7)
        return
      end if 
    end if 
    
    ! check ITCNT
    call i_scalar_nancheck(ITCNT,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"ITCNT is invalid",INFO,-8)
      return
    end if
    call i_scalar_infcheck(ITCNT,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"ITCNT is invalid",INFO,-8)
      return
    end if
    if (ITCNT < 0) then
      INFO = -8
      call u_infocode_check(__FILE__,__LINE__,"ITCNT must be non-negative",INFO,INFO)
      return
    end if 

  end if
  
  ! compute a nonzero shift
  ! random shift
  if(mod(ITCNT+1,11) == 0)then
    call random_number(s1)
    call random_number(s2)
    shift = cmplx(s1,s2,kind=8)
          
  ! wilkinson shift
  else
    ! get 2x2 block
    call ZUFTDB(N,STP,Q,D,block,INFO) 
      
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"ZUFTDB failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
        
    ! compute eigenvalues and eigenvectors
    call z_2x2array_eig(block,eigs,temp,INFO)
      
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_2x2array_eig failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
          
    ! choose wikinson shift
    if(abs(block(2,2)-eigs(1)) < abs(block(2,2)-eigs(2)))then
      shift = eigs(1)
    else
      shift = eigs(2)
    end if

  end if
        
  ! project shift onto unit circle
  if (abs(shift) .EQ. 0d0) then
    shift = cmplx(1d0,0d0,kind=8)
  end if
  shift = shift/abs(shift)

  ! build bulge
  call z_unifact_buildbulge(N,STR,Q,D,shift,bulge,INFO)
        
  ! check INFO in debug mode
  if (DEBUG) then
    call u_infocode_check(__FILE__,__LINE__,"z_unifact_buildbulge failed",INFO,INFO)
    if (INFO.NE.0) then 
      return 
    end if 
  end if
  
  ! bulge inverse
  binv(1) = bulge(1)
  binv(2) = -bulge(2)
  binv(3) = -bulge(3)
  
  ! fusion at top
  call z_unifact_mergebulge('T',N,STR,STP,Q,D,binv,INFO)
  
  ! check INFO in debug mode
  if (DEBUG) then
    call u_infocode_check(__FILE__,__LINE__,"z_unifact_mergebulge failed",INFO,INFO)
    if (INFO.NE.0) then 
      return 
    end if 
  end if
  
  ! main chasing loop
  do ii=STR,(STP-1)
  
    ! update eigenvectors
    if (COMPZ .NE. 'N')then
      temp(1,1) = cmplx(bulge(1),bulge(2),kind=8)
      temp(2,1) = cmplx(bulge(3),0d0,kind=8)
      temp(1,2) = -temp(2,1)
      temp(2,2) = conjg(temp(1,1))
      Z(:,ii:(ii+1)) = matmul(Z(:,ii:(ii+1)),temp)
    end if
     
    ! set indices
    ind1 = 2*(ii-1) + 1
    ind2 = ind1+3
     
    ! through diag
    call z_rot3_swapdiag('R',D(ind1:ind2),bulge,INFO)

    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_swapdiag failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! set indices
    ind1 = 3*(ii-1) + 1
    ind2 = ind1+2
     
    ! through Q
    call z_rot3_turnover(Q(ind1:ind2),Q((ind1+3):(ind2+3)),bulge,INFO)

    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_turnover failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if

  end do
  
  ! update eigenvectors
  if (COMPZ .NE. 'N')then
    temp(1,1) = cmplx(bulge(1),bulge(2),kind=8)
    temp(2,1) = cmplx(bulge(3),0d0,kind=8)
    temp(1,2) = -temp(2,1)
    temp(2,2) = conjg(temp(1,1))
    Z(:,STP:(STP+1)) = matmul(Z(:,STP:(STP+1)),temp)
  end if
  
  ! fusion at bottom
  call z_unifact_mergebulge('B',N,STR,STP,Q,D,bulge,INFO)
  
  ! check INFO in debug mode
  if (DEBUG) then
    call u_infocode_check(__FILE__,__LINE__,"z_unifact_mergebulge failed",INFO,INFO)
    if (INFO.NE.0) then 
      return 
    end if 
  end if
  
  ! update ITCNT
  ITCNT = ITCNT + 1

end subroutine z_unifact_singlestep
