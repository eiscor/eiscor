#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DOFDFI (Double Orthogonal hessenberg Factored Doubleshift Francis Iteration)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' doubleshift 
! algorithm on an orthogonal upper hessenberg matrix that is stored as a 
! product of givens rotations and a diagonal matrix. 
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
!  Q               REAL(8) array of dimension (2*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  Z               REAL(8) array of dimension (N,N)
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
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies COMPZ is invalid
!                    INFO = -3 implies STR is invalid
!                    INFO = -4 implies STP is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DOFDFI(COMPZ,N,STR,STP,Q,D,Z,ITCNT,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: COMPZ
  integer, intent(in) :: N, STR, STP
  integer, intent(inout) :: ITCNT, INFO
  real(8), intent(inout) :: Q(2*(N-1)), D(2*N), Z(N,N)
  
  ! compute variables
  integer :: ii, ind
  real(8) :: s1, s2
  real(8) :: b1(2), b2(2), b3(2), temp(2)
  real(8) :: block(2,2) 
  complex(8) :: eigs(2), h(2,2)
  
  ! initialize INFO
  INFO = 0
 
  ! check input in debug mode
  if (DEBUG) then
  
    ! check COMPZ
    if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
      INFO = -1
      call UARERR(__FILE__,__LINE__,"COMPZ must be 'N', 'I' or 'V'",INFO,INFO)
      return
    end if
    
    ! check STR
    if ((STR < 1).OR.(STR > N-1)) then
      INFO = -3
      call UARERR(__FILE__,__LINE__,"STR must be 1 <= STR <= N-1",INFO,INFO)
      return
    end if
    
    ! check STP
    if ((STP < STR).OR.(STP > N-1)) then
      INFO = -4
      call UARERR(__FILE__,__LINE__,"STP must be STR <= STP <= N-1",INFO,INFO)
      return
    end if
  
  end if
  
  ! compute a nonzero shift
  ! shifts +1
  if((mod(ITCNT+1,11) == 0).AND.(ITCNT.NE.0))then
    eigs(1) = cmplx(1d0,0d0,kind=8)
    eigs(2) = cmplx(1d0,0d0,kind=8)
    
  ! shifts -1
  else if((mod(ITCNT+1,16) == 0).AND.(ITCNT.NE.0))then
    eigs(1) = cmplx(-1d0,0d0,kind=8)
    eigs(2) = cmplx(-1d0,0d0,kind=8)
    
  ! random shift
  else if((mod(ITCNT+1,21) == 0).AND.(ITCNT.NE.0))then
    call random_number(s1)
    call random_number(s2)
    eigs(1) = cmplx(s1,s2,kind=8)
    eigs(2) = cmplx(s1,-s2,kind=8)
          
  ! wilkinson shifts
  else
    ! get 2x2 block
    call DOFTDB(N,STP,Q,D,block,INFO)
      
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DOFTDB failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
        
    ! compute eigenvalues and eigenvectors
    call DTTEEV(block,eigs,h,INFO)
      
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DTTEEV failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if

  end if
  
  ! two real shifts
  if ((aimag(eigs(1)).EQ.0d0).AND.(aimag(eigs(2)).EQ.0d0)) then
 
    ! first shift
    eigs(1) = cmplx(sign(1d0,dble(eigs(1))),0d0,kind=8)
    
    ! build first bulge
    temp(1) = dble(eigs(1))
    temp(2) = aimag(eigs(1))
    call DOFCFT('S',N,STR,Q,D,temp,b1,b2,INFO)
          
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DOFCFT failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if

    ! fusion to initialize first bulge
    ind = 2*(STR-1)
    b3(1) = b1(1)
    b3(2) = -b1(2)
    call DARFGR('R',b3,Q((ind+1):(ind+2)),INFO)
     
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARFGR failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! first bulge update eigenvectors
    if (COMPZ .NE. 'N')then
      h(1,1) = b1(1)
      h(2,1) = b1(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,STR:(STR+1)) = matmul(Z(:,STR:(STR+1)),h)
    end if 
    
    ! first bulge through D
    call DARGTD('R',D((ind+1):(ind+4)),b1,INFO)
      
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARGTD failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! first bulge through Q
    call DARGTO(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b1,INFO)
      
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARGTO failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! second shift
    eigs(2) = cmplx(sign(1d0,dble(eigs(2))),0d0,kind=8)
  
    ! build bulge
    temp(1) = dble(eigs(2))
    temp(2) = aimag(eigs(2))
    call DOFCFT('S',N,STR,Q,D,temp,b2,b3,INFO)
          
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DOFCFT failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
   
    ! fusion to initialize second bulge
    b3(1) = b2(1)
    b3(2) = -b2(2)
    call DARFGR('R',b3,Q((ind+1):(ind+2)),INFO)
     
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARFGR failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! main chasing loop
    do ii=STR,(STP-2)
       
      ! set ind
      ind = 2*(ii-1)
       
      ! first bulge update eigenvectors
      if (COMPZ .NE. 'N')then
        h(1,1) = b1(1)
        h(2,1) = b1(2)
        h(1,2) = -h(2,1)
        h(2,2) = h(1,1)
        Z(:,(ii+1):(ii+2)) = matmul(Z(:,(ii+1):(ii+2)),h)
      end if 
        
      ! first bulge through D
      call DARGTD('R',D((ind+3):(ind+6)),b1,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call UARERR(__FILE__,__LINE__,"DARGTD failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
    
      ! first bulge through Q
      call DARGTO(Q((ind+3):(ind+4)),Q((ind+5):(ind+6)),b1,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call UARERR(__FILE__,__LINE__,"DARGTO failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
         
      ! second bulge update eigenvectors
      if (COMPZ .NE. 'N')then
        h(1,1) = b2(1)
        h(2,1) = b2(2)
        h(1,2) = -h(2,1)
        h(2,2) = h(1,1)
        Z(:,ii:(ii+1)) = matmul(Z(:,ii:(ii+1)),h)
      end if 
        
      ! second bulge through D
      call DARGTD('R',D((ind+1):(ind+4)),b2,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call UARERR(__FILE__,__LINE__,"DARGTD failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
    
      ! second bulge through Q
      call DARGTO(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b2,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call UARERR(__FILE__,__LINE__,"DARGTO failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
       
    end do
    
    ! set ind
    ind = 2*(stp-2)
    
    ! first bulge update eigenvectors
    if (COMPZ .NE. 'N')then
      h(1,1) = b1(1)
      h(2,1) = b1(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,STP:(STP+1)) = matmul(Z(:,STP:(STP+1)),h)
    end if 
       
    ! first bulge through D
    call DARGTD('R',D((ind+3):(ind+6)),b1,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARGTD failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! first bulge fuse with Q
    call DARFGR('L',Q((ind+3):(ind+4)),b1,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARFGR failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
        
    ! second bulge update eigenvectors
    if (COMPZ .NE. 'N')then
      h(1,1) = b2(1)
      h(2,1) = b2(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,(STP-1):STP) = matmul(Z(:,(STP-1):STP),h)
    end if
        
    ! second bulge through D
    call DARGTD('R',D((ind+1):(ind+4)),b2,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARGTD failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! second bulge through Q  
    call DARGTO(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b2,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARGTO failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
        
    ! set index
    ind = 2*(stp-1)
    
    ! last bulge update eigenvectors
    if (COMPZ .NE. 'N')then
      h(1,1) = b2(1)
      h(2,1) = b2(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,STP:(STP+1)) = matmul(Z(:,STP:(STP+1)),h)
    end if 
  
    ! last bulge through D
    call DARGTD('R',D((ind+1):(ind+4)),b2,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARGTD failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
  
    ! last bulge fuse with Q
    call DARFGR('L',Q((ind+1):(ind+2)),b2,INFO)
  
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARFGR failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
 
    ! update ITCNT
    ITCNT = ITCNT + 1

  ! complex conjugate pair
  else
  
    ! normalize shifts
    if (abs(eigs(1)).EQ.0d0) then
      eigs(1) = cmplx(0d0,1d0,kind=8)
    end if
    eigs(1) = eigs(1)/abs(eigs(1)) 

    ! build bulge
    temp(1) = dble(eigs(1))
    temp(2) = aimag(eigs(1))
    call DOFCFT('D',N,STR,Q,D,temp,b1,b2,INFO)
          
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DOFCFT failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if

    ! turnover to initialize bulge
    ind = 2*(STR-1)
    temp(1) = b2(1)
    temp(2) = -b2(2)
    b3(1) = b1(1)
    b3(2) = -b1(2)
    call DARGTO(temp,b3,Q((ind+1):(ind+2)),INFO)
      
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARGTO failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if

    ! fusion to finish initialization
    call DARFGR('R',b3,Q((ind+3):(ind+4)),INFO)
     
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARFGR failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
      
    ! update b3
    b3 = Q((ind+1):(ind+2))
    
    ! update Q
    Q((ind+1):(ind+2)) = temp
    
    ! main chasing loop
    do ii=str,(stp-2)
       
      ! set ind
      ind = 2*(ii-1)
       
      ! first bulge update eigenvectors
      if (COMPZ .NE. 'N')then
        h(1,1) = b1(1)
        h(2,1) = b1(2)
        h(1,2) = -h(2,1)
        h(2,2) = h(1,1)
        Z(:,(ii+1):(ii+2)) = matmul(Z(:,(ii+1):(ii+2)),h)
      end if 
        
      ! first bulge through D
      call DARGTD('R',D((ind+3):(ind+6)),b1,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call UARERR(__FILE__,__LINE__,"DARGTD failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
    
      ! first bulge through Q
      call DARGTO(Q((ind+3):(ind+4)),Q((ind+5):(ind+6)),b1,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call UARERR(__FILE__,__LINE__,"DARGTO failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
         
      ! second bulge update eigenvectors
      if (COMPZ .NE. 'N')then
        h(1,1) = b2(1)
        h(2,1) = b2(2)
        h(1,2) = -h(2,1)
        h(2,2) = h(1,1)
        Z(:,ii:(ii+1)) = matmul(Z(:,ii:(ii+1)),h)
      end if 
        
      ! second bulge through D
      call DARGTD('R',D((ind+1):(ind+4)),b2,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call UARERR(__FILE__,__LINE__,"DARGTD failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
    
      ! second bulge through Q
      call DARGTO(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b2,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call UARERR(__FILE__,__LINE__,"DARGTO failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
      
      ! push b3 down
      call DARGTO(b3,b1,b2,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call UARERR(__FILE__,__LINE__,"DARGTO failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
  
    ! update bulges
      temp = b2
      b2 = b3
      b3 = b1
      b1 = temp
       
    end do
    
    ! set ind
    ind = 2*(stp-2)
    
    ! first bulge update eigenvectors
    if (COMPZ .NE. 'N')then
      h(1,1) = b1(1)
      h(2,1) = b1(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,STP:(STP+1)) = matmul(Z(:,STP:(STP+1)),h)
    end if 
       
    ! first bulge through D
    call DARGTD('R',D((ind+3):(ind+6)),b1,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARGTD failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! first bulge fuse with Q
    call DARFGR('L',Q((ind+3):(ind+4)),b1,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARFGR failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
        
    ! second bulge update eigenvectors
    if (COMPZ .NE. 'N')then
      h(1,1) = b2(1)
      h(2,1) = b2(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,(STP-1):STP) = matmul(Z(:,(STP-1):STP),h)
    end if
        
    ! second bulge through D
    call DARGTD('R',D((ind+1):(ind+4)),b2,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARGTD failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! second bulge through Q  
    call DARGTO(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b2,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARGTO failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
        
    ! fuse b2 and b3
    call DARFGR('L',b3,b2,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARFGR failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
        
    ! set index
    ind = 2*(stp-1)
    
    ! last bulge update eigenvectors
    if (COMPZ .NE. 'N')then
      h(1,1) = b3(1)
      h(2,1) = b3(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,STP:(STP+1)) = matmul(Z(:,STP:(STP+1)),h)
    end if 
  
    ! last bulge through D
    call DARGTD('R',D((ind+1):(ind+4)),b3,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARGTD failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
  
    ! last bulge fuse with Q
    call DARFGR('L',Q((ind+1):(ind+2)),b3,INFO)
  
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"DARFGR failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
 
    ! update ITCNT
    ITCNT = ITCNT + 1
  
  end if ! end of doubleshift step

end subroutine DOFDFI
