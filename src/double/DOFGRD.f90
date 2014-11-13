#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DOFGRD (Double Orthogonal hessenberg Factored Givens Rotation Deflation)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks for deflations in an orthogonal upper hessenberg 
! matrix that is stored as a product of givens rotations and a complex 
! diagonal matrix. When a deflation occurs the corresponding rotation
! is set to the identity matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  STR             INTEGER
!                    index of the top most givens rotation where 
!                    deflation should be checked
!
!  STP             INTEGER
!                    index of the bottom most givens rotation where 
!                    deflation should be checked
!
!  Q               REAL(8) array of dimension (2*(N-1))
!                    array of generators for givens rotations
!                    generators must be orthogonal to working precision
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    generators must be orthogonal to working precision
!
!  ITCNT           INTEGER
!                    number of iterations since last deflation
!
! OUTPUT VARIABLES:
!
!  ZERO            INTEGER
!                     index of the last known deflation
!                     on output contains index of newest deflation
!
!  ITS             INTEGER array of dimension (N-1)
!                    Contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies N is invalid
!                    INFO = -2 implies STR is invalid
!                    INFO = -3 implies STP is invalid
!                    INFO = -7 implies ITCNT is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DOFGRD(N,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: STR, STP, ZERO, ITCNT, INFO
  real(8), intent(inout) :: Q(2*(N-1)), D(2*N)
  integer, intent(inout) :: ITS(N-1)

  ! compute variables
  integer :: ii, ind, ll, jj, k
  real(8), parameter :: tol = epsilon(1d0)
  real(8) :: td, tc, ts, nrm
  
  ! initialize info
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then  
  
    ! check N
    if (N < 2) then
      INFO = -1
      call UARERR(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
      return
    end if
    
    ! check STR
    if ((STR < 1).OR.(STR > N-1)) then
      INFO = -2
      call UARERR(__FILE__,__LINE__,"STR must be 1 <= STR <= N-1",INFO,INFO)
      return
    end if
    
    ! check STP
    if ((STP < STR).OR.(STP > N-1)) then
      INFO = -3
      call UARERR(__FILE__,__LINE__,"STP must be STR <= STP <= N-1",INFO,INFO)
      return
    end if
    
    ! check ITCNT
    if (ITCNT < 0) then
      INFO = -7
      call UARERR(__FILE__,__LINE__,"ITCNT must be positive",INFO,INFO)
      return
    end if
  
  end if
  
  ! check for deflation
  do ii=1,(STP-STR+1)
  
    ! one less than the index of the rotation being checked
    ind = 2*(STP-ii)
     
    ! deflate if subdiagonal is small enough
    nrm = abs(Q(ind+2))
    if(nrm < tol)then
     
      ! set sub-diagonal to 0
      Q(ind+2) = 0d0
        
      ! update first Q
      tc = sign(1d0,Q(ind+1))
      Q(ind+1) = 1d0
    
      ! return if tc == -1
      if (tc < 0d0) then
      
        ! update first diagonal  
        D(ind+1) = -D(ind+1)
        
        ! update second diagonal
        D(ind+3) = -D(ind+3)

        ! update rest of Q
        if ((ind+4) <= 2*(N-1)) then
          Q(ind+4) = -Q(ind+4)
        end if
        
      end if
      
      ! update indices
      ZERO = STP+1-ii
      STR = ZERO + 1
        
      ! store it_count
      ITS(ZERO) = ITCNT
      ITCNT = 0
      
      ! exit loop  
      exit

    end if
  end do
  

  
end subroutine DOFGRD
