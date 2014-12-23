#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_deflationcheck 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks for deflations in a unitary upper hessenberg 
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
!  Q               REAL(8) array of dimension (3*(N-1))
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
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies N is invalid
!                   INFO = -2 implies STR is invalid
!                   INFO = -3 implies STP is invalid
!                   INFO = -4 implies ZERO is invalid
!                   INFO = -8 implies ITCNT is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_deflationcheck(N,STR,STP,ZERO,P,Q,D,ITCNT,ITS,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: STR, STP, ZERO, ITCNT, INFO
  logical, intent(in) :: P(N-1)
  real(8), intent(inout) :: Q(3*(N-1)), D(2,2*N)
  integer, intent(inout) :: ITS(N-1)

  ! compute variables
  integer :: ii, ind1, ll, jj
  real(8), parameter :: tol = epsilon(1d0)
  real(8) :: qr, qi, d1r, d1i, c1r, c1i, s, nrm
  
  ! initialize info
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
    
    ! check N
    if (N < 2) then
      INFO = -1
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
      return
    end if
    
    ! check STR
    if ((STR < 1).OR.(STR > N-1)) then
      INFO = -2
      call u_infocode_check(__FILE__,__LINE__,"STR must 1 <= STR <= N-1",INFO,INFO)
      return
    end if 
    
    ! check STP
    if ((STP < STR).OR.(STP > N-1)) then
      INFO = -3
      call u_infocode_check(__FILE__,__LINE__,"STP must STR <= STP <= N-1",INFO,INFO)
      return
    end if  
    
    ! check ZERO
    if ((ZERO >= STR).OR.(ZERO < 0)) then
      INFO = -4
      call u_infocode_check(__FILE__,__LINE__,"ZERO must 0 <= ZERO < STR",INFO,INFO)
      return
    end if  
    
! check P
    
! check Q

! check D
    
    ! check ITCNT
    if (ITCNT < 0) then
      INFO = -7
      call u_infocode_check(__FILE__,__LINE__,"ITCNT must be non-negative",INFO,INFO)
      return
    end if 

  end if
  
  ! check for deflation
  do ii=1,(STP-STR+1)
  
    ! one less than the index of the rotaion being checked
    ind1 = 3*(STP-ii)
     
    ! deflate if subdiagonal is small enough
    nrm = abs(Q(ind1+3))
    if(nrm < tol)then
      
      ! extract diagonal
      qr = Q(ind1+1)
      qi = Q(ind1+2)
                
      ! set rotation to identity
      Q(ind1+1) = 1d0
      Q(ind1+2) = 0d0
      Q(ind1+3) = 0d0
      
        
      ! deflate upward
      do ll = 
        
                  
        ! update first diagonal
        ind1 = 2*(STP-ii)
        d1r = D(ind1+1)
        d1i = D(ind1+2)
          
        nrm = c1r*d1r - c1i*d1i
        d1i = c1r*d1i + c1i*d1r
        d1r = nrm
        nrm = sqrt(d1r*d1r + d1i*d1i)
        if (nrm.NE.0) then
          d1r = d1r/nrm
          d1i = d1i/nrm
        else
          d1r = 0d0
          d1i = 0d0
        end if
          
        D(ind1+1) = d1r
        D(ind1+2) = d1i
   
        
        ! exit loop if P == .FALSE.
        if ( .EQ..FALSE.)
         
        end if
        
      end do
        
      ! deflate downward
      do ll = 
        
                ! 1x1 deflation
        if(ii == 1)then
          
          ! update second diagonal
          ind1 = 2*(STP-ii)
          d1r = D(ind1+3)
          d1i = D(ind1+4)
             
          nrm = c1r*d1r + c1i*d1i
          d1i = c1r*d1i - c1i*d1r
          d1r = nrm
          nrm = sqrt(d1r*d1r + d1i*d1i)
          d1r = d1r/nrm
          d1i = d1i/nrm
             
          D(ind1+3) = d1r
          D(ind1+4) = d1i
             
          ! 2x2 or bigger
          else

            ! update Q
            do ll=(STP+1-ii),(STP-1)
              ind1 = 3*(ll)
              d1r = Q(ind1+1)
              d1i = Q(ind1+2)
              s = Q(ind1+3)
                
              nrm = c1r*d1r + c1i*d1i
              d1i = c1r*d1i - c1i*d1r
              d1r = nrm
              nrm = sqrt(d1r*d1r + d1i*d1i + s*s)
              d1r = d1r/nrm
              d1i = d1i/nrm
              s = s/nrm
                
              Q(ind1+1) = d1r
              Q(ind1+2) = d1i
              Q(ind1+3) = s
            end do
             
            ! update second diagonal
            ind1 = 2*(STP)
            d1r = D(ind1+1)
            d1i = D(ind1+2)
             
            nrm = c1r*d1r + c1i*d1i
            d1i = c1r*d1i - c1i*d1r
            d1r = nrm
            nrm = sqrt(d1r*d1r + d1i*d1i)
            d1r = d1r/nrm
            d1i = d1i/nrm
             
            D(ind1+1) = d1r
            D(ind1+2) = d1i
          end if
   
        ! exit loop if P == .TRUE.
        if ( .EQ..TRUE.)
          
        end if
        
      end do
       
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
  
end subroutine z_upr1fact_deflationcheck
