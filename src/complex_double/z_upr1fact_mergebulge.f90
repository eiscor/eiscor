#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_mergebulge 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine fuses a Givens rotation represented by 3 real numbers: 
! the real and imaginary parts of a complex cosine and a scrictly real 
! sine into the top or bottom of a unitary hessenberg matrix that is 
! stored as a product of givens rotations and a complex diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  JOB             CHARACTER
!                    'L': fuse from the left
!                    'R': fuse from the right
!
!  N               INTEGER
!                    dimension of matrix
!
!  STR             INTEGER
!                    index of the top most givens rotation where 
!                    fusion could happen
!
!  STP             INTEGER
!                    index of the bottom most givens rotation where 
!                    fusion could happen
!
!  K               INTEGER
!                    index where fusion should happen
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N+1)
!                    array of generators for complex diagonal matrix
!
!  B               REAL(8) array of dimension (3)
!                    generators for rotation that will be fused
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                   INFO = 1 implies a core transformation interferes
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies JOB is invalid
!                   INFO = -2 implies N is invalid
!                   INFO = -3 implies STR is invalid
!                   INFO = -4 implies STP is invalid
!                   INFO = -7 implies B is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_mergebulge(JOB,N,STR,STP,K,P,Q,D,G,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: JOB
  integer, intent(in) :: N, STR, STP, K
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)),D(2*(N+1))
  real(8), intent(in) :: G(3)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii
  
  ! initialize INFO
  INFO = 0

  ! check input in debug mode
  if (DEBUG) then
  
    ! check JOB
    if ((JOB.NE.'L').AND.(JOB.NE.'R')) then
      INFO = -1
      call u_infocode_check(__FILE__,__LINE__,"JOB must be 'L' or 'R'",INFO,INFO)
      return
    end if
    
    ! check N
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
    
    ! check K
    if ((K < STR).OR.(STP < K)) then
      INFO = -5
      call u_infocode_check(__FILE__,__LINE__,"STP must STR <= K <= STP",INFO,INFO)
      return
    end if   
    
    ! check Q
    call d_1Darray_check(3*(N-1),Q,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,-7)
      return
    end if 
    
    ! check D
    call d_1Darray_check(2*(N+1),D,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,-8)
      return
    end if 
    
    ! check G
    call d_1Darray_check(3,G,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"G is invalid",INFO,-9)
      return
    end if  
    
    ! check for valid fusion
    ! JOB == L
    if (JOB.EQ.'L') then
      
      ! check upper core transformation
      if (K > 1) then
        if (P(K-1).EQV..FALSE.) then
          INFO = 1 
          call u_infocode_check(__FILE__,__LINE__,"Upper core transformation interferes",INFO,INFO)
          return
        end if
      end if
      
      ! check lower core transformation
      if (K < (N-1)) then
        if (P(K).EQV..TRUE.) then
          INFO = 1 
          call u_infocode_check(__FILE__,__LINE__,"Lower core transformation interferes",INFO,INFO)
          return
        end if
      end if
      
     ! JOB == R
     else
     
      ! check upper core transformation
      if (K > 1) then
        if (P(K-1).EQV..TRUE.) then
          INFO = 1 
          call u_infocode_check(__FILE__,__LINE__,"Upper core transformation interferes",INFO,INFO)
          return
        end if
      end if
      
      ! check lower core transformation
      if (K < (N-1)) then
        if (P(K).EQV..FALSE.) then
          INFO = 1 
          call u_infocode_check(__FILE__,__LINE__,"Lower core transformation interferes",INFO,INFO)
          return
        end if
      end if
     
     end if     

  end if
  
  ! fusion from left
  if(JOB.EQ.'L')then
  

     
  ! fusion from right
  else
  

     
  end if

end subroutine z_upr1fact_mergebulge
