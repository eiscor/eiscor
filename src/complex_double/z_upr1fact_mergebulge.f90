#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_mergebulge 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine fuses a Givens rotation represented by 3 real numbers: 
! the real and imaginary parts of a complex cosine and a scrictly real 
! sine into the unitary extended hessenberg part of a upr1 pencil.
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
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N+1)
!                    array of generators for complex diagonal matrix
!
!  G               REAL(8) array of dimension (3)
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
!                   INFO = -5 implies K is invalid
!                   INFO = -7 implies Q is invalid
!                   INFO = -8 implies D is invalid
!                   INFO = -9 implies G is invalid
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
  integer :: ii, jj, ind, up, down
  real(8) :: c1r, c1i, s1
  real(8) :: c2r, c2i, s2
  real(8) :: c3r, c3i, s3r, s3i
  real(8) :: d1r, d1i
  real(8) :: d2r, d2i
  real(8) :: phr, phi
  real(8) :: nrm
  
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
      if (K > STR) then
        if (P(K-1).EQV..FALSE.) then
          INFO = 1 
          call u_infocode_check(__FILE__,__LINE__,"Upper core transformation interferes",INFO,INFO)
          return
        end if
      end if
      
      ! check lower core transformation
      if (K < STP) then
        if (P(K).EQV..TRUE.) then
          INFO = 1 
          call u_infocode_check(__FILE__,__LINE__,"Lower core transformation interferes",INFO,INFO)
          return
        end if
      end if
      
     ! JOB == R
     else
     
      ! check upper core transformation
      if (K > STR) then
        if (P(K-1).EQV..TRUE.) then
          INFO = 1 
          call u_infocode_check(__FILE__,__LINE__,"Upper core transformation interferes",INFO,INFO)
          return
        end if
      end if
      
      ! check lower core transformation
      if (K < STP) then
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
  
    ! set inputs  
    c2r = G(1)
    c2i = G(2)
    s2 = G(3)
  
    ! retrieve Q(K)  
    ind = 3*(K-1)
    c1r = Q(ind+1)
    c1i = Q(ind+2)
    s1 = Q(ind+3)

    ! compute product GQ(K)    
    c3r = c1r*c2r - c1i*c2i - s1*s2
    c3i = c1r*c2i + c1i*c2r
    s3r = s1*c2r + s2*c1r
    s3i = -(s1*c2i - s2*c1i)
    
    ! compute phase
    call d_rot2_vec2gen(s3r,s3i,phr,phi,nrm,INFO)

    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
      
    ! update Q
    c2r = c3r*phr + c3i*phi
    c2i = -c3r*phi + c3i*phr
    s2 = nrm
    
    call z_rot3_vec3gen(c2r,c2i,s2,Q(ind+1),Q(ind+2),Q(ind+3),nrm,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec3gen failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! initialize up
    up = K
  
    ! move phase upward
    do jj = 1,(K-STR)
        
      ! set upward index
      up = K-jj
        
      ! exit loop if P == .FALSE.
      if (P(up).EQV..FALSE.) then
        up = up + 1
        exit    
      end if
   
      ! update Q
      c1r = Q(3*up-2)
      c1i = Q(3*up-1)
      s1 = Q(3*up)
                
      nrm = phr*c1r + phi*c1i
      c1i = phr*c1i - phi*c1r
      c1r = nrm
        
      call z_rot3_vec3gen(c1r,c1i,s1,Q(3*up-2),Q(3*up-1),Q(3*up),nrm,INFO)
        
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec3gen failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
        
    end do
       
    ! update upward diagonal
    d1r = D(2*up-1)
    d1i = D(2*up)
            
    nrm = phr*d1r - phi*d1i
    d1i = phr*d1i + phi*d1r
    d1r = nrm

    call d_rot2_vec2gen(d1r,d1i,D(2*up-1),D(2*up),nrm,INFO)

    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if 
    
    ! initialize downward index
    down = K
        
    ! move phase downward
    do jj = 1,(STP-K)
        
      ! exit if P == .TRUE.
      if (P(down).EQV..TRUE.) then
        exit
      end if
                
      ! set downward index
      down = K + jj
        
      ! update Q
      c2r = Q(3*down-2)
      c2i = Q(3*down-1)
      s2 = Q(3*down)
                
      nrm = phr*c2r + phi*c2i
      c2i = phr*c2i - phi*c2r
      c2r = nrm

      call z_rot3_vec3gen(c2r,c2i,s2,Q(3*down-2),Q(3*down-1),Q(3*down),nrm,INFO)
        
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec3gen failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
                    
    end do
      
    ! update downward index
    down = down + 1
      
    ! update downward diagonal
    d2r = D(2*down-1)
    d2i = D(2*down)
            
    nrm = phr*d2r + phi*d2i
    d2i = phr*d2i - phi*d2r
    d2r = nrm
      
    call d_rot2_vec2gen(d2r,d2i,D(2*down-1),D(2*down),nrm,INFO)

    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if 
     
  ! fusion from right
  else
  
    ! set inputs  
    c2r = G(1)
    c2i = G(2)
    s2 = G(3)
  
    ! retrieve Q(K)  
    ind = 3*(K-1)
    c1r = Q(ind+1)
    c1i = Q(ind+2)
    s1 = Q(ind+3)
    
    ! compute product Q(K)G
    c3r = c1r*c2r - c1i*c2i - s1*s2
    c3i = c1r*c2i + c1i*c2r
    s3r = s1*c2r + s2*c1r
    s3i = s1*c2i - s2*c1i
    
    ! compute phase
    call d_rot2_vec2gen(s3r,s3i,phr,phi,nrm,INFO)

    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
      
    ! update Q
    c2r = c3r*phr + c3i*phi
    c2i = -c3r*phi + c3i*phr
    s2 = nrm
    
    call z_rot3_vec3gen(c2r,c2i,s2,Q(ind+1),Q(ind+2),Q(ind+3),nrm,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec3gen failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! retrieve D
    ind = 2*(K-1)
    d1r = D(ind+1)
    d1i = D(ind+2)
    d2r = D(ind+3)
    d2i = D(ind+4)
     
    ! update first entry of D
    c1r = phr*d1r - phi*d1i
    c1i = phr*d1i + phi*d1r
    
    call d_rot2_vec2gen(c1r,c1i,D(ind+1),D(ind+2),nrm,INFO)

    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! update second entry of D
    c2r = phr*d2r + phi*d2i
    c2i = phr*d2i - phi*d2r

    call d_rot2_vec2gen(c2r,c2i,D(ind+3),D(ind+4),nrm,INFO)

    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
  end if

end subroutine z_upr1fact_mergebulge
