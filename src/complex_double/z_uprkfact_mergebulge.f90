#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkfact_mergebulge 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine fuses a Givens rotation represented by 3 real numbers: 
! the real and imaginary parts of a complex cosine and a scrictly real 
! sine into the unitary extended hessenberg part of a upr1 pencil.
!
! There are only 4 possibilities for this fusion:
!
! 1)    top, left  <=>      TOP.AND..NOT.P(1)
!
! 2)    top, right <=>      TOP.AND.P(1)
!
! In case (1) and (2) it is assumed that the rotation in Q above is an
! identity, i.e., that there was a deflation.
!
! 3) bottom, left  <=> .NOT.TOP.AND.P(N-2)
!
! 4) bottom, right <=> .NOT.TOP.AND..NOT.P(N-2)
!
! In case (3) and (4) there is assumed that the rotation in Q below the
! merge is an identity. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  TOP             LOGICAL
!                    .TRUE.: fuse at top
!                    .FALSE.: fuse at bottom
!
!  N               INTEGER
!                    dimension of matrix
!
!  STR             INTEGER
!                    first non-identity rotation
!
!  STP             INTEGER
!                    last non-identity rotation
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  G               REAL(8) array of dimension (3)
!                    generators for rotation that will be fused
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkfact_mergebulge(TOP,N,STR,STP,P,Q,D,G)

  implicit none
  
  ! input variables
  logical, intent(in) :: TOP
  integer, intent(in) :: N, STR, STP
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)),D(2*N)
  real(8), intent(in) :: G(3)
  
  ! compute variables
  integer :: jj, ind, up, down
  real(8) :: c1r, c1i, s1
  real(8) :: c2r, c2i, s2
  real(8) :: c3r, c3i, s3r, s3i
  real(8) :: d1r, d1i
  real(8) :: d2r, d2i
  real(8) :: phr, phi
  real(8) :: nrm
  
  ! fusion at top from left
  if (TOP.AND..NOT.P(STR)) then
  
    ! set inputs  
    c2r = G(1)
    c2i = G(2)
    s2 = G(3)
  
    ! retrieve Q  
    c1r = Q(3*STR-2)
    c1i = Q(3*STR-1)
    s1 = Q(3*STR)

    ! compute product GQ    
    c3r = c1r*c2r - c1i*c2i - s1*s2
    c3i = c1r*c2i + c1i*c2r
    s3r = s1*c2r + s2*c1r
    s3i = -(s1*c2i - s2*c1i)
    
    ! compute phase
    call d_rot2_vec2gen(s3r,s3i,phr,phi,nrm)

    ! update Q
    c2r = c3r*phr + c3i*phi
    c2i = -c3r*phi + c3i*phr
    s2 = nrm
    
    call z_rot3_vec3gen(c2r,c2i,s2,Q(3*STR-2),Q(3*STR-1),Q(3*STR),nrm)
    
    ! update D
    d1r = D(2*STR-1)
    d1i = D(2*STR)
            
    nrm = phr*d1r - phi*d1i
    d1i = phr*d1i + phi*d1r
    d1r = nrm

    call d_rot2_vec2gen(d1r,d1i,D(2*STR-1),D(2*STR),nrm)

    ! initialize downward index
    down = STR
        
    ! move phase downward
    do jj = STR,(STP-1)
        
      ! exit if P == .TRUE.
      if (P(down)) then
        exit
      end if
                
      ! set downward index
      down = 1 + jj
        
      ! update Q
      c2r = Q(3*down-2)
      c2i = Q(3*down-1)
      s2 = Q(3*down)
                
      nrm = phr*c2r + phi*c2i
      c2i = phr*c2i - phi*c2r
      c2r = nrm

      call z_rot3_vec3gen(c2r,c2i,s2,Q(3*down-2),Q(3*down-1),Q(3*down),nrm)
        
    end do
      
    ! update downward index
    down = down + 1
      
    ! update downward diagonal
    d2r = D(2*down-1)
    d2i = D(2*down)
            
    nrm = phr*d2r + phi*d2i
    d2i = phr*d2i - phi*d2r
    d2r = nrm
      
    call d_rot2_vec2gen(d2r,d2i,D(2*down-1),D(2*down),nrm)

  ! fusion at top from right
  else if (TOP.AND.P(STR)) then
 
    ! set inputs  
    c1r = G(1)
    c1i = G(2)
    s1 = G(3)
  
    ! retrieve Q  
    c2r = Q(3*STR-2)
    c2i = Q(3*STR-1)
    s2 = Q(3*STR)

    ! compute product QG    
    c3r = c1r*c2r - c1i*c2i - s1*s2
    c3i = c1r*c2i + c1i*c2r
    s3r = s1*c2r + s2*c1r
    s3i = -(s1*c2i - s2*c1i)
    
    ! compute phase
    call d_rot2_vec2gen(s3r,s3i,phr,phi,nrm)

    ! update Q
    c2r = c3r*phr + c3i*phi
    c2i = -c3r*phi + c3i*phr
    s2 = nrm
    
    call z_rot3_vec3gen(c2r,c2i,s2,Q(3*STR-2),Q(3*STR-1),Q(3*STR),nrm)
    
    ! update D
    d1r = D(2*STR-1)
    d1i = D(2*STR)
            
    nrm = phr*d1r - phi*d1i
    d1i = phr*d1i + phi*d1r
    d1r = nrm

    call d_rot2_vec2gen(d1r,d1i,D(2*STR-1),D(2*STR),nrm)

    ! update downward diagonal
    d2r = D(2*STR+1)
    d2i = D(2*STR+2)
            
    nrm = phr*d2r + phi*d2i
    d2i = phr*d2i - phi*d2r
    d2r = nrm
      
    call d_rot2_vec2gen(d2r,d2i,D(2*STR+1),D(2*STR+2),nrm)
 
  ! fusion at bottom from right
  else if (.NOT.TOP.AND..NOT.P(STP-1)) then
 
    ! set inputs  
    c1r = G(1)
    c1i = G(2)
    s1 = G(3)
  
    ! retrieve Q  
    c2r = Q(3*STP-2)
    c2i = Q(3*STP-1)
    s2 = Q(3*STP)

    ! compute product QG
    c3r = c1r*c2r - c1i*c2i - s1*s2
    c3i = c1r*c2i + c1i*c2r
    s3r = s1*c2r + s2*c1r
    s3i = -(s1*c2i - s2*c1i)
   
    ! compute phase
    call d_rot2_vec2gen(s3r,s3i,phr,phi,nrm)

    ! update Q
    c2r = c3r*phr + c3i*phi
    c2i = -c3r*phi + c3i*phr
    s2 = nrm
    
    call z_rot3_vec3gen(c2r,c2i,s2,Q(3*STP-2),Q(3*STP-1),Q(3*STP),nrm)
    
    ! update D
    d1r = D(2*STP-1)
    d1i = D(2*STP)
            
    nrm = phr*d1r - phi*d1i
    d1i = phr*d1i + phi*d1r
    d1r = nrm

    call d_rot2_vec2gen(d1r,d1i,D(2*STP-1),D(2*STP),nrm)

    ! update downward diagonal
    d2r = D(2*STP+1)
    d2i = D(2*STP+2)
            
    nrm = phr*d2r + phi*d2i
    d2i = phr*d2i - phi*d2r
    d2r = nrm
 
    call d_rot2_vec2gen(d2r,d2i,D(2*STP+1),D(2*STP+2),nrm)
  
  ! fusion at bottom from left
  else
  
    ! set inputs  
    c1r = G(1)
    c1i = G(2)
    s1 = G(3)
  
    ! retrieve Q  
    ind = 3*(STP-1)
    c2r = Q(ind+1)
    c2i = Q(ind+2)
    s2 = Q(ind+3)
    
    ! compute product QG
    c3r = c1r*c2r - c1i*c2i - s1*s2
    c3i = c1r*c2i + c1i*c2r
    s3r = s1*c2r + s2*c1r
    s3i = s1*c2i - s2*c1i
    
    ! compute phase
    call d_rot2_vec2gen(s3r,s3i,phr,phi,nrm)

    ! update Q
    c2r = c3r*phr + c3i*phi
    c2i = -c3r*phi + c3i*phr
    s2 = nrm
    
    call z_rot3_vec3gen(c2r,c2i,s2,Q(ind+1),Q(ind+2),Q(ind+3),nrm)
    
    ! retrieve D
    ind = 2*STP
    d2r = D(ind+1)
    d2i = D(ind+2)
     
    ! update lower diagonal 
    c2r = phr*d2r + phi*d2i
    c2i = phr*d2i - phi*d2r

    call d_rot2_vec2gen(c2r,c2i,D(ind+1),D(ind+2),nrm)
  
    ! initialize upward index
    up = STP
        
    ! move phase upward
    do jj = STR,(STP-1)
        
      ! exit if P == .FALSE.
      if (.NOT.P(up-1)) then
        exit
      end if
                
      ! set upward index
      up = up - 1
        
      ! update Q
      c2r = Q(3*up-2)
      c2i = Q(3*up-1)
      s2 = Q(3*up)
                
      nrm = phr*c2r + phi*c2i
      c2i = phr*c2i - phi*c2r
      c2r = nrm

      call z_rot3_vec3gen(c2r,c2i,s2,Q(3*up-2),Q(3*up-1),Q(3*up),nrm)
        
    end do
      
    ! update upward diagonal
    d1r = D(2*up-1)
    d1i = D(2*up)
            
    nrm = phr*d1r - phi*d1i
    d1i = phr*d1i + phi*d1r
    d1r = nrm
      
    call d_rot2_vec2gen(d1r,d1i,D(2*up-1),D(2*up),nrm)

  end if

end subroutine z_uprkfact_mergebulge
