#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_mergebulge 
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
!  TOP             LOGICAL
!                    .TRUE.: fuses at the Top from the left
!                    .FALSE.: fuses at the Bottom from the right
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  B               REAL(8) array of dimension (3)
!                    generators for rotation that will be fused
!                    if TOP = .TRUE. contains diagonal rotation
!                    needed to adjust the phase
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_mergebulge(TOP,N,Q,D,B)

  implicit none
  
  ! input variables
  logical, intent(in) :: TOP
  integer, intent(in) :: N
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N), B(3)
  
  ! compute variables
  integer :: k,ii
  real(8) :: c1r, c1i, s1
  real(8) :: c2r, c2i, s2
  real(8) :: c3r, c3i, s3r, s3i
  real(8) :: d1r, d1i
  real(8) :: d2r, d2i
  real(8) :: phr, phi
  real(8) :: nrm
 
  ! fusion at top
  if (TOP) then

    ! set inputs  
    c2r = B(1)
    c2i = B(2)
    s2 = B(3)
    
    ! retrieve Q  
    c1r = Q(1)
    c1i = Q(2)
    s1 = Q(3)
     
    ! compute givens product
    c3r = c1r*c2r - c1i*c2i - s1*s2
    c3i = c1r*c2i + c1i*c2r
    s3r = s1*c2r + s2*c1r
    s3i = -(s1*c2i - s2*c1i)
     
    ! compute phase
    call d_rot2_vec2gen(s3r,s3i,phr,phi,nrm)

    ! store in B
    B(1) = phr
    B(2) = phi
    B(3) = 0d0
! not correct from here down
    ! update Q
    c2r = c3r*phr + c3i*phi
    c2i = -c3r*phi + c3i*phr
    s2 = nrm
    call z_rot3_vec3gen(c2r,c2i,s2,Q(1),Q(2),Q(3),nrm)

    do ii = 2,(N-1)
     
       k = 3*(ii-1)
      c1r = Q(k+1)
      c1i = Q(k+2)
      s1 = Q(k+3)
      
      nrm = c1r*phr + c1i*phi
      c1i = -c1r*phi + c1i*phr
      c1r = nrm
      call z_rot3_vec3gen(c1r,c1i,s1,Q(k+1),Q(k+2),Q(k+3),nrm)
    
    end do
     
    ! retrieve D
    d1r = D(1)
    d1i = D(2)

    ! update D
    c1r = phr*d1r - phi*d1i
    c1i = phr*d1i + phi*d1r
    call d_rot2_vec2gen(c1r,c1i,D(1),D(2),nrm)
     
  ! fusion at bottom
  else
  
    ! pass through diag
    call z_rot3_swapdiag(.FALSE.,D((2*N-3):(2*N)),B)
    
    ! set inputs  
    c2r = B(1)
    c2i = B(2)
    s2 = B(3)
  
    ! retrieve Q  
    c1r = Q(3*N-5)
    c1i = Q(3*N-4)
    s1 = Q(3*N-3)
     
    ! compute givens product
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
    call z_rot3_vec3gen(c2r,c2i,s2,Q(3*N-5),Q(3*N-4),Q(3*N-3),nrm)    
     
    ! retrieve D
    d1r = D(2*N-3)
    d1i = D(2*N-2)

    ! update D
    c1r = phr*d1r - phi*d1i
    c1i = phr*d1i + phi*d1r
    call d_rot2_vec2gen(c1r,c1i,D(2*N-3),D(2*N-2),nrm)

  end if

  ! retrieve D
  d2r = D(2*N-1)
  d2i = D(2*N)
     
  ! update D
  c2r = phr*d2r + phi*d2i
  c2i = phr*d2i - phi*d2r
  call d_rot2_vec2gen(c2r,c2i,D(2*N-1),D(2*N),nrm)
     
end subroutine z_unifact_mergebulge
