#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot3_swapdiag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine passes the generator for a Givens rotation represented
! by 3 real numbers: the real and imaginary parts of a complex cosine
! and a scrictly real sine through a 2x2 complex diagonal matrix 
! represented by 4 real numbers.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  D               REAL(8) array of dimension (4)
!                    array of generators for complex diagonal matrix
!
!  G               REAL(8) array of dimension (3)
!                    generator for a Givens rotation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3_swapdiag(D,G)

  implicit none
  
  ! input variables
  real(8), intent(inout) :: D(4), G(3)
  
  ! compute variables
  real(8) :: c1r, c1i, s1
  real(8) :: d1r, d1i
  real(8) :: d2r, d2i
  real(8) :: nrm
  
  ! set inputs
  c1r = G(1)
  c1i = G(2)
  s1 = G(3)

  ! retrieve D
  d1r = D(1)
  d1i = D(2)
  d2r = D(3)
  d2i = D(4)  
     

  if (s1.EQ.0d0) then
     ! fuse diagonals

     nrm =  d1r*c1r - d1i*c1i
     d1i =  d1r*c1i + d1i*c1r
     d1r =  nrm
     
     ! renormalize
     call d_rot2_vec2gen(d1r,d1i,D(1),D(2),nrm) 

     nrm =  d2r*c1r + d2i*c1i
     d2i = -d2r*c1i + d2i*c1r
     d2r =  nrm
     
     ! renormalize
     call d_rot2_vec2gen(d2r,d2i,D(3),D(4),nrm)      

     G(1) = 1d0
     G(2) = 0d0
     G(3) = 0d0
     
  else
       
     ! pass through diagonal
     nrm = (d1r*d2r + d1i*d2i)*c1r - (-d1r*d2i + d1i*d2r)*c1i
     c1i = (d1r*d2r + d1i*d2i)*c1i + (-d1r*d2i + d1i*d2r)*c1r
     c1r = nrm
     
     ! renormalize
     call z_rot3_vec3gen(c1r,c1i,s1,G(1),G(2),G(3),nrm) 
    
     ! set D
     D(1) = d2r 
     D(2) = d2i
     D(3) = d1r
     D(4) = d1i

  end if
  
end subroutine z_rot3_swapdiag
