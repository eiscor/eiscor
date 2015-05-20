#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot4_turnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine uses the generators for three Givens rotations represented
! by 3 real numbers: the real and imaginary parts of a complex cosine
! and a scrictly real sine and performs a turnover. The input arrays should be 
! ordered:
!
!    G1  G3
!      G2
!
! The new generators are stored as:
!
!      G1
!    G3  G2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  G1, G2, G3       REAL(8) arrays of dimension (3)
!                    generators for givens rotations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot4_turnover(G1,G2,G3)

  implicit none

  ! input variables
  real(8), intent(inout) :: G1(4), G2(4), G3(3)

  ! compute variables
  real(8) :: nrm, Hr(3,2), Hi(3,2)

  real(8) :: c1r, c1i, s1r, s1i
  real(8) :: c2r, c2i, s2r, s2i
  real(8) :: c3r, c3i, s3

  real(8) :: c4r, c4i, s4
  real(8) :: c5r, c5i, s5r, s5i
  real(8) :: c6r, c6i, s6r, s6i

  ! set inputs
  c1r = G1(1)
  c1i = G1(2)
  s1r = G1(3)
  s1i = G1(4)

  c2r = G2(1)
  c2i = G2(2)
  s2r = G2(3)
  s2i = G2(4)

  c3r = G3(1)
  c3i = G3(2)
  s3 = G3(3)
 
  ! first two columns
  Hr(1,1) = c1r*c3r - c1i*c3i - c2i*s3*s1i - c2r*s3*s1r
  Hr(2,1) = c3r*s1r - c3i*s1i + c1i*c2i*s3 + c1r*c2r*s3
  Hr(3,1) = s3*s2r

  Hi(1,1) = c1i*c3r + c3i*c1r - c2i*s3*s1r + c2r*s3*s1i
  Hi(2,1) = c3i*s1r + c3r*s1i - c1i*c2r*s3 + c2i*c1r*s3
  Hi(3,1) = s3*s2i

  Hr(1,2) = -c1r*s3 - c3i*(c2i*s1r - c2r*s1i) - c3r*(c2i*s1i + c2r*s1r)
  Hr(2,2) = c3r*(c1i*c2i + c1r*c2r) - c3i*(c1i*c2r - c2i*c1r) - s3*s1r
  Hr(3,2) = c3i*s2i + c3r*s2r

  Hi(1,2) = c3i*(c2i*s1i + c2r*s1r) - c1i*s3 - c3r*(c2i*s1r - c2r*s1i)
  Hi(2,2) = -s3*s1i - c3i*(c1i*c2i + c1r*c2r) - c3r*(c1i*c2r - c2i*c1r)
  Hi(3,2) = c3r*s2i - c3i*s2r

  ! first rotation
  call z_rot3_vec4gen(Hr(2,1),Hi(2,1),Hr(3,1),Hi(3,1),c4r,c4i,s4,nrm)

  ! update H
  nrm = Hi(2,1)*c4i + Hr(2,1)*c4r + Hr(3,1)*s4
  Hi(2,1) = Hi(2,1)*c4r - Hr(2,1)*c4i + Hi(3,1)*s4
  Hr(2,1) = nrm

  c6r = Hi(2,2)*c4i + Hr(2,2)*c4r + Hr(3,2)*s4
  s6r = Hr(3,2)*c4r - Hi(3,2)*c4i - Hr(2,2)*s4
  c6i = Hi(2,2)*c4r - Hr(2,2)*c4i + Hi(3,2)*s4
  s6i = Hi(3,2)*c4r + Hr(3,2)*c4i - Hi(2,2)*s4
  Hr(2,2) = c6r
  Hi(2,2) = c6i
  Hr(3,2) = s6r
  Hi(3,2) = s6i

  ! second rotation
  call z_rot4_vec4gen(Hr(1,1),Hi(1,1),Hr(2,1),Hi(2,1),c5r,c5i,s5r,s5i,nrm)
  
  ! update H
  nrm = Hr(2,2)*c5r - Hi(2,2)*c5i + Hi(1,2)*s5i - Hr(1,2)*s5r
  Hi(2,2) = Hi(2,2)*c5r + Hr(2,2)*c5i - Hi(1,2)*s5r - Hr(1,2)*s5i
  Hr(2,2) = nrm

  ! third rotation
  call z_rot4_vec4gen(Hr(2,2),Hi(2,2),Hr(3,2),Hi(3,2),c6r,c6i,s6r,s6i,nrm)
  
  ! set outputs
  G1(1) = c5r
  G1(2) = c5i
  G1(3) = s5r
  G1(4) = s5i

  G2(1) = c6r
  G2(2) = c6i
  G2(3) = s6r
  G2(4) = s6i

  G3(1) = c4r
  G3(2) = c4i
  G3(3) = s4
  
end subroutine z_rot4_turnover
