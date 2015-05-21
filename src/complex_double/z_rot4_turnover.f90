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
  real(8) :: nrm, Hr(3,3), Hi(3,3), T(4)

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


!!$  ! v1
!!$  ! first column
!!$  Hr(1,1) = c1r*c3r - c1i*c3i + (-c2i*s1i - c2r*s1r)*s3
!!$  Hr(2,1) = c3r*s1r - c3i*s1i + (c1i*c2i + c1r*c2r)*s3
!!$  Hr(3,1) = s3*s2r
!!$
!!$  Hi(1,1) = c1i*c3r + c3i*c1r + (-c2i*s1r + c2r*s1i)*s3
!!$  Hi(2,1) = c3i*s1r + c3r*s1i + (-c1i*c2r + c2i*c1r)*s3
!!$  Hi(3,1) = s3*s2i
!!$
!!$  ! first rotation
!!$  call z_rot3_vec4gen(Hr(2,1),Hi(2,1),Hr(3,1),Hi(3,1),c4r,c4i,s4,nrm)
!!$
!!$  ! update first column of H
!!$  nrm = Hi(2,1)*c4i + Hr(2,1)*c4r + Hr(3,1)*s4
!!$  Hi(2,1) = Hi(2,1)*c4r - Hr(2,1)*c4i + Hi(3,1)*s4
!!$  Hr(2,1) = nrm
!!$
!!$  ! second rotation
!!$  call z_rot4_vec4gen(Hr(1,1),Hi(1,1),Hr(2,1),Hi(2,1),c5r,c5i,s5r,s5i,nrm)
!!$  
!!$  ! last column
!!$  Hr(1,3) = s1r*s2r - s1i*s2i
!!$  Hr(2,3) = c1i*s2i - c1r*s2r
!!$  Hr(3,3) = c2r
!!$
!!$  Hi(1,3) = -s1i*s2r - s2i*s1r
!!$  Hi(2,3) = c1i*s2r + c1r*s2i
!!$  Hi(3,3) = -c2i
!!$
!!$  ! update third column of H
!!$  c6r = Hi(2,3)*c4i + Hr(2,3)*c4r + Hr(3,3)*s4
!!$  s6r = Hr(3,3)*c4r - Hi(3,3)*c4i - Hr(2,3)*s4
!!$  c6i = Hi(2,3)*c4r - Hr(2,3)*c4i + Hi(3,3)*s4
!!$  s6i = Hi(3,3)*c4r + Hr(3,3)*c4i - Hi(2,3)*s4
!!$  Hr(2,3) = c6r
!!$  Hi(2,3) = c6i
!!$  Hr(3,3) = -s6r
!!$  Hi(3,3) = s6i
!!$
!!$  ! update third column of H
!!$  nrm = Hr(2,3)*c5r - Hi(2,3)*c5i + Hi(1,3)*s5i - Hr(1,3)*s5r
!!$  Hi(2,3) = Hi(2,3)*c5r + Hr(2,3)*c5i - Hi(1,3)*s5r - Hr(1,3)*s5i
!!$  Hr(2,3) = -nrm
!!$
!!$  ! third rotation
!!$  call z_rot4_vec4gen(Hr(3,3),Hi(3,3),Hr(2,3),Hi(2,3),c6r,c6i,s6r,s6i,nrm)
!!$  
!!$  print*, c4r,c4i,s4
!!$  print*, c5r,c5i,s5r,s5i
!!$  print*, c6r,c6i,s6r,s6i


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! v2
  T(1) = c3r*s1r - c3i*s1i + (c1i*c2i + c1r*c2r)*s3
  T(2) = c3i*s1r + c3r*s1i + (-c1i*c2r + c2i*c1r)*s3
  T(3) = s3*s2r
  T(4) = s3*s2i

  ! compute first rotation
  call z_rot3_vec4gen(T(1),T(2),T(3),T(4),c4r,c4i,s4,nrm)

  nrm = T(2)*c4i + T(1)*c4r + T(3)*s4
  T(4) = T(2)*c4r - T(1)*c4i + T(4)*s4
  T(3) = nrm
  T(1) = c1r*c3r - c1i*c3i - (s1r*c2r + c2i*s1i)*s3
  T(2) = c1r*c3i + c1i*c3r + (-s1r*c2i + c2r*s1i)*s3
  
  ! compute second rotation
  call z_rot4_vec4gen(T(1),T(2),T(3),T(4),c5r,c5i,s5r,s5i,nrm)
!!$
!!$  ! last column
!!$  Hr(1,3) = s1r*s2r - s1i*s2i
!!$  Hr(2,3) = c1i*s2i - c1r*s2r
!!$  !Hr(3,3) = c2r
!!$
!!$  Hi(1,3) = -s1i*s2r - s2i*s1r
!!$  Hi(2,3) = c1i*s2r + c1r*s2i
!!$  !Hi(3,3) = -c2i
!!$
!!$
!!$  c6r = Hi(2,3)*c4i + Hr(2,3)*c4r + c2r*s4
!!$  s6r = -c2r*c4r - c2i*c4i + Hr(2,3)*s4
!!$  c6i = Hi(2,3)*c4r - Hr(2,3)*c4i -c2i*s4
!!$  s6i = -c2i*c4r + c2r*c4i - Hi(2,3)*s4
!!$
!!$
!!$  !Hr(2,3) = c6r
!!$  !Hi(2,3) = c6i
!!$  !Hr(3,3) = s6r
!!$  !Hi(3,3) = s6i
!!$
!!$  ! update third column of H
!!$  !nrm = -Hr(2,3)*c5r + Hi(2,3)*c5i - Hi(1,3)*s5i + Hr(1,3)*s5r
!!$  !Hi(2,3) = Hi(2,3)*c5r + Hr(2,3)*c5i - Hi(1,3)*s5r - Hr(1,3)*s5i
!!$  !Hr(2,3) = nrm
!!$  
!!$  Hr(2,3) = -c6r*c5r + c6i*c5i - Hi(1,3)*s5i + Hr(1,3)*s5r
!!$  Hi(2,3) = c6i*c5r + c6r*c5i - Hi(1,3)*s5r - Hr(1,3)*s5i
!!$  !Hr(2,3) = nrm



!!$  c6r = (c1i*s2r + c1r*s2i)*c4i + (c1i*s2i - c1r*s2r)*c4r + c2r*s4
!!$  T(1) = -c2r*c4r - c2i*c4i + (c1i*s2i - c1r*s2r)*s4
!!$  c6i = (c1i*s2r + c1r*s2i)*c4r - (c1i*s2i - c1r*s2r)*c4i -c2i*s4
!!$  T(2) = -c2i*c4r + c2r*c4i - (c1i*s2r + c1r*s2i)*s4
!!$
!!$  T(3) = -c6r*c5r + c6i*c5i + (s1i*s2r + s2i*s1r)*s5i + (s1r*s2r - s1i*s2i)*s5r
!!$  T(4) = c6i*c5r + c6r*c5i + (s1i*s2r + s2i*s1r)*s5r - (s1r*s2r - s1i*s2i)*s5i


  !c6r = (c1i*s2r + c1r*s2i)*c4i + (c1i*s2i - c1r*s2r)*c4r + c2r*s4
  T(1) = c2r*c4r + c2i*c4i + (-c1i*s2i + c1r*s2r)*s4
  !c6i = (c1i*s2r + c1r*s2i)*c4r - (c1i*s2i - c1r*s2r)*c4i -c2i*s4
  T(2) = c2i*c4r - c2r*c4i + (c1i*s2r + c1r*s2i)*s4

  T(3) = -((c1i*s2r + c1r*s2i)*c4i + (c1i*s2i - c1r*s2r)*c4r + c2r*s4)*c5r &
       &+ ((c1i*s2r + c1r*s2i)*c4r + (-c1i*s2i + c1r*s2r)*c4i - c2i*s4)*c5i &
       &+ (s1i*s2r + s2i*s1r)*s5i + (s1r*s2r - s1i*s2i)*s5r
  T(4) = ((c1i*s2r + c1r*s2i)*c4r + (-c1i*s2i + c1r*s2r)*c4i - c2i*s4)*c5r &
       &+ ((c1i*s2r + c1r*s2i)*c4i + (c1i*s2i - c1r*s2r)*c4r + c2r*s4)*c5i &
       &+ (s1i*s2r + s2i*s1r)*s5r + (-s1r*s2r + s1i*s2i)*s5i

  nrm = T(1)*T(1) + T(2)*T(2) + T(3)*T(3) + T(4)*T(4)
  if (abs(nrm-1)<3d0*EISCOR_DBL_EPS) then
     c6r = T(1)
     c6i = T(2)
     s6r = T(3)
     s6i = T(4)
  else
     ! third rotation
     call z_rot4_vec4gen(T(1),T(2),T(3),T(4),c6r,c6i,s6r,s6i,nrm)
  end if

!!$  print*, c4r,c4i,s4  
!!$  print*, c5r,c5i,s5r,s5i
!!$  print*, c6r,c6i,s6r,s6i
!!$  print*, ""
  
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
