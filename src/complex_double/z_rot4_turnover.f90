#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot4_turnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine uses the generators for three Givens rotations. Two of these
! rotations are represented by 4 real numbers: the real and imaginary parts 
! of a complex cosine and sine and the third is represented  by 3 real 
! numbers: the real and imaginary parts of a complex cosine and a scrictly 
! real sine and performs a turnover. The input arrays should be 
! ordered:
!
!    G1  G3
!      G2
!
! It is always assumed that G3 is the 3 parameter rotation. There is a flag 
! that dictates whether the 3 parameter rotation appears on the right or the
! left of the factorization after the turnover. If FLAG == 0 then the new 
! generators are stored as:
!
!      G1
!    G3  G2
!
! and if FLAG == 1 then the new generators are stored as:
!
!      G1
!    G2  G3
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  G1, G2, G3       REAL(8) arrays of dimension (3)
!                    generators for givens rotations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot4_turnover(flag,G1,G2,G3)

  implicit none

  ! input variables
  logical, intent(in) :: flag
  real(8), intent(inout) :: G1(4), G2(4), G3(3)

  ! compute variables
  real(8) :: nrm, C1real(3), C1imag(3), C2real(3), C2imag(3)

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

  ! flag == 0 case
  if (.NOT.flag) then
  
    ! compute first column of G1*G2*G3
    C1real(1) = c1r*c3r - c1i*c3i - c2i*s3*s1i - c2r*s3*s1r
    C1real(2) = c3r*s1r - c3i*s1i + c1i*c2i*s3 + c1r*c2r*s3
    C1real(3) = s3*s2r

    C1imag(1) = c1i*c3r + c3i*c1r - c2i*s3*s1r + c2r*s3*s1i
    C1imag(2) = c3i*s1r + c3r*s1i - c1i*c2r*s3 + c2i*c1r*s3
    C1imag(3) = s3*s2i

    ! compute second column of G1*G2*G3
    C2real(1) = -c1r*s3 - c3i*(c2i*s1r - c2r*s1i) - c3r*(c2i*s1i + c2r*s1r)
    C2real(2) = c3r*(c1i*c2i + c1r*c2r) - c3i*(c1i*c2r - c2i*c1r) - s3*s1r
    C2real(3) = c3i*s2i + c3r*s2r

    C2imag(1) = c3i*(c2i*s1i + c2r*s1r) - c1i*s3 - c3r*(c2i*s1r - c2r*s1i) 
    C2imag(2) = -s3*s1i - c3i*(c1i*c2i + c1r*c2r) - c3r*(c1i*c2r - c2i*c1r) 
    C2imag(3) = c3r*s2i - c3i*s2r
 
    ! compute 3 parameter rotation
    call z_rot3vec4gen(C1real(2),C1imag(2),C1real(3),C1imag(3),c4r,c4i,s4,nrm)

    ! update first column using c4r, c4i and s4
    nrm = Ci2_1*c4i + Cr2_1*c4r + Cr3_1*s4
    Ci2_1 = Ci2_1*c4r - Cr2_1*c4i + Ci3_1*s4
    Cr2_1 = nrm
    

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

  end if
  
end subroutine z_rot4_turnover
