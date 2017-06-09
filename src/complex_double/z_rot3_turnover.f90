#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot3_turnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  G1, G2, G3       REAL(8) arrays of dimension (3)
!                    generators for givens rotations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3_turnover(G1,G2,G3)

  implicit none

  ! input variables
  real(8), intent(inout) :: G1(3), G2(3), G3(3)

  ! compute variables
  real(8) :: nrm, T(3), D(4)
  
  real(8) :: c1r, c1i, s1
  real(8) :: c2r, c2i, s2
  real(8) :: c3r, c3i, s3
  
  real(8) :: c4r, c4i, s4
  real(8) :: c5r, c5i, s5
  real(8) :: c6r, c6i, s6
  
  ! set local variables
  c1r = G1(1)
  c1i = G1(2)
  s1 = G1(3)
  c2r = G2(1)
  c2i = G2(2)
  s2 = G2(3)
  c3r = G3(1)
  c3i = G3(2)
  s3 = G3(3)

  !print*, ""
  !print*, s1, s2, s3

  
  if ((s1.EQ.0d0) .AND. (s3.EQ.0d0)) then
     ! the case s1=s3=0d0 is special
     ! using the procedure for the generic case results in a rotation with
     ! c=1 and s=0 and in a rotation with a not necessarily real sine
     
     ! initialize c4r, c4i and s4
     T(1) = c1r*c2r + c1i*c2i
     T(2) = -c1i*c2r + c1r*c2i
     T(3) = s2
     
     ! compute first rotation
     call z_rot3_vec3gen(T(1),T(2),T(3),c4r,c4i,s4,nrm)
     
     ! initialize c5r, c5i and s5
     T(1) = c1r*c3r - c1i*c3i 
     T(2) = c1r*c3i + c1i*c3r 
     T(3) = 0d0
     
     ! compute second rotation
     call z_rot3_vec3gen(T(1),T(2),T(3),c5r,c5i,s5,nrm)
     
     ! initialize c6r, c6i and s6
     T(1) = c2r*c4r + c2i*c4i + c1r*s2*s4
     T(2) = c2i*c4r - c2r*c4i + c1i*s2*s4
     T(3) = s1*s2*s5 - (c5r*c2r + c5i*c2i)*s4 + s2*(c1r*(c4r*c5r + c4i*c5i) &
     + c1i*(c4r*c5i - c4i*c5r))
     
     ! compute third rotation
     call z_rot3_vec3gen(T(1),T(2),T(3),c6r,c6i,s6,nrm)
     
  else
     if (s1.EQ.0d0) then
        D(1) = 1d0
        D(2) = 0d0
        call d_rot2_vec2gen(c1r,c1i,D(3),D(4),nrm)
        T(1) = c2r
        T(2) = c2i
        T(3) = s2
        !print*, D
        call z_rot3_swapdiag(D,T)
        c4r = T(1)
        c4i = T(2)
        s4  = T(3)
        !print*, D

        c6r = D(1)
        c6i = D(2)
        s6 = 0d0
               
        D(1) = 1d0
        D(2) = 0d0
        call d_rot2_vec2gen(c1r,-c1i,D(3),D(4),nrm)
        T(1) = c3r
        T(2) = c3i
        T(3) = s3
        !print*, D
        call z_rot3_swapdiag(D,T)
        c5r = T(1)
        c5i = T(2)
        s5  = T(3)
        !print*, D

        !print*, "1", c1r, c1i, s1  
        !print*, "2", c2r, c2i, s2
        !print*, "3", c3r, c3i, s3

        !print*, "4", c4r, c4i, s4  
        !print*, "5", c5r, c5i, s5
        !print*, "6", c6r, c6i, s6

        
     else
        if (s3.EQ.0d0) then
           D(3) = 1d0
           D(4) = 0d0
           call d_rot2_vec2gen(c3r,-c3i,D(1),D(2),nrm)
           T(1) = c2r
           T(2) = c2i
           T(3) = s2
           !print*, D
           call z_rot3_swapdiag(D,T)
           c6r = T(1)
           c6i = T(2)
           s6  = T(3)
           !print*, D
           
           c4r = D(3)
           c4i = -D(4)
           s4 = 0d0
           
           D(3) = 1d0
           D(4) = 0d0
           call d_rot2_vec2gen(c3r,c3i,D(1),D(2),nrm)
           T(1) = c1r
           T(2) = c1i
           T(3) = s1
           !print*, D
           call z_rot3_swapdiag(D,T)
           c5r = T(1)
           c5i = T(2)
           s5  = T(3)
           !print*, D
           
           !print*, "b1", c1r, c1i, s1  
           !print*, "b2", c2r, c2i, s2
           !print*, "b3", c3r, c3i, s3
           
           !print*, "b4", c4r, c4i, s4  
           !print*, "b5", c5r, c5i, s5
           !print*, "b6", c6r, c6i, s6
                   
        else
           
           ! initialize c4r, c4i and s4
           T(1) = s1*c3r + (c1r*c2r + c1i*c2i)*s3
           T(2) = s1*c3i + (-c1i*c2r + c1r*c2i)*s3
           T(3) = s2*s3
           
           ! compute first rotation
           call z_rot3_vec3gen(T(1),T(2),T(3),c4r,c4i,s4,nrm)
           
           ! initialize c5r, c5i and s5
           T(1) = c1r*c3r - c1i*c3i - s1*c2r*s3
           T(2) = c1r*c3i + c1i*c3r - s1*c2i*s3
           T(3) = nrm
           
           ! compute second rotation
           call z_rot3_vec3gen(T(1),T(2),T(3),c5r,c5i,s5,nrm)
           
           ! initialize c6r, c6i and s6
           T(1) = c2r*c4r + c2i*c4i + c1r*s2*s4
           T(2) = c2i*c4r - c2r*c4i + c1i*s2*s4
           T(3) = s1*s2*s5 - (c5r*c2r + c5i*c2i)*s4 + s2*(c1r*(c4r*c5r + c4i*c5i) &
                + c1i*(c4r*c5i - c4i*c5r))
           
           ! compute third rotation
           call z_rot3_vec3gen(T(1),T(2),T(3),c6r,c6i,s6,nrm)
        end if
     end if
  end if
  
  ! store output
  G1(1) = c5r
  G1(2) = c5i
  G1(3) = s5
  G2(1) = c6r
  G2(2) = c6i
  G2(3) = s6
  G3(1) = c4r
  G3(2) = c4i
  G3(3) = s4

  
  
!!$  if ((abs(s4).LT.EISCOR_DBL_EPS).OR.(abs(s5).LT.EISCOR_DBL_EPS)&
!!$       &.OR.(abs(s6).LT.EISCOR_DBL_EPS)) then
!!$     print*, "almost diagonal rotation in turnover"
!!$     print*, s4, s5, s6
!!$     !stop
!!$  end if

  
end subroutine z_rot3_turnover
