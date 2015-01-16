#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot3_turnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine uses the generators for three Givens rotations represented
! by 3 real numbers: the real and imaginary parts of a complex cosine
! and a scrictly real sine and performs a turnover.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  G1, G2, G3       REAL(8) arrays of dimension (3)
!                    generators for givens rotations
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies G1 is invalid
!                    INFO = -2 implies G2 is invalid
!                    INFO = -3 implies G3 is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3_turnover(G1,G2,G3,INFO)

  implicit none

  ! input variables
  integer, intent(inout) :: INFO
  real(8), intent(inout) :: G1(3), G2(3), G3(3)

  ! compute variables
  real(8), parameter :: tol = epsilon(1d0)
  real(8) :: nrm, T(3)
  
  real(8) :: c1r, c1i, s1
  real(8) :: c2r, c2i, s2
  real(8) :: c3r, c3i, s3
  
  real(8) :: c4r, c4i, s4
  real(8) :: c5r, c5i, s5
  real(8) :: c6r, c6i, s6
  
  ! initialize INFO
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
  
    ! check G1
    call d_1Darray_check(3,G1,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"G2 is invalid",INFO,-1)
      return
    end if
    nrm = sqrt(G1(1)**2 + G1(2)**2 + G1(3)**2)
    if (abs(nrm-1d0) > 1d1*tol) then
      INFO = -1      
      call u_infocode_check(__FILE__,__LINE__,"Q1 is not orthogonal to working precision",INFO,INFO)
      return
    end if
  
    ! check G2
    call d_1Darray_check(3,G2,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"G2 is invalid",INFO,-2)
      return
    end if
    nrm = sqrt(G2(1)**2 + G2(2)**2 + G2(3)**2)
    if (abs(nrm-1d0) > 1d1*tol) then
      INFO = -2      
      call u_infocode_check(__FILE__,__LINE__,"Q1 is not orthogonal to working precision",INFO,INFO)
      return
    end if
  
    ! check G3
    call d_1Darray_check(3,G3,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"G3 is invalid",INFO,-3)
      return
    end if
    nrm = sqrt(G3(1)**2 + G3(2)**2 + G3(3)**2)
    if (abs(nrm-1d0) > 1d1*tol) then
      INFO = -3
      call u_infocode_check(__FILE__,__LINE__,"Q1 is not orthogonal to working precision",INFO,INFO)
      return
    end if

  end if
  
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

  if ((s1.EQ.0d0) .AND. (s3.EQ.0d0)) then
     ! the case s1=s3=0d0 is special
     ! using the procedure for the generic case results in a rotation with
     ! c=1 and s=0 and in a rotation with a not necessarily real sine
     
     ! initialize c4r, c4i and s4
     T(1) = c1r*c2r + c1i*c2i
     T(2) = -c1i*c2r + c1r*c2i
     T(3) = s2
     
     ! compute first rotation
     call z_rot3_vec3gen(T(1),T(2),T(3),c4r,c4i,s4,nrm,INFO)
     
     ! check INFO in debug mode
     if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec3gen failed",INFO,INFO)
        if (INFO.NE.0) then 
           return 
        end if
     end if
     
     ! initialize c5r, c5i and s5
     T(1) = c1r*c3r - c1i*c3i 
     T(2) = c1r*c3i + c1i*c3r 
     T(3) = 0d0
     
     ! compute second rotation
     call z_rot3_vec3gen(T(1),T(2),T(3),c5r,c5i,s5,nrm,INFO)
     
     
     ! check INFO in debug mode
     if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec3gen failed",INFO,INFO)
        if (INFO.NE.0) then 
           return 
        end if
     end if
     
     ! initialize c6r, c6i and s6
     T(1) = c2r*c4r + c2i*c4i + c1r*s2*s4
     T(2) = c2i*c4r - c2r*c4i + c1i*s2*s4
     T(3) = s1*s2*s5 - (c5r*c2r + c5i*c2i)*s4 + s2*(c1r*(c4r*c5r + c4i*c5i) + c1i*(c4r*c5i - c4i*c5r))
     
     ! compute third rotation
     call z_rot3_vec3gen(T(1),T(2),T(3),c6r,c6i,s6,nrm,INFO)
     
     
     ! check INFO in debug mode
     if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec3gen failed",INFO,INFO)
        if (INFO.NE.0) then 
           return 
        end if
     end if
  else
     ! initialize c4r, c4i and s4
     T(1) = s1*c3r + (c1r*c2r + c1i*c2i)*s3
     T(2) = s1*c3i + (-c1i*c2r + c1r*c2i)*s3
     T(3) = s2*s3
     
     ! compute first rotation
     call z_rot3_vec3gen(T(1),T(2),T(3),c4r,c4i,s4,nrm,INFO)
     
     ! check INFO in debug mode
     if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec3gen failed",INFO,INFO)
        if (INFO.NE.0) then 
           return 
        end if
     end if
     
     ! initialize c5r, c5i and s5
     T(1) = c1r*c3r - c1i*c3i - s1*c2r*s3
     T(2) = c1r*c3i + c1i*c3r - s1*c2i*s3
     T(3) = nrm
     
     ! compute second rotation
     call z_rot3_vec3gen(T(1),T(2),T(3),c5r,c5i,s5,nrm,INFO)
     
     
     ! check INFO in debug mode
     if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec3gen failed",INFO,INFO)
        if (INFO.NE.0) then 
           return 
        end if
     end if
     
     ! initialize c6r, c6i and s6
     T(1) = c2r*c4r + c2i*c4i + c1r*s2*s4
     T(2) = c2i*c4r - c2r*c4i + c1i*s2*s4
     T(3) = s1*s2*s5 - (c5r*c2r + c5i*c2i)*s4 + s2*(c1r*(c4r*c5r + c4i*c5i) + c1i*(c4r*c5i - c4i*c5r))
     
     ! compute third rotation
     call z_rot3_vec3gen(T(1),T(2),T(3),c6r,c6i,s6,nrm,INFO)
     
     
     ! check INFO in debug mode
     if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec3gen failed",INFO,INFO)
        if (INFO.NE.0) then 
           return 
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
  
end subroutine z_rot3_turnover
