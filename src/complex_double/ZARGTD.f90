#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZARGTD (Zomplex Auxiliary Routine Givens Through Diagonal)
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
!  JOB             CHARACTER
!                    'L': pass rotation from Left of diagonal to right
!                    'R': pass rotation from Right of diagonal to left
!
!  D               REAL(8) array of dimension (4)
!                    array of generators for complex diagonal matrix
!
!  B               REAL(8) array of dimension (3)
!                    generator for a Givens rotation
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies JOB is invalid
!                    INFO = -2 implies D is invalid
!                    INFO = -3 implies B is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZARGTD(JOB,D,B,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: JOB
  integer, intent(inout) :: INFO
  real(8), intent(inout) :: D(4), B(3)
  
  ! compute variables
  integer :: ii
  real(8) :: c1r, c1i, s1
  real(8) :: d1r, d1i
  real(8) :: d2r, d2i
  real(8) :: nrm
  
  ! initialize INFO
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
  
    ! check JOB
    if ((JOB.NE.'L').AND.(JOB.NE.'R')) then
      INFO = -1
      call UARERR(__FILE__,__LINE__,"JOB must be 'L' or 'R'",INFO,INFO)
      return
    end if
  
    ! check D
    call DARACH1(4,D,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"D is invalid",INFO,-2)
      return
    end if
  
    ! check B
    call DARACH1(3,B,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"B is invalid",INFO,-3)
      return
    end if

  end if
  
  ! set inputs
  c1r = B(1)
  c1i = B(2)
  s1 = B(3)
  
  ! retrieve D
  d1r = D(1)
  d1i = D(2)
  d2r = D(3)
  d2i = D(4)  
  
  !if JOB == L
  if (JOB.EQ.'L')then
    ! pass through diagonal
    nrm = (d1r*d2r + d1i*d2i)*c1r - (-d1r*d2i + d1i*d2r)*c1i
    c1i = (d1r*d2r + d1i*d2i)*c1i + (-d1r*d2i + d1i*d2r)*c1r
    c1r = nrm
  
    ! renormalize
    nrm = sqrt(c1r*c1r + c1i*c1i + s1*s1) ! not a good normalization, is that a problem?
    c1r = c1r/nrm
    c1i = c1i/nrm
    s1 = s1/nrm  
  ! if JOB == R
  else
    ! pass through diagonal
    nrm = (d1r*d2r + d1i*d2i)*c1r - (-d1r*d2i + d1i*d2r)*c1i
    c1i = (d1r*d2r + d1i*d2i)*c1i + (-d1r*d2i + d1i*d2r)*c1r
    c1r = nrm
  
    ! renormalize
    nrm = sqrt(c1r*c1r + c1i*c1i + s1*s1) ! not a good normalization, is that a problem?
    c1r = c1r/nrm
    c1i = c1i/nrm
    s1 = s1/nrm
  end if
  
  ! set B
  B(1) = c1r
  B(2) = c1i
  B(3) = s1
  
  ! set D
  D(1) = d2r 
  D(2) = d2i
  D(3) = d1r
  D(4) = d1i
  
  
end subroutine ZARGTD
