#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DARGTO (Double Auxiliary Routine Givens TurnOver)
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
!  Q1, Q2, Q3       REAL(8) arrays of dimension (2)
!                    generators for givens rotations
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies Q1 is invalid
!                    INFO = -2 implies Q2 is invalid
!                    INFO = -3 implies Q3 is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DARGTO(Q1,Q2,Q3,INFO)

  implicit none
  
  ! input variables
  integer, intent(inout) :: INFO
  real(8), intent(inout) :: Q1(2), Q2(2), Q3(2)

  ! compute variables
  real(8), parameter :: tol = epsilon(1d0)
  real(8) :: nrm
  real(8) :: a, b 
  real(8) :: c1, s1
  real(8) :: c2, s2
  real(8) :: c3, s3
  real(8) :: c4, s4
  real(8) :: c5, s5
  real(8) :: c6, s6
  
  ! initialize INFO
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
  
    ! check Q1
    call DARACH1(2,Q1,INFO)
    call UARERR(__FILE__,__LINE__,"Q1 is invalid",INFO,-1)
    if (INFO.NE.0) then 
      return 
    end if 
    
    ! check Q2
    call DARACH1(2,Q2,INFO)
    call UARERR(__FILE__,__LINE__,"Q2 is invalid",INFO,-2)
    if (INFO.NE.0) then 
      return 
    end if   
    
    ! check Q3
    call DARACH1(2,Q3,INFO)
    call UARERR(__FILE__,__LINE__,"Q3 is invalid",INFO,-3)
    if (INFO.NE.0) then 
      return 
    end if  
  
  end if

  ! set local variables
  c1 = Q1(1)
  s1 = Q1(2)
  c2 = Q2(1)
  s2 = Q2(2)
  c3 = Q3(1)
  s3 = Q3(2)
  
  ! initialize c4 and s4
  a = s1*c3 + c1*c2*s3
  b = s2*s3
  
  ! compute first rotation
  call DARCG22(a,b,c4,s4,nrm,INFO)

  ! check INFO in debug mode
  if (DEBUG) then
    call UARERR(__FILE__,__LINE__,"DARCG22 failed",INFO,INFO)
    if (INFO.NE.0) then 
      return 
    end if 
  end if

  ! initialize c5 and s5
  a = c1*c3 - s1*c2*s3
  b = nrm

  ! compute second rotation
  call DARCG22(a,b,c5,s5,nrm,INFO)

  ! check INFO in debug mode
  if (DEBUG) then
    call UARERR(__FILE__,__LINE__,"DARCG22 failed",INFO,INFO)
    if (INFO.NE.0) then 
      return 
    end if 
  end if

  ! second column
  a = -(-c1*s3-s1*c2*c3)*s5 + ((-s1*s3+c1*c2*c3)*c4 + s2*c3*s4)*c5
  b = -(-s1*s3+c1*c2*c3)*s4 + s2*c3*c4
  
  ! compute first rotation
  call DARCG22(a,b,c6,s6,nrm,INFO)

  ! check INFO in debug mode
  if (DEBUG) then
    call UARERR(__FILE__,__LINE__,"DARCG22 failed",INFO,INFO)
    if (INFO.NE.0) then 
      return 
    end if 
  end if

  ! set output
  Q1(1) = c5
  Q1(2) = s5
  Q2(1) = c6
  Q2(2) = s6
  Q3(1) = c4
  Q3(2) = s4
       
end subroutine DARGTO
