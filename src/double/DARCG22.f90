#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DARCG22 (Double Auxiliary Routine Compute Givens generators 2->2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generator for A Givens rotation represented
! by 2 real numbers: A strictly real cosine and A scrictly real sine. 
! The first column is constructed to be parallel with the real vector 
! [A,B]^T. The (2->2) refers to the 2 double inputs and 2 double outputs 
! (excluding the norm).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  A               REAL(8) 
!                    the first component of the real vector [A,B]^T
!
!  B               REAL(8) 
!                    the second component of the real vector [A,B]^T
!
! OUTPUT VARIABLES:
!
!  C               REAL(8)
!                    on exit contains the generator for the cosine
!                    component of the Givens rotation
!
!  S               REAL(8)
!                    on exit contains the generator for the sine
!                    component of the Givens rotation
!
!  NRM             REAL(8)
!                    on exit contains the norm of the vector [A,B]^T
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies A is invalid
!                    INFO = -2 implies B is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DARCG22(A,B,C,S,NRM,INFO)

  implicit none
  
  ! input variables
  integer, intent(inout) :: INFO
  real(8), intent(in) :: A,B
  real(8), intent(inout) :: C,S,NRM
  
  ! compute variables
  real(8), parameter :: tol = epsilon(1d0)

  ! initialize INFO
  INFO = 0

  ! check input in debug mode
  if (DEBUG) then   
  
    ! Check A
    call DARACH1(1,A,INFO)
    call UARERR(__FILE__,__LINE__,"A is invalid",INFO,-1)
    if (INFO.NE.0) then 
      return 
    end if   
    
    ! Check B
    call DARACH1(1,B,INFO)
    call UARERR(__FILE__,__LINE__,"B is invalid",INFO,-2)
    if (INFO.NE.0) then 
      return 
    end if 
    
  end if 
  
  ! construct rotation
  NRM = 1d0  
  if (B == 0) then
     if (A<0) then 
        C = -1d0
        S = 0d0
        NRM = -A
     else
        C = 1d0
        S = 0d0
        NRM = A
     endif
  else if (abs(A) >= abs(B)) then
     S = B/A
     NRM = sqrt(1.d0 + S**2)
     if (A<0) then
        C = -1.d0/NRM
        S =  S*C
        NRM = -A*NRM
     else
        C =  1.d0/NRM
        S =  S*C
        NRM =  A*NRM
     end if
  else
     C = A/B;
     NRM = sqrt(1.d0 + C**2)
     if (B<0) then
        S = -1.d0/NRM
        C =  C*S
        NRM = -B*NRM
     else
        S =  1.d0/NRM
        C =  C*S
        NRM =  B*NRM
     end if
  end if
           
end subroutine DARCG22
