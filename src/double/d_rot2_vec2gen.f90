#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_rot2_vec2gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generator for a Givens rotation represented by 2
! real numbers: A strictly real cosine and a scrictly real sine.  The first
! column is constructed to be parallel with the real vector [A,B]^T. 
!
! The rot2 refers to two double outputs and the vec2 to two double inputs;
! excluding the norm.
! 
! This routine serves the purpose of eiscor. The routine is fast and accurate
! enough. However, the routine is not robust against under- and
! overflow. Further, the routine does not compute continous rotations. For an
! algorithm with these properties, see "On computing Givens rotations reliably
! and efficiently" by D. Bindel, J. Demmel, W. Kahan, and O. Marques.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! EXCEPTIONAL CASES
!
!    A    |    B    |    C    |    S    |   NRM 
! ------- | ------- | ------- | ------- | -------
!   INF   |   0d0   |   1d0   |   0d0   |   INF
!   INF   |   XdX   |   1d0   |   0d0   |   INF
!   INF   |  -XdX   |   1d0   |   0d0   |   INF
!  -INF   |   0d0   |  -1d0   |   0d0   |   INF
!  -INF   |   XdX   |  -1d0   |   0d0   |   INF
!  -INF   |  -XdX   |  -1d0   |   0d0   |   INF
! ------- | ------- | ------- | ------- | -------
!   0d0   |   INF   |   0d0   |   1d0   |   INF
!   XdX   |   INF   |   0d0   |   1d0   |   INF
!  -XdX   |   INF   |   0d0   |   1d0   |   INF
!   0d0   |  -INF   |   0d0   |  -1d0   |   INF
!   XdX   |  -INF   |   0d0   |  -1d0   |   INF
!  -XdX   |  -INF   |   0d0   |  -1d0   |   INF
! ------- | ------- | ------- | ------- | -------
!   INF   |   INF   |   NAN   |   NAN   |   NAN
!  -INF   |   INF   |   NAN   |   NAN   |   NAN
!   INF   |  -INF   |   NAN   |   NAN   |   NAN
!  -INF   |  -INF   |   NAN   |   NAN   |   NAN
! ------- | ------- | ------- | ------- | -------
!   NAN   | +-XdX   |   NAN   |   NAN   |   NAN
! +-XdX   |   NAN   |   NAN   |   NAN   |   NAN
!   NAN   |   NAN   |   NAN   |   NAN   |   NAN
!   NAN   |   0d0   |   NAN   |   NAN   |   NAN
!   0d0   |   NAN   |   NAN   |   NAN   |   NAN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_rot2_vec2gen(A,B,C,S,NRM)

  implicit none
  
  ! input variables
  real(8), intent(in) :: A,B
  real(8), intent(inout) :: C,S,NRM
  
  ! compute variables
  real(8), parameter :: tol = epsilon(1d0)
 
  ! construct rotation
  if (B == 0) then
     if (A.NE.A) then
        C = A
        S = A
        NRM = A
     else if (A<0) then 
        C = -1d0
        S = 0d0
        NRM = -A
     else
        C = 1d0
        S = 0d0
        NRM = A
     endif
  else if (A**2 + B**2 - 1 < 1d1*tol) then
     C = A
     S = B
     NRM = 1d0
  else if (abs(A) >= abs(B)) then
     S = B/A

     if (A<0) then
        NRM = -sqrt(1.d0 + S**2)
     else
        NRM = sqrt(1.d0 + S**2)
     end if
     C =  1.d0/NRM
     S =  S*C
     NRM =  A*NRM

!!$     NRM = sqrt(1.d0 + S**2)
!!$     if (A<0) then
!!$        C = -1.d0/NRM
!!$        S =  S*C
!!$        NRM = -A*NRM
!!$     else
!!$        C =  1.d0/NRM
!!$        S =  S*C
!!$        NRM =  A*NRM
!!$     end if

  else
     C = A/B;

     if (B<0) then
        NRM = -sqrt(1.d0 + C**2)
     else
        NRM = sqrt(1.d0 + C**2)
     end if
     S =  1.d0/NRM
     C =  C*S
     NRM =  B*NRM

!!$     NRM = sqrt(1.d0 + C**2)
!!$     if (B<0) then
!!$        S = -1.d0/NRM
!!$        C =  C*S
!!$        NRM = -B*NRM
!!$     else
!!$        S =  1.d0/NRM
!!$        C =  C*S
!!$        NRM =  B*NRM
!!$     end if

  end if
           
end subroutine d_rot2_vec2gen
