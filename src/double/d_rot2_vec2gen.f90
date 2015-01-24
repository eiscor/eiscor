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
! The rot2 refers to a rotation desribed by two double and the vec2 to a vector
! of lenght two described by two double.
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
!     0d0 !     0d0 !     1d0 !     0d0 |   0d0
! +/- INF | +/- XdX | +/- 1d0 |     0d0 |   INF
! +/- XdX | +/- INF |     0d0 | +/- 1d0 |   INF
! +/- INF | +/- INF |     NAN |     NAN |   NAN
!     NAN | +/- XdX |     NAN |     NAN |   NAN
! +/- XdX |     NAN |     NAN |     NAN |   NAN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_rot2_vec2gen(A,B,C,S,NRM)

  implicit none
  
  ! input variables
  real(8), intent(in) :: A,B
  real(8), intent(inout) :: C,S,NRM
  
  ! compute variables
  real(8) :: tB
  
  ! call DROTG from oenblas
  NRM = A; tB = B
  call drotg(NRM,tB,C,S)  
           
end subroutine d_rot2_vec2gen
