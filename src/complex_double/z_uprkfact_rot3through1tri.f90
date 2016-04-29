#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkfact_rot3through1tri
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine passes a rotation through one upper-triangular matrix 
! stored as the product of a diagonal matrix and 2 sequences of 
! rotations.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  DIR             LOGICAL
!                    .TRUE.: pass rotation from left to right
!                    .FALSE.: pass rotation from right to left
!
!  N               INTEGER
!                    size of the matrix
!
!  K               INTEGER
!                    rank, i.e., number of upper triangulars
!
!  KC              INTEGER
!                    number of the triangular through which to pass
!
!  D               REAL(8) array of dimension (2*K*(N+1))
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C               REAL(8) array of dimension (3*K*N)
!                    array of generators for upper-triangular parts
!
!  B               REAL(8) array of dimension (3*K*N)
!                    array of generators for upper-triangular parts
!
!  G               REAL(8) array of dimension (3)
!                    generator for rotation
! 
!  ROW             INTEGER
!                    first row of rotation G
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkfact_rot3through1tri(DIR,N,K,KC,D,C,B,G,ROW)

  implicit none
  
  ! input variables
  logical, intent(in) :: DIR
  integer, intent(in) :: N, K, KC, ROW
  real(8), intent(inout) :: D(2*K*(N+1)), C(3*K*N), B(3*K*N), G(3)
 

  ! compute variables
  integer :: ind0,ind1

  ind0 = (KC-1)*2*(N+1)+(ROW-1)*2
  ind1 = (KC-1)*3*N+(ROW-1)*3

  print*, KC, ROW

  call z_upr1fact_rot3throughtri(DIR,&
       &D((ind0+1):(ind0+4)),&
       &C((ind1+1):(ind1+6)),&
       &B((ind1+1):(ind1+6)),G)
  
end subroutine z_uprkfact_rot3through1tri
