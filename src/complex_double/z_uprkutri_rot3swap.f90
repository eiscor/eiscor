#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkutri_rot3swap
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine passes a rotation through abs(k2-k1) unitary plus rank one
! upper-triangular matrix (upr1utri) stored as the product of a diagonal 
! matrix and 2 sequences of rotations.
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
!  K1,K2           INTEGER
!                    thr rotation is passed through the triangulars
!                    K1, K1+1, K1+2, ..., K2 if DIR.EQ..TRUE
!                    K2, K2-1, K2-2, ..., K1 if DIR.EQ..FALSE
!
!  D               REAL(8) array of dimension (2*K*N)
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C               REAL(8) array of dimension (3*K*N)
!                    first array of generators for upper-triangular part
!
!  B               REAL(8) array of dimension (3*K*N)
!                    second array of generators for upper-triangular part
!
!  G               REAL(8) array of dimension (3)
!                    generator for rotation
! 
!  ROW             INTEGER
!                    first row of rotation G
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkutri_rot3swap(DIR,N,K,K1,K2,D,C,B,G,ROW)

  implicit none
  
  logical, intent(in) :: DIR
  integer, intent(in) :: N, K, K1, K2, ROW
  real(8), intent(inout) :: D(2*K*N), C(3*K*N), B(3*K*N), G(3)
 

  ! compute variables
  integer :: ind0,ind1, kk

  
  if (DIR) then
     do kk=k1,k2
  
        ind0 = (kk-1)*2*N+(ROW-1)*2
        ind1 = (kk-1)*3*N+(ROW-1)*3
        
        call z_upr1utri_rot3swap(DIR,&
             &D((ind0+1):(ind0+4)),&
             &C((ind1+1):(ind1+6)),&
             &B((ind1+1):(ind1+6)),G)
     end do
  else 
     do kk=k2,k1,-1
        
        ind0 = (kk-1)*2*N+(ROW-1)*2
        ind1 = (kk-1)*3*N+(ROW-1)*3
        
        call z_upr1utri_rot3swap(DIR,&
             &D((ind0+1):(ind0+4)),&
             &C((ind1+1):(ind1+6)),&
             &B((ind1+1):(ind1+6)),G)
     end do
  end if
  
end subroutine z_uprkutri_rot3swap
