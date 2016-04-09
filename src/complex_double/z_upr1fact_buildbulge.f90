#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_buildbulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the first transformation for a sinlge shift
! iteration on a upr1 matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  P               LOGICAL 
!                    position flag for Q
!
!  Q               REAL(8) array of dimension (6)
!                    array of generators for first sequence of rotations
!
!  D               REAL(8) arrays of dimension (4)
!                    array of generators for complex diagonal matrix
!                    in the upper-triangular factor
!
!  C,B             REAL(8) arrays of dimension (6)
!                    array of generators for upper-triangular parts
!
!  SHFT            COMPLEX(8) 
!                    contains the shift needed for the first transformation
!
! OUTPUT VARIABLES:
!
!  G               REAL(8) array of dimension (3)
!                    on exit contains the generators for the first
!                    transformation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_buildbulge(P,Q,D,C,B,SHFT,G)
  
  implicit none
  
  ! input variables
  logical, intent(in) :: P
  real(8), intent(in) :: Q(6), D(4), C(6), B(6)
  complex(8), intent(in) :: SHFT
  real(8), intent(inout) :: G(3)
  
  ! compute variables
  real(8) :: nrm
  complex(8) :: A(2,2), vec1(2), vec2(2)
  
  ! get leading 2x2 block of upper-triangular part
  call z_upr1fact_extracttri(.FALSE.,2,D,C,B,A)
     
  ! compute first column
  ! hessenberg case
  if (.NOT.P) then
   
    ! first column of A
    vec1(1) = cmplx(Q(1),Q(2),kind=8)*A(1,1)
    vec1(2) = cmplx(Q(3),0d0,kind=8)*A(1,1)
    
    ! vec2
    vec2(1) = cmplx(1d0,0d0,kind=8)
    vec2(2) = cmplx(0d0,0d0,kind=8)
    
  ! inverse hessenberg case
  else
    
    ! Q*e1
    vec2(1) = cmplx(Q(1),-Q(2),kind=8)
    vec2(2) = cmplx(-Q(3),0d0,kind=8)
    
    ! back solve with A
    vec2(2) = vec2(2)/A(2,2)
    vec2(1) = (vec2(1) - A(1,2)*vec2(2))/A(1,1)
    
    ! vec1
    vec1(1) = cmplx(1d0,0d0,kind=8)
    vec1(2) = cmplx(0d0,0d0,kind=8)
  
  end if
  
  ! insert shift
  vec1 = vec1 - SHFT*vec2
  
  ! compute eliminator
  call z_rot3_vec4gen(dble(vec1(1)),aimag(vec1(1)),dble(vec1(2)),&
                      aimag(vec1(2)),G(1),G(2),G(3),nrm)
      
end subroutine z_upr1fact_buildbulge
