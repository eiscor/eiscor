#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_buildbulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the first transformation for a sinlge shift
! iteration on a upr1 pencil.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  P               LOGICAL 
!                    position flag for Q
!
!  Q               REAL(8) array of dimension (3)
!                    array of generators for first sequence of rotations
!
!  D               REAL(8) arrays of dimension (4)
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C,B             REAL(8) arrays of dimension (6)
!                    array of generators for upper-triangular parts
!                    of the pencil
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
  real(8), intent(in) :: Q(3), D(4), C(6), B(6)
  complex(8), intent(in) :: SHFT
  real(8), intent(inout) :: G(3)
  
  ! compute variables
  real(8) :: nrm
  complex(8) :: R1(2,2), R2(2,2), vec1(2), vec2(2)
  
  ! get R1
  call z_upr1utri_decompress(.FALSE.,2,D,C,B,R1)
     
  ! set R2 to I
  R2 = cmplx(0d0,0d0,kind=8)
  R2(1,1) = cmplx(1d0,0d0,kind=8); R2(2,2) = R2(1,1)

  ! compute first columns
  ! P == FALSE
  if (.NOT.P) then
   
    ! first column of R1
    vec1(1) = cmplx(Q(1),Q(2),kind=8)*R1(1,1)
    vec1(2) = cmplx(Q(3),0d0,kind=8)*R1(1,1)
    
    ! first column of R2
    vec2(1) = R2(1,1)
    vec2(2) = R2(2,1)
  
  ! P == TRUE
  else
    
    ! Q*e1
    vec2(1) = cmplx(Q(1),-Q(2),kind=8)
    vec2(2) = cmplx(-Q(3),0d0,kind=8)
    
    ! back solve with R1
    vec2(2) = vec2(2)/R1(2,2)
    vec2(1) = (vec2(1) - R1(1,2)*vec2(2))/R1(1,1)
    
    ! R2^-1 e1
    vec1(1) = 1d0/R2(1,1)
    vec1(2) = 0d0
  
  end if
  
  ! insert shift
  vec1 = vec1 - SHFT*vec2
  
  ! compute eliminator
  call z_rot3_vec4gen(dble(vec1(1)),aimag(vec1(1)),dble(vec1(2)),&
                      aimag(vec1(2)),G(1),G(2),G(3),nrm)
      
end subroutine z_upr1fact_buildbulge
