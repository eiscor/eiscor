#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkfpen_buildbulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the first transformation for a single shift
! iteration on a factored unitary plus rank one (uprkfpen) matrix pencil.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  P               LOGICAL
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (6)
!                    array of generators for first sequence of rotations
!
!  D1,D2           REAL(8) arrays of dimension (4)
!                    arrays of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1,C2,B1,B2     REAL(8) arrays of dimension (6)
!                    arrays of generators for unitary plus rank one
!                    upper-trinagular matrices
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkfpen_buildbulge(N,K,STR,STP,P,Q,D1,C1,B1,D2,C2,B2,SHFT,G)
  
  implicit none
  
  ! input variables
  integer, intent(in) :: N, K, STR, STP
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*N*K), C1(3*N*K), B1(3*N*K)
  real(8), intent(inout) :: D2(2*N*K), C2(3*N*K), B2(3*N*K)
  complex(8), intent(in) :: SHFT
  real(8), intent(inout) :: G(3)
  
  ! compute variables
  real(8) :: nrm
  complex(8) :: R1(2,2), R2(2,2), vec1(2), vec2(2)
  integer :: ind

  !print*, "z_uprkfpen_buildbulge", shft, N, STR, STP
  
  R1 = cmplx(0d0,0d0,kind=8)
  R2 = cmplx(0d0,0d0,kind=8)
  ! extract R1
  call z_uprkutri_decompress(.FALSE.,N,K,STR,STR,D1,C1,B1,R1)
  
  ! extract R2
  call z_uprkutri_decompress(.FALSE.,N,K,STR,STR,D2,C2,B2,R2)
     
  ! compute first columns
  if (STR.EQ.STP) then
     ind = 3*STR-2
     ! first column of R1
     vec1(1) = cmplx(Q(ind),Q(ind+1),kind=8)*R1(1,1)
     vec1(2) = cmplx(Q(ind+2),0d0,kind=8)*R1(1,1)
     
     ! first column of R2
     vec2(1) = R2(1,1)
     vec2(2) = R2(2,1)
  ! P == FALSE
  elseif (.NOT.P(STR)) then
   
     ind = 3*STR-2
     ! first column of R1
     vec1(1) = cmplx(Q(ind),Q(ind+1),kind=8)*R1(1,1)
     vec1(2) = cmplx(Q(ind+2),0d0,kind=8)*R1(1,1)
     
     ! first column of R2
     vec2(1) = R2(1,1)
     vec2(2) = R2(2,1)
  
  ! P == TRUE
  else
    
     ind = 3*STR-2
     ! Q*e1
     vec2(1) = cmplx(Q(ind),-Q(ind+1),kind=8)
     vec2(2) = cmplx(-Q(ind+2),0d0,kind=8)
     
     ! back solve with R1
     vec2(2) = vec2(2)/R1(2,2)
     vec2(1) = (vec2(1) - R1(1,2)*vec2(2))/R1(1,1)
     
     ! multiply by R2
     vec2 = matmul(R2,vec2)
     
     ! R2^-1 e1
     vec1(1) = cmplx(1d0,0d0,kind=8)
     vec1(2) = cmplx(0d0,0d0,kind=8)
  
  end if
  
  ! insert shift
  vec1 = vec1 - SHFT*vec2

  ! compute eliminator
  call z_rot3_vec4gen(dble(vec1(1)),aimag(vec1(1)),dble(vec1(2)),&
                      aimag(vec1(2)),G(1),G(2),G(3),nrm)
  
end subroutine z_uprkfpen_buildbulge
