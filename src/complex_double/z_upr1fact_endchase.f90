#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_endchase 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine initializes one iteration of Francis' singleshift 
! algorithm for a upr1 factorization. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE. update schurvectors
!                    .FALSE. no schurvectors
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) arrays of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  C,B             REAL(8) arrays of dimension (3*N)
!                    array of generators for upper-triangular part
!
!  M               INTEGER
!                    leading dimesnion of V and W
!
!  V               COMPLEX(8) array of dimension (M,N)
!                    right schur vectors
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. update V to store right schurvectors 
!
!  FLAG            LOGICAL                     
!                    position flag for merging the misfit at the bottom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_endchase(VEC,N,P,Q,D,C,B,M,V,G,FLAG)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, FLAG
  integer, intent(in) :: M, N
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N), C(3*N), B(3*N), G(3)
  complex(8), intent(inout) :: V(M,N)
  
  ! compute variables
!  integer :: 
!  real(8) :: 
!  complex(8) :: 
  
  ! if N > 2 we do final turnover based on P(N-2) and FLAG
  if (N.GT.2) then

print*,""
print*,"WARNING: This case is not yet supported!"
print*,""
stop

  ! if N == 2 we do a single fusion
  else

    ! fuse Q and G, G is now a diagonal rotation
    call z_rot3_fusion(.TRUE.,Q(1:3),G)

    ! G scales the rows of the upper-triangular part
    call z_upr1utri_unimodscale(.TRUE.,D(1:2),C(1:3),B(1:3), &
                                cmplx(G(1),G(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D(3:4),C(4:6),B(4:6), &
                                cmplx(G(1),-G(2),kind=8))

  end if
  
end subroutine z_upr1fact_endchase
