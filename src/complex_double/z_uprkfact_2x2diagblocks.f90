#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkfact_2x2diagblocks
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a set of two by two diagonal blocks of a 
! unitary plus rank one matrix pencil stored in factored form.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  TOP             LOGICAL
!                    .TRUE.: top block is computed
!                    .FALSE.: bottom block is computed
!
!  HESS            LOGICAL
!                    .TRUE.: the extended hessenberg part is included
!                    .FALSE.: only upper triangular parts are returned
!
!  QZ              LOGICAL
!                    .TRUE.: second triangular factor is assumed nonzero
!                    .FALSE.: second triangular factor is assumed to be identity
!
!  N               INTEGER
!                    dimension of matrix
!
!  K               INTEGER
!                    number of triangulars/rank
!
!  ROW             INTEGER
!                    first row of 2x2 block, it is assumed that in case
!                    of TOP that there is no Q-rotation above
!                    and in case of not TOP that there is no Q-rotation below
!
!  P               LOGICAL
!                    position flag for Q
!
!  Q               REAL(8) array of dimension (6)
!                    array of generators for first sequence of rotations
!                    if HESS = .FALSE., unused
!
!  D1,D2           REAL(8) arrays of dimension (4)
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!                    if QZ = .FALSE., D2 is unused
!
!  C1,B1,C2,B2     REAL(8) arrays of dimension (6)
!                    array of generators for upper-triangular parts of the pencil
!                    if QZ = .FALSE., C2 and B2 are unused
!
! OUTPUT VARIABLES:
!
!  A,B             COMPLEX(8) array of dimension (2,2)
!                    on exit contains the desired 2x2 block from
!                    the extended hessenberg matrix 
!                    if QZ = .FALSE., B is unused
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkfact_2x2diagblocks(TOP,HESS,QZ,N,K,ROW,P,Q,D1,C1,B1,D2,C2,B2,A,B)
  
  implicit none
  
  ! input variables
  logical, intent(in) :: TOP, HESS, QZ, P(N-2)
  integer, intent(in) :: N, K, ROW
  real(8), intent(inout) :: Q(3*K*(N-1)), D1(2*K*(N+1)), C1(3*K*N), B1(3*K*N)
  real(8), intent(inout) :: D2(2*K*(N+1)), C2(3*K*N), B2(3*K*N)
  complex(8), intent(inout) :: A(2,2), B(2,2)
  
  ! compute variables
  complex(8) :: AA(2,2),BB(2,2)
  integer :: ir1,ir2,id1,id2, ir3,ir4,id3,id4, ii

  ir2 = 3*ROW+3; ir1 = ir2-5
  id2 = 2*ROW+2; id1 = id2-3

  call z_upr1fact_2x2diagblocks(TOP,HESS,QZ,P(ROW-1),Q((ir1-3):(ir2-3)),&
       &D1(id1:id2),C1(ir1:ir2),B1(ir1:ir2),&
       &D2(id1:id2),C2(ir1:ir2),B2(ir1:ir2),A,B)

  do ii=2,K
     ir4 = 3*(ii-1)*N+3*ROW+3; ir3 = ir4-5
     id4 = 2*(ii-1)*(N+1)+2*ROW+2; id3 = id4-3
     
     call z_upr1fact_2x2diagblocks(TOP,.FALSE.,QZ,P(ROW-1),Q((ir1-3):(ir2-3)),&
          &D1(id3:id4),C1(ir3:ir4),B1(ir3:ir4),&
          &D2(id3:id4),C2(ir3:ir4),B2(ir3:ir4),AA,BB)
     
     A = matmul(A,AA)
     
     if (QZ) then
        B = matmul(B,BB)
     end if
  end do
  
end subroutine z_uprkfact_2x2diagblocks
