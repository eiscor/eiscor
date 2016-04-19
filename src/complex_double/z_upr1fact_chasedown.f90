#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_chasedown 
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
!  P               LOGICAL array of dimension (2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (6)
!                    array of generators for givens rotations
!
!  D               REAL(8) arrays of dimension (4)
!                    array of generators for complex diagonal matrix
!
!  C,B             REAL(8) arrays of dimension (6)
!                    array of generators for upper-triangular part
!
!  M               INTEGER
!                    leading dimesnion of V and W
!
!  V               COMPLEX(8) array of dimension (M,2)
!                    right schur vectors
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. update V to store right schurvectors 
!
!  G               REAL(8) array of dimension (3)
!                    array of generators for misfit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_chasedown(VEC,P,Q,D,C,B,M,V,G)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: M
  logical, intent(inout) :: P(2)
  real(8), intent(inout) :: Q(6), D(4), C(6), B(6), G(3)
  complex(8), intent(inout) :: V(M,2)
  
  ! compute variables
!  integer :: 
  real(8) :: G1(3), G2(3), G3(3)
  complex(8) :: A(2,2) 
  
  ! set cores for turnover based on P(1)
  ! hess 
  if (.NOT.P(1)) then
 
    G1 = Q(1:3)
    G2 = Q(4:6)
    G3 = G

  ! invhess
  else

    G1 = G
    G2 = Q(1:3)
    G3 = Q(4:6)

  end if

  ! compute turnover
  call z_rot3_turnover(G1,G2,G3)

  ! bottom fusion based on FLAG
  ! hess
  if (.NOT.P(2)) then
 
    ! update Q
    Q(1:3) = G1
    Q(4:6) = G2

    ! update V
    if (VEC) then
     
      A(1,1) = cmplx(G3(1),G3(2),kind=8)
      A(2,1) = cmplx(G3(3),0d0,kind=8)
      A(1,2) = cmplx(-G3(3),0d0,kind=8)
      A(2,2) = cmplx(G3(1),-G3(2),kind=8)

      V = matmul(V,A)

    end if

    ! pass G3 through upper-triangular part
    call z_upr1utri_rot3swap(.FALSE.,D,C,B,G3)

    ! fuse G3 with Q
    call z_rot3_fusion(.TRUE.,Q(4:6),G3)

    ! scale rows of upper-triangular part
    call z_upr1utri_unimodscale(.TRUE.,D(1:2),C(1:3), & 
                                B(1:3),cmplx(G3(1),G3(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D(3:4),C(4:6), & 
                                B(4:6),cmplx(G3(1),-G3(2),kind=8))

  ! invhess
  else

    ! update Q
    Q(1:3) = G1
    Q(4:6) = G3

    ! pass G2 through upper-triangular part
    call z_upr1utri_rot3swap(.TRUE.,D,C,B,G2)

    ! update V using G2inv
    if (VEC) then
     
      A(1,1) = cmplx(G2(1),-G2(2),kind=8)
      A(2,1) = cmplx(-G2(3),0d0,kind=8)
      A(1,2) = cmplx(G2(3),0d0,kind=8)
      A(2,2) = cmplx(G2(1),G2(2),kind=8)

      V = matmul(V,A)

    end if

    ! fuse G2 with Q
    call z_rot3_fusion(.FALSE.,G2,Q(4:6))

    ! move G2 to the otherside and update V
    if (VEC) then
     
      V(:,1) = V(:,1)*cmplx(G2(1),G2(2),kind=8)
      V(:,2) = V(:,2)*cmplx(G2(1),-G2(2),kind=8)

    end if

    ! scale columns of upper-triangular part
    call z_upr1utri_unimodscale(.FALSE.,D(1:2), &
                                C(1:3), & 
                                B(1:3), &
                                cmplx(G2(1),G2(2),kind=8))
    call z_upr1utri_unimodscale(.FALSE.,D(3:4),C(4:6), & 
                                B(4:6),cmplx(G2(1),-G2(2),kind=8))

  end if

  ! update position flag
  P(1) = P(2)

end subroutine z_upr1fact_chasedown
