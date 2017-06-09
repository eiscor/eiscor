#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkfpen_endchase 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine finaalizes one iteration of Francis' singleshift 
! algorithm for a factored unitary plus rank one (uprkfpen) matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VECR            LOGICAL
!                    .TRUE.: compute right schurvectors (V)
!                    .FALSE.: no schurvectors
!
!  VECL            LOGICAL
!                    .TRUE.: compute left schurvectors (W)
!                    .FALSE.: no schurvectors
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D1,D2           REAL(8) arrays of dimension (2*N)
!                    arrays of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1,C2,B1,B2     REAL(8) arrays of dimension (3*N)
!                    arrays of generators for unitary plus rank one
!                    upper-trinagular matrices
!
!  M               INTEGER
!                    leading dimension of V and W
!
!  V,W             COMPLEX(8) array of dimension (M,2)
!                    right and left schurvectors 
!
!  FLAG            LOGICAL                     
!                    position flag for merging the misfit at the bottom
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkfpen_endchase(VECR,VECL,N,K,STR,STP,P,Q,D1,C1,B1,D2,C2,B2,M,V,W,G,FLAG)

  implicit none
  
  ! input variables
  logical, intent(in) :: VECR, VECL, FLAG
  integer, intent(in) :: M, N, K, STR, STP
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*N), C1(3*N), B1(3*N)
  real(8), intent(inout) :: D2(2*N), C2(3*N), B2(3*N), G(3)
  complex(8), intent(inout) :: V(M,2),W(M,2)
  
  ! compute variables
  real(8) :: G1(3), G2(3), G3(3)
  complex(8) :: A(2,2) 

  !print*, "z_uprkfpen_endchase start"
  
  ! if only one rotation in Q we do a single fusion
  if ( STP.EQ.STR ) then

    ! fuse Q and G, G is now a diagonal rotation
    call z_rot3_fusion(.TRUE.,Q(3*STR-2:3*STR),G)

    ! G scales the rows of the upper-triangular part
    call z_upr1utri_unimodscale(.TRUE.,D1(2*STR-1:2*STR),C1(3*STR-2:3*STR),B1(3*STR-2:3*STR), &
                                cmplx(G(1),G(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D1(2*STR+1:2*STR+2),C1(3*STR+1:3*STR+3),B1(3*STR+1:3*STR+3), &
                                cmplx(G(1),-G(2),kind=8))
    ! does not touch B1 or C1

    ! return
    return

  end if
  
  !print*, "z_uprkfpen_endchase part 2"
  
  ! set cores for turnover based on P(STP-1)
  ! hess 
  if (.NOT.P(STP-1)) then
 
    G1 = Q(3*(STP-1)-2:3*(STP-1))
    G2 = Q(3*STP-2:3*STP)
    G3 = G

  ! invhess
  else

    G1 = G
    G2 = Q(3*STP-2:3*STP)
    G3 = Q(3*(STP-1)-2:3*(STP-1))

  end if

  ! compute turnover
  call z_rot3_turnover(G1,G2,G3)

  ! bottom fusion based on FLAG
  ! hess
  if (.NOT.FLAG) then
 
    ! update Q
    Q(3*(STP-1)-2:3*(STP-1)) = G1
    Q(3*STP-2:3*STP) = G2

    ! update left schurvectors with G3
    if (VECL) then 
    
      A(1,1) = cmplx(G3(1),G3(2),kind=8)
      A(2,1) = cmplx(G3(3),0d0,kind=8)
      A(1,2) = cmplx(-G3(3),0d0,kind=8)
      A(2,2) = cmplx(G3(1),-G3(2),kind=8)

      W = matmul(W,A)
   
    end if

    ! invert G3
    G3(2) = -G3(2)
    G3(3) = -G3(3)
 
    ! pass G3 through R2
    !call z_upr1utri_rot3swap(.TRUE.,D2(2*N-3:2*N), &
    !                         C2(3*N-5:3*N),B2(3*N-5:3*N),G3)
    call z_uprkutri_rot3swap(.TRUE.,N,K,1,K,D2,C2,B2,G3,STP)

    ! invert G3
    G3(2) = -G3(2)
    G3(3) = -G3(3)
 
    ! update right schurvectors with G3
    if (VECR) then
    
      A(1,1) = cmplx(G3(1),G3(2),kind=8)
      A(2,1) = cmplx(G3(3),0d0,kind=8)
      A(1,2) = cmplx(-G3(3),0d0,kind=8)
      A(2,2) = cmplx(G3(1),-G3(2),kind=8)

      V = matmul(V,A)
   
    end if

    ! pass G3 through R1
    !call z_upr1utri_rot3swap(.FALSE.,D1(2*N-3:2*N), &
    !                         C1(3*N-5:3*N),B1(3*N-5:3*N),G3)
    call z_uprkutri_rot3swap(.FALSE.,N,K,1,K,D1,C1,B1,G3,STP)

    ! fuse G3 with Q
    call z_rot3_fusion(.TRUE.,Q(3*STP-2:3*STP),G3)

    ! scale rows of R1
    call z_upr1utri_unimodscale(.TRUE.,D1(2*STP-1:2*STP), &
                                C1(3*STP-2:3*STP), & 
                                B1(3*STP-2:3*STP), &
                                cmplx(G3(1),G3(2),kind=8))
    ! does not touch C1 and B1
    call z_upr1utri_unimodscale(.TRUE.,D1(2*STP+1:2*STP+2), &
                                C1(3*STP+1:3*(STP+1)), & 
                                B1(3*STP+1:3*(STP+1)), &
                                cmplx(G3(1),-G3(2),kind=8))
    ! does not touch C1 and B1

  ! invhess
  else

    ! update Q
    Q(3*(STP-1)-2:3*(STP-1)) = G1
    Q(3*STP-2:3*STP) = G3

    ! pass G2 through R1
    !call z_upr1utri_rot3swap(.TRUE.,D1(2*(N-1)-1:2*N), &
    !                         C1(3*(N-1)-2:3*N),B1(3*(N-1)-2:3*N),G2)
    call z_uprkutri_rot3swap(.TRUE.,N,K,1,K,D1,C1,B1,G2,STP)

    ! invert G2
    G2(2) = -G2(2)
    G2(3) = -G2(3)

    ! update right schurvectors using G2
    if (VECR) then
     
      A(1,1) = cmplx(G2(1),G2(2),kind=8)
      A(2,1) = cmplx(G2(3),0d0,kind=8)
      A(1,2) = cmplx(-G2(3),0d0,kind=8)
      A(2,2) = cmplx(G2(1),-G2(2),kind=8)

      V = matmul(V,A)

    end if

    ! pass G2 through R2
    !call z_upr1utri_rot3swap(.FALSE.,D2(2*(N-1)-1:2*N), &
    !                         C2(3*(N-1)-2:3*N),B2(3*(N-1)-2:3*N),G2)
    call z_uprkutri_rot3swap(.FALSE.,N,K,1,K,D2,C2,B2,G2,STP)

    ! update left schurvectors using G2
    if (VECL) then
     
      A(1,1) = cmplx(G2(1),G2(2),kind=8)
      A(2,1) = cmplx(G2(3),0d0,kind=8)
      A(1,2) = cmplx(-G2(3),0d0,kind=8)
      A(2,2) = cmplx(G2(1),-G2(2),kind=8)

      W = matmul(W,A)

    end if

    ! invert G2
    G2(2) = -G2(2)
    G2(3) = -G2(3)

    ! fuse G2 with Q, G2 is now diagonal
    call z_rot3_fusion(.FALSE.,G2,Q(3*STP-2:3*STP))

    ! update left schurvectors with G2
    if (VECL) then
     
      W(:,1) = W(:,1)*cmplx(G2(1),G2(2),kind=8)
      W(:,2) = W(:,2)*cmplx(G2(1),-G2(2),kind=8)

    end if

    ! scale rows of R2
    call z_upr1utri_unimodscale(.TRUE.,D2(2*STP-1:2*STP), &
                                C2(3*STP-2:3*STP), & 
                                B2(3*STP-2:3*STP), &
                                cmplx(G2(1),-G2(2),kind=8))
    ! does not touch C1 and B1
    call z_upr1utri_unimodscale(.TRUE.,D2(2*STP+1:2*STP+2), &
                                C2(3*STP+1:3*(STP+1)), & 
                                B2(3*STP+1:3*(STP+1)), &
                                cmplx(G2(1),G2(2),kind=8))
    ! does not touch C1 and B1

  end if

  ! update position flag
  if (STP.GT.2) then
    P(STP-2) = P(STP-1)
  end if

  ! update P
  P(STP-1) = FLAG

end subroutine z_uprkfpen_endchase
