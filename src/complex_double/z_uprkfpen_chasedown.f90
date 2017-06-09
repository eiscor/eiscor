#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkfpen_chasedown 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine chases the a single misfit core transformation down 
! one row in a factored unitary plus rank one (uprkfpen) matrix. 
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
!  P               LOGICAL array of dimension (2)
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
!  M               INTEGER
!                    leading dimension of V and W
!
!  V,W             COMPLEX(8) array of dimension (M,2)
!                    right schur vectors
!                    if VECR/L = .FALSE. unused
!                    if VECR/L = .TRUE. update V/W to store right/left schurvectors 
!
!  MISFIT          REAL(8) array of dimension (3)
!                    array of generators for misfit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkfpen_chasedown(VECR,VECL,N,K,ii,P,Q,D1,C1,B1,D2,C2,B2,M,V,W,MISFIT)

  implicit none
  
  ! input variables
  logical, intent(in) :: VECR, VECL
  integer, intent(in) :: M, N, K, ii
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*N*K), C1(3*N*K), B1(3*N*K)
  real(8), intent(inout) :: D2(2*N*K), C2(3*N*K), B2(3*N*K), MISFIT(3)
  complex(8), intent(inout) :: V(M,2), W(M,2)
  
  ! compute variables
  real(8) :: G1(3), G2(3), G3(3)
  complex(8) :: A(2,2) 
  integer ::  ind1, ind2

  ! compute indices
  ind1 = 3*ii-2
  ind2 = 3*ii+1

  
  ! set cores for turnover based on P(ii)
  ! hess 
  if (.NOT.P(ii)) then
 
    G1 = Q(ind1:ind1+2)
    G2 = Q(ind2:ind2+2)
    G3 = MISFIT

  ! invhess
  else

    G1 = MISFIT
    G2 = Q(ind2:ind2+2)
    G3 = Q(ind1:ind1+2)

  end if

  ! compute turnover
  call z_rot3_turnover(G1,G2,G3)

  ! move misfit to appropriate side based on P(2)
  ! hess
  if (.NOT.P(ii+1)) then
 
    ! update Q
    Q(ind1:ind1+2) = G1
    Q(ind2:ind2+2) = G2

    ! update left schurvectors with G3
    if (VECL) then 
    
      A(1,1) = cmplx(G3(1),G3(2),kind=8)
      A(2,1) = cmplx(G3(3),0d0,kind=8)
      A(1,2) = cmplx(-G3(3),0d0,kind=8)
      A(2,2) = cmplx(G3(1),-G3(2),kind=8)

      W = matmul(W,A)
   
    end if

    ! set MISFIT as inverse of G3
    MISFIT(1) = G3(1)
    MISFIT(2) = -G3(2)
    MISFIT(3) = -G3(3)

    ! pass MISFIT through R2
    !call z_upr1utri_rot3swap(.TRUE.,D2,C2,B2,MISFIT)
    call z_uprkutri_rot3swap(.TRUE.,N,K,1,K,D2,C2,B2,MISFIT,ii+1)

    ! invert MISFIT
    MISFIT(2) = -MISFIT(2)
    MISFIT(3) = -MISFIT(3)
    
    ! update right schurvectors with MISFIT
    if (VECR) then
    
      A(1,1) = cmplx(MISFIT(1),MISFIT(2),kind=8)
      A(2,1) = cmplx(MISFIT(3),0d0,kind=8)
      A(1,2) = cmplx(-MISFIT(3),0d0,kind=8)
      A(2,2) = cmplx(MISFIT(1),-MISFIT(2),kind=8)

      V = matmul(V,A)
   
    end if

    ! pass MISFIT through R1
    !call z_upr1utri_rot3swap(.FALSE.,D1,C1,B1,MISFIT)
    call z_uprkutri_rot3swap(.FALSE.,N,K,1,K,D1,C1,B1,MISFIT,ii+1)

  ! invhess
  else

    ! update Q
    Q(ind1:ind1+2) = G1
    Q(ind2:ind2+2) = G3

    ! pass G2 through R1
    !call z_upr1utri_rot3swap(.TRUE.,D1,C1,B1,G2)
    call z_uprkutri_rot3swap(.TRUE.,N,K,1,K,D1,C1,B1,G2,ii+1)

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
    !call z_upr1utri_rot3swap(.FALSE.,D2,C2,B2,G2)
    call z_uprkutri_rot3swap(.FALSE.,N,K,1,K,D2,C2,B2,G2,ii+1)

    ! update left schurvectors using G2
    if (VECL) then
     
      A(1,1) = cmplx(G2(1),G2(2),kind=8)
      A(2,1) = cmplx(G2(3),0d0,kind=8)
      A(1,2) = cmplx(-G2(3),0d0,kind=8)
      A(2,2) = cmplx(G2(1),-G2(2),kind=8)

      W = matmul(W,A)

    end if

    ! copy inverse of G2 to MISFIT
    MISFIT(1) = G2(1)
    MISFIT(2) = -G2(2)
    MISFIT(3) = -G2(3)

  end if

  ! update position flag
  P(ii) = P(ii+1)

end subroutine z_uprkfpen_chasedown
