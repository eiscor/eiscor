#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_2x2deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine deflates a 2x2 block in a upr1 pencil
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  QZ              LOGICAL
!                    .TRUE.: second triangular factor is assumed nonzero
!                    .FALSE.: second triangular factor is assumed to be identity
!
!  VEC             LOGICAL
!                    .TRUE.: schurvectors, assume V and W already initialized
!                    .FALSE.: no schurvectors
!
!  Q               REAL(8) array of dimension (3)
!                    array of generators for first sequence of rotations
!
!  D1,D2           REAL(8) arrays of dimension (4)
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1,B1,C2,B2     REAL(8) arrays of dimension (6)
!                    array of generators for upper-triangular parts
!                    of the pencil
!
!  M               INTEGER
!                    leading dimension of V and W
!
! OUTPUT VARIABLES:
!
!  V               COMPLEX(8) array of dimension (M,2)
!                    right schur vectors
!                    if QZ = .FALSE. unused
!                    if VEC = .TRUE. stores right schurvectors in V 
!                    if VEC = .FALSE. unused
!
!  W               COMPLEX(8) array of dimension (M,2)
!                    left schur vectors
!                    if QZ = .FALSE. unused
!                    if VEC = .TRUE. stores right schurvectors in V 
!                    if VEC = .FALSE. unused
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_2x2deflation(QZ,VEC,Q,D1,C1,B1,D2,C2,B2,M,V,W)

  implicit none
  
  ! input variables
  logical, intent(in) :: QZ, VEC
  real(8), intent(inout) :: Q(3), D1(4), C1(6), B1(6)
  real(8), intent(inout) :: D2(4), C2(6), B2(6)
  integer, intent(in) :: M 
  complex(8), intent(inout) :: V(M,2), W(M,2)
  
  ! compute variables
  real(8) :: nrm, G1(3), G2(3), G3(3), Qt(6)
  complex(8) :: A(2,2), B(2,2), Vt(2,2), Wt(2,2)
  
  ! compute 2x2 blocks
  Qt = 0d0; Qt(1:3) = Q; Qt(4) = 1d0
  call z_upr1fact_2x2diagblocks(.TRUE.,.TRUE.,QZ,.FALSE.,Qt,D1,C1,B1,D2,C2,B2,A,B)

  ! compute standard Schur decomposition
  if (.NOT.QZ) then
  
    ! standard schur decomposition
    call z_2x2array_eig(QZ,A,B,Vt,Wt)

    ! replace Vt with rotation G1
    call z_rot3_vec4gen(dble(Vt(1,1)),aimag(Vt(1,1)),dble(Vt(2,1)),aimag(Vt(2,1)),G1(1),G1(2),G1(3),nrm)
    
    ! pass G1 through triangular factor
    G3 = G1
    call z_upr1fact_rot3throughtri(.FALSE.,D1,C1,B1,G3)
    
    ! similarity transform of Q
    A(1,1) = cmplx(Q(1),Q(2),kind=8)
    A(2,1) = cmplx(Q(3),0d0,kind=8)
    A(1,2) = -A(2,1)
    A(2,2) = conjg(A(1,1))
    
    B(1,1) = cmplx(G3(1),G3(2),kind=8)
    B(2,1) = cmplx(G3(3),0d0,kind=8)
    B(1,2) = -B(2,1)
    B(2,2) = conjg(B(1,1))
    
    A = matmul(A,B)
    
    ! update V
    if (VEC) then
    
      B(1,1) = cmplx(G1(1),G1(2),kind=8)
      B(2,1) = cmplx(G1(3),0d0,kind=8)
      B(1,2) = -B(2,1)
      B(2,2) = conjg(B(1,1))
      
      V = matmul(V,B)
    
    end if

    B(1,1) = cmplx(G1(1),-G1(2),kind=8)
    B(2,1) = cmplx(-G1(3),0d0,kind=8)
    B(1,2) = -B(2,1)
    B(2,2) = conjg(B(1,1))    

    A = matmul(B,A)
    
    Q(1) = 1d0
    Q(2) = 0d0
    Q(3) = 0d0
        
    ! deflate into D1
    A(1,1) = A(1,1)*cmplx(D1(1),D1(2),kind=8)
    call d_rot2_vec2gen(dble(A(1,1)),aimag(A(1,1)),D1(1),D1(2),nrm)
    
    A(2,2) = A(2,2)*cmplx(D1(3),D1(4),kind=8)
    call d_rot2_vec2gen(dble(A(2,2)),aimag(A(2,2)),D1(3),D1(4),nrm)
    
  ! compute generalized Schur decomposition
  else
  
    ! generalized schur decomposition
    call z_2x2array_eig(QZ,A,B,Vt,Wt)
      
    ! replace Vt with rotation G1
    call z_rot3_vec4gen(dble(Wt(1,1)),aimag(Wt(1,1)),dble(Wt(2,1)),aimag(Wt(2,1)),G1(1),G1(2),G1(3),nrm)
    
    ! replace Wt with rotation G2
    call z_rot3_vec4gen(dble(Vt(1,1)),aimag(Vt(1,1)),dble(Vt(2,1)),aimag(Vt(2,1)),G2(1),G2(2),G2(3),nrm)
    
    ! pass G1 through right triangular factor
    G3 = G1
    call z_upr1fact_rot3throughtri(.FALSE.,D2,C2,B2,G3)
    
    ! equivalence transform
    A(1,1) = cmplx(G2(1),G2(2),kind=8)
    A(2,1) = cmplx(G2(3),0d0,kind=8)
    A(1,2) = -A(2,1)
    A(2,2) = conjg(A(1,1))
    
    B(1,1) = cmplx(G3(1),G3(2),kind=8)
    B(2,1) = cmplx(G3(3),0d0,kind=8)
    B(1,2) = -B(2,1)
    B(2,2) = conjg(B(1,1))
    
    A = matmul(transpose(conjg(A)),B)
    
    ! deflate into D2
    A(1,1) = A(1,1)*cmplx(D2(1),D2(2),kind=8)
    call d_rot2_vec2gen(dble(A(1,1)),aimag(A(1,1)),D2(1),D2(2),nrm)
    
    A(2,2) = A(2,2)*cmplx(D2(3),D2(4),kind=8)
    call d_rot2_vec2gen(dble(A(2,2)),aimag(A(2,2)),D2(3),D2(4),nrm)
    
    ! pass G1 through left triangular factor
    G3 = G1
    call z_upr1fact_rot3throughtri(.FALSE.,D1,C1,B1,G3)
    
    ! equivalence transform of Q
    A(1,1) = cmplx(Q(1),Q(2),kind=8)
    A(2,1) = cmplx(Q(3),0d0,kind=8)
    A(1,2) = -A(2,1)
    A(2,2) = conjg(A(1,1))
    
    B(1,1) = cmplx(G3(1),G3(2),kind=8)
    B(2,1) = cmplx(G3(3),0d0,kind=8)
    B(1,2) = -B(2,1)
    B(2,2) = conjg(B(1,1))
    
    A = matmul(A,B)
    
    B(1,1) = cmplx(G2(1),-G2(2),kind=8)
    B(2,1) = cmplx(-G2(3),0d0,kind=8)
    B(1,2) = -B(2,1)
    B(2,2) = conjg(B(1,1))
    
    A = matmul(B,A)
    
    Q(1) = 1d0
    Q(2) = 0d0
    Q(3) = 0d0
    
    ! deflate into D1
    A(1,1) = A(1,1)*cmplx(D1(1),D1(2),kind=8)
    call d_rot2_vec2gen(dble(A(1,1)),aimag(A(1,1)),D1(1),D1(2),nrm)
    
    A(2,2) = A(2,2)*cmplx(D1(3),D1(4),kind=8)
    call d_rot2_vec2gen(dble(A(2,2)),aimag(A(2,2)),D1(3),D1(4),nrm)
    
    ! update vecs
    if (VEC) then
    
      ! update V
      B(1,1) = cmplx(G1(1),G1(2),kind=8)
      B(2,1) = cmplx(G1(3),0d0,kind=8)
      B(1,2) = -B(2,1)
      B(2,2) = conjg(B(1,1))
      
      W = matmul(W,B) 
      
      ! update W
      B(1,1) = cmplx(G2(1),G2(2),kind=8)
      B(2,1) = cmplx(G2(3),0d0,kind=8)
      B(1,2) = -B(2,1)
      B(2,2) = conjg(B(1,1))
      
      V = matmul(V,B)
    
    end if
  
  end if

end subroutine z_upr1fact_2x2deflation
