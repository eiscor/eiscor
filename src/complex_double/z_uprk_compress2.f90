#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprk_compress2.f90
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This functions computes the factorized representation of a matrix polynomial 
! B_1 \lambda^d + (B_2 + A_1) \lambda^d-1 + ... + (A_(d-1) + B_d \lambda + A_d
! of degree d with B_1 lower triangular and A_d upper triangular 
!
! If not QZ then B_1 = I and B_i = 0 for i>1.
!
! After factorizing the problem in d upper Hessenberg matrices for A and
! d upper triangular matrices for B, the problem is transformed to a product 
! eigenvalue problem as a product of one Hessenberg matrix and 2*d-1 upper 
! trinagular matrices. This is done by similarity transformations.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  QZ              LOGICAL
!                    .TRUE.: second triangular factor is assumed nonzero
!                    .FALSE.: second triangular factor is assumed to be identity
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvector
!                    .FALSE.: no schurvectors
!
!  ID              LOGICAL
!                    .TRUE.: initialize V and W to I
!                    .FALSE.: assume V and W already initialized
!
!  N               INTEGER
!                    dimension of matrix (K*d)
!
!  K               INTEGER
!                    rank, i.e., number of upper triangulars
!
!  Ain,Bin         COMPLEX(8) array of dimension (N,K)
!                    coefficients matrix polynomial
! 
! OUTPUT VARIABLES:
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D1,D2           REAL(8) arrays of dimension (2*K*N)
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1,B1,C2,B2     REAL(8) arrays of dimension (3*K*N)
!                    array of generators for upper-triangular parts
!                    of the pencil
!
!  V               COMPLEX(8) array of dimension (N,N)
!                    right schur vectors
!
!  W               COMPLEX(8) array of dimension (N,N)
!                    left schur vectors
!
!  INFO            INTEGER
!                    TBA   
!                    INFO = 1 implies no convergence 
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies N, Q, D1, C1, B1, D2, C2 or B2 is invalid
!                    INFO = -9 implies V is invalid
!                    INFO = -10 implies W is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprk_compress2(QZ,VEC,ID,N,K,Ain,Bin,P,Q,D1,&
     &C1,B1,D2,C2,B2,V,W,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: QZ, VEC, ID
  integer, intent(in) :: N,K
  complex(8), intent(in) :: Ain(N,K), Bin(N,K)
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*K*(N-1)), D1(2*K*N), C1(3*N*K), B1(3*N*K)
  real(8), intent(inout) :: D2(2*K*N), C2(3*N*K), B2(3*N*K)
  complex(8), intent(inout) :: V(N,N), W(N,N)
  integer, intent(inout) :: INFO

  ! compute variables
  complex(8) :: A(N,K), B(N,K), H(2,2)
!  complex(8) :: MV(N,N), MW(N,N)
  real(8) :: bulge(3), D0(2*K*N),dr,di,dr1,di1,dr2,di2,nrm
  complex(8) :: hc(K)!, Gf(2,2)
  !complex(8) :: REIGS2(N), WORK(5*N)
  !real(8) :: RWORK(2*N)
  integer :: ii, jj, ll, row, col, ind, ind1, ind0 !, lwork

  ! initialize info
  INFO = 0

  A = Ain
  B = Bin
  
  ! fill P
  P = .FALSE.

  ! initialize Q
  Q = 0d0
  do ii = 1,K*(N-1)
    Q(3*ii) = 1d0
  end do

  ! initialize D0
  D0 = 0d0
  do ii = 1,k*N
    D0(2*ii-1) = 1d0
  end do


  ! initalize V
  if (VEC .AND. ID) then
     V = cmplx(0d0,0d0,kind=8)
     do ii=1,N
        V(ii,ii) = cmplx(1d0,0d0,kind=8)
     end do
  end if

  ! initalize W
  if (QZ .AND. VEC .AND. ID) then
     W = cmplx(0d0,0d0,kind=8)
     do ii=1,N
        W(ii,ii) = cmplx(1d0,0d0,kind=8)
     end do
  end if

  if (DEBUG) then
     print*, "updated A"
     do ii=1,N
        print*, A(ii,:)
     end do
  end if

  ! initalize unitary-plus-rank-1 upper triangulars A
  do ll = 1,K
     hc(1:K) = A(1:K,ll)
     A(1:N-K,ll) = A(K+1:N,ll) 
     A(N-K+1:N,ll) = hc(1:K)*(-1)**(N+1)
  end do

  do ll = 1,K
     call z_idpspike_compress(.TRUE.,N,N+1-ll,A(1:N,K+1-ll),&
          &D1(((ll-1)*2*N+1):(ll*2*N)),&
          &C1(((ll-1)*3*N+1):(ll*3*N)),&
          &B1(((ll-1)*3*N+1):(ll*3*N)),INFO)
  end do

  ! initalize unitary-plus-rank-1 upper triangulars B
  if (QZ) then
     do ll = 1,K
        call z_idpspike_compress(.TRUE.,N,N+1-ll,B(1:N,K+1-ll),&
             &D2(((ll-1)*2*N+1):(ll*2*N)),&
             &C2(((ll-1)*3*N+1):(ll*3*N)),&
             &B2(((ll-1)*3*N+1):(ll*3*N)),INFO)
     end do
  end if

!!$  if (DEBUG) then
!!$     ! Check V
!!$     do ll = 1,K
!!$        call z_upr1fact_extracttri(.FALSE.,N,&
!!$             &D1(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),&
!!$             &C1(((ll-1)*3*N+1):(ll*3*N)),&
!!$             &B1(((ll-1)*3*N+1):(ll*3*N)),MV)
!!$        do jj = N-1,1,-1
!!$           ! Apply Q(jj)
!!$           Gf(1,1) = cmplx(Q(3*jj-2),Q(3*jj-1),kind=8)
!!$           Gf(2,1) = cmplx(Q(3*jj),0d0,kind=8)
!!$           Gf(1,2) = -Gf(2,1)
!!$           Gf(2,2) = conjg(Gf(1,1))
!!$           
!!$           MV((jj):(jj+1),:) = matmul(Gf,MV((jj):(jj+1),:))
!!$        end do
!!$        if (ll.EQ.1) then
!!$           V = MV
!!$        else
!!$           V = matmul(V,MV)
!!$        end if
!!$     end do
!!$     print*, "V before reduction"
!!$     do ii=1,N
!!$        write (*,"(ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3)") V(ii,:)
!!$        print*, ""
!!$     end do
!!$     
!!$     ! Check W
!!$     do ll = 1,K
!!$        call z_upr1fact_extracttri(.FALSE.,N,&
!!$             &D2(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C2(((ll-1)*3*N+1):(ll*3*N)),&
!!$             &B2(((ll-1)*3*N+1):(ll*3*N)),V)
!!$        if (ll.EQ.1) then
!!$           W = V
!!$        else
!!$           W = matmul(W,V)
!!$        end if
!!$     end do
!!$     print*, "W"
!!$     do ii=1,N
!!$        write (*,"(ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3)") V(ii,:)
!!$        print*, ""
!!$        !write (*,"(ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3)") W(ii,:)
!!$     end do
!!$
!!$     ! check factorization
!!$     ! Form upper triangulars and multiply them together
!!$     MV = cmplx(0d0,0d0,kind=8)
!!$     MW = cmplx(0d0,0d0,kind=8)
!!$     do ll = 1,N
!!$        MV(ll,ll) = cmplx(1d0,0d0,kind=8)
!!$     end do
!!$     do ll = 1,K
!!$        call z_upr1fact_extracttri(.false.,N,&
!!$             &D1(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C1(((ll-1)*3*N+1):(ll*3*N)),&
!!$             &B1(((ll-1)*3*N+1):(ll*3*N)),MW)
!!$        ! Apply Q
!!$        do jj = N-1,1,-1
!!$           ! Apply Q(jj)
!!$           Gf(1,1) = cmplx(Q((ll-1)*(N-1)+3*jj-2),Q((ll-1)*(N-1)+3*jj-1),kind=8)
!!$           Gf(2,1) = cmplx(Q((ll-1)*(N-1)+3*jj),0d0,kind=8)
!!$           Gf(1,2) = -Gf(2,1)
!!$           Gf(2,2) = conjg(Gf(1,1))
!!$           
!!$           MW((jj):(jj+1),:) = matmul(Gf,MW((jj):(jj+1),:))
!!$        end do
!!$        
!!$        ! print*, "MW", ll
!!$        ! do jj = 1,N
!!$        !    print*, jj, MW(jj,:)
!!$        ! end do
!!$        
!!$        MV = matmul(MV,MW)
!!$     end do
!!$     
!!$     print*, "MV"
!!$     do jj = 1,N
!!$        print*, jj, MV(jj,:)
!!$     end do
!!$  end if
     
     
  ! reduce to Hessenberg * Triangular * ... * Triangular form
  ! remove rotations in row jj
  do jj = 1,(N-1)
     ! remove rotation from Q_ii starting with K
     do ii = K,2,-1
        col = ii-1 ! rotation right of col th upper triangular
        bulge = Q((3*col*(N-1)+3*(jj-1)+1):(3*col*(N-1)+3*jj))
        row = jj
        do while ((row.LE.N-2).AND.(bulge(3).NE.0d0)) 
           !print*, "row", row, "col", col
           !do row = jj,N-2 ! first row of the rotation
           ! pass through diagonal D0
           ind0 = (col-1)*2*N+(row-1)*2
           call z_rot3_swapdiag(D0((ind0+1):(ind0+4)),bulge)
           if (bulge(3).NE.0d0) then
              ! turnover through Q
              call z_rot3_turnover(&
                   &Q((3*(col-1)*(N-1)+3*(row-1)+1):(3*(col-1)*(N-1)+3*row)),&
                   &Q((3*(col-1)*(N-1)+3*(row-1)+4):(3*(col-1)*(N-1)+3*row+3)),bulge)
              ! update col
              col = col - 1
              ! update col if left most sequence was passed
              if (col.EQ.0) then
                 if (QZ) then
                    ! update (left) eigenvectors 
                    ! update W
                    if (VEC) then
                       
                       H(1,1) = cmplx(bulge(1),bulge(2),kind=8)
                       H(2,1) = cmplx(bulge(3),0d0,kind=8)
                       H(1,2) = -H(2,1)
                       H(2,2) = conjg(H(1,1))
                    
                       W(:,row+1:row+2) = matmul(W(:,row+1:row+2),H)
                       
                    end if
                    ! invert bulge
                    bulge(1) = bulge(1)
                    bulge(2) = -bulge(2) 
                    bulge(3) = -bulge(3)
                    ! through B
                    call z_uprkutri_rot3swap(.TRUE.,N,K,1,K,D2,C2,B2,bulge,row+1)
                    ! invert bulge
                    bulge(1) = bulge(1)
                    bulge(2) = -bulge(2) 
                    bulge(3) = -bulge(3)
                 end if
                 ! update (right) eigenvectors 
                 ! update V
                 if (VEC) then
                    
                    H(1,1) = cmplx(bulge(1),bulge(2),kind=8)
                    H(2,1) = cmplx(bulge(3),0d0,kind=8)
                    H(1,2) = -H(2,1)
                    H(2,2) = conjg(H(1,1))
                    
                    V(:,row+1:row+2) = matmul(V(:,row+1:row+2),H)
                    
                 end if
                 ! pass through triangulars without rotations in front
                 ! through A
                 call z_uprkutri_rot3swap(.FALSE.,N,K,1,K,D1,C1,B1,bulge,row+1)
                 col = K                 
              end if
           end if
           row = row + 1
        end do

        if (bulge(3).NE.0d0) then
           ! fuse rotation into current rotation
           ind0 = (col-1)*2*N+(row-1)*2
           call z_rot3_swapdiag(D0((ind0+1):(ind0+4)),bulge)
           call z_rot3_fusion(.TRUE., Q((col*3*(N-1)-2):(col*3*(N-1))), bulge)
           ! scale rows of R1
           ind = col*3*N
           ind1 = col*2*N
           call z_upr1utri_unimodscale(.TRUE.,D0(ind1-3:ind1-2), &
                C1(ind-5:ind-3), & 
                B1(ind-5:ind-3), &
                cmplx(bulge(1),bulge(2),kind=8))
           call z_upr1utri_unimodscale(.TRUE.,D0(ind1-1:ind1), &
                C1(ind-2:ind), & 
                B1(ind-2:ind), &
                cmplx(bulge(1),-bulge(2),kind=8))
        end if

     end do
  end do

  ! fuse D0 into D1
  do ii=1,N
     dr = D0(2*ii-1)
     di = D0(2*ii)
     do col = 2,K
        dr1 = dr
        di1 = di
        dr2 = D0((col-1)*2*N+2*ii-1)
        di2 = D0((col-1)*2*N+2*ii)
     
        dr = dr1*dr2 - di1*di2
        di = dr1*di2 + di1*dr2
     end do
     dr1 = dr
     di1 = di
     dr2 = D1(2*ii-1)
     di2 = D1(2*ii)
     
     dr = dr1*dr2 - di1*di2
     di = dr1*di2 + di1*dr2
     call d_rot2_vec2gen(dr,di,D1(2*ii-1),D1(2*ii),nrm) 
  end do

  
end subroutine z_uprk_compress2
