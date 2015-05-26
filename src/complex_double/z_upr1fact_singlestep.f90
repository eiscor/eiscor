#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_singlestep 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift 
! algorithm on a upr1 pencil. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  QZ              LOGICAL
!                    .TRUE. second triangular factor is assumed nonzero
!                    .FALSE. second triangular factor is assumed to be identity
!
!  VEC             LOGICAL
!                    .TRUE. update schurvectors
!                    .FALSE. no schurvectors
!
!  FUN             LOGICAL FUNCTION FUN(N,P)
!                    takes integer N and logical array P of 
!                    dimension N-2 and outputs a logical 
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
!  D1,D2           REAL(8) arrays of dimension (2*(N+1))
!                    array of generators for complex diagonal matrices
!                    If QZ = .FALSE., D2 is unused.
!
!  C1,B1,C2,B2     REAL(8) arrays of dimension (3*N)
!                    array of generators for upper-triangular parts of the pencil
!                    If QZ = .FALSE., C2 and B2 are unused.
!
!  M               INTEGER
!                    leading dimesnion of V and W
!
!  V               COMPLEX(8) array of dimension (M,N)
!                    right schur vectors
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. update V to store right schurvectors 
!
!  W               COMPLEX(8) array of dimension (M,N)
!                    left schur vectors
!                    if QZ = .FALSE. unused
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. update W to store left schurvectors
!
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_singlestep(QZ,VEC,FUN,N,P,Q,D1,C1,B1,D2,C2,B2,M,V,W,ITCNT)

  implicit none
  
  ! input variables
  logical, intent(in) :: QZ, VEC
  integer, intent(in) :: M, N
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*(N+1)), C1(3*N), B1(3*N)
  real(8), intent(inout) :: D2(2*(N+1)), C2(3*N), B2(3*N)
  complex(8), intent(inout) :: V(M,N), W(M,N)
  integer, intent(in) :: ITCNT
  interface
    function FUN(m,flags)
      logical :: FUN
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function FUN
  end interface
  
  ! compute variables
  integer :: ii, ir1, ir2, id1, id2
  logical :: final_flag
  real(8) :: G1(3), G2(3), G3(3)
  complex(8) :: shift, rho, A(2,2), B(2,2), Vt(2,2), Wt(2,2)
  
  ! compute final_flag
  final_flag = FUN(N,P)
  
  ! compute shift
  ! random shift
  if(mod(ITCNT+1,11) == 0)then
    call random_number(G1(1))
    call random_number(G1(2))
    shift = cmplx(G1(1),G1(2),kind=8)
          
  ! wilkinson shift
  else
  
    ! get 2x2 blocks
    ir2 = 3*N; ir1 = ir2-5
    id2 = 2*N; id1 = id2-3
    call z_upr1fact_2x2diagblocks(.FALSE.,.TRUE.,QZ,P(N-2) &
    ,Q((ir1-3):(ir2-3)),D1(id1:id2),C1(ir1:ir2),B1(ir1:ir2) &
    ,D2(id1:id2),C2(ir1:ir2),B2(ir1:ir2),A,B)
  
    ! if not QZ
    if (.NOT.QZ) then
      B = cmplx(0d0,0d0,kind=8)
      B(1,1) = cmplx(1d0,0d0,kind=8); B(2,2) = B(1,1)
    end if
    
    ! store ratio of bottom right entries
    rho = A(2,2)/B(2,2) 
      
    ! compute eigenvalues and eigenvectors
    call z_2x2array_eig(QZ,A,B,Vt,Wt)
          
    ! wilkinson shift
    if(abs(A(2,2)/B(2,2)-rho) < abs(A(1,1)/B(1,1)-rho))then
      shift = A(2,2)/B(2,2)
    else
      shift = A(1,1)/B(1,1)
    end if
    
    ! avoid INFs and NANs
    if ((shift.NE.shift).OR.(abs(shift) > EISCOR_DBL_INF)) then
      shift = cmplx(1d0,0d0,kind=8) ! not sure if this is a good idea?
    end if

  end if

  ! build bulge
  call z_upr1fact_buildbulge(QZ,P(1),Q(1:6),D1(1:4),C1(1:6),B1(1:6) &
  ,D2(1:4),C2(1:6),B2(1:6),shift,G2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! iteration for QZ
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (QZ) then

    ! update W
    if (VEC) then
      
      A(1,1) = cmplx(G2(1),G2(2),kind=8)
      A(2,1) = cmplx(G2(3),0d0,kind=8)
      A(1,2) = -A(2,1)
      A(2,2) = conjg(A(1,1))
      
      W(:,1:2) = matmul(W(:,1:2),A)
      
    end if

    ! initialize turnover
    ! hess
    if (.NOT.P(1)) then
    
      ! set G2 as G2^-1 
      G2(2) = -G2(2)
      G2(3) = -G2(3)
    
      ! merge from left
      call z_upr1fact_mergebulge(.TRUE.,N,P,Q,D1(1:(2*N)),G2)
      
      ! set G1 for turnover
      G1(1) = Q(1)
      G1(2) = Q(2)
      G1(3) = Q(3)

      ! pass G2 through R2
      call z_upr1fact_rot3throughtri(.TRUE.,D2(1:4),C2(1:6),B2(1:6),G2)
    
      ! set G2 as G2^-1 
      G2(2) = -G2(2)
      G2(3) = -G2(3)
    
      ! update V
      if (VEC) then
      
        A(1,1) = cmplx(G2(1),G2(2),kind=8)
        A(2,1) = cmplx(G2(3),0d0,kind=8)
        A(1,2) = -A(2,1)
        A(2,2) = conjg(A(1,1))
      
        V(:,1:2) = matmul(V(:,1:2),A)
      
      end if
    
      ! pass G2 through R1
      call z_upr1fact_rot3throughtri(.FALSE.,D1(1:4),C1(1:6),B1(1:6),G2)
      
      ! set G3 for turnover
      G3 = G2

      ! set G2 for turnover
      G2 = Q(4:6)
  
    ! inverse hess
    else

      ! set G2 as G2^-1 
      G2(2) = -G2(2)
      G2(3) = -G2(3)
    
      ! set G1 for turnover
      G1 = G2

      ! pass G2 through R2
      call z_upr1fact_rot3throughtri(.TRUE.,D2(1:4),C2(1:6),B2(1:6),G2)
    
      ! set G2 as G2^-1 
      G2(2) = -G2(2)
      G2(3) = -G2(3)
    
      ! update V
      if (VEC) then
      
        A(1,1) = cmplx(G2(1),G2(2),kind=8)
        A(2,1) = cmplx(G2(3),0d0,kind=8)
        A(1,2) = -A(2,1)
        A(2,2) = conjg(A(1,1))
      
        V(:,1:2) = matmul(V(:,1:2),A)
      
      end if
    
      ! pass G2 through R1
      call z_upr1fact_rot3throughtri(.FALSE.,D1(1:4),C1(1:6),B1(1:6),G2)
      
      ! merge from right
      call z_upr1fact_mergebulge(.TRUE.,N,P,Q,D1(1:(2*N)),G2)
      
      ! set G3 for turnover
      G3 = Q(1:3)

      ! set G2 for turnover
      G2 = Q(4:6)
  
    end if

    ! chase bulge
    do ii=1,(N-3)
    
      ! execute turnover of G1*G2*G3
      call z_rot3_turnover(G1,G2,G3)
      
      ! prepare for next turnover based on P(ii+1)
      ! hess
      if (.NOT.P(ii+1)) then
      
        ! set P(ii)
        P(ii) = P(ii+1)
        
        ! set Q(ii)
        Q(3*ii-2) = G1(1)
        Q(3*ii-1) = G1(2)
        Q(3*ii) = G1(3)
      
        ! set Q(ii+1)
        Q(3*ii+1) = G2(1)
        Q(3*ii+2) = G2(2)
        Q(3*ii+3) = G2(3)        
 
        ! set G1 for turnover
        G1 = G2     
        
        ! set G2 for turnover
        G2(1) = Q(3*ii+4)
        G2(2) = Q(3*ii+5)
        G2(3) = Q(3*ii+6)
        
        ! update W
        if (VEC) then
          
          A(1,1) = cmplx(G3(1),G3(2),kind=8)
          A(2,1) = cmplx(G3(3),0d0,kind=8)
          A(1,2) = -A(2,1)
          A(2,2) = conjg(A(1,1))
          
          W(:,(ii+1):(ii+2)) = matmul(W(:,(ii+1):(ii+2)),A)
          
        end if

        ! G3 = G3^-1
        G3(2) = -G3(2)
        G3(3) = -G3(3)
        
        ! pass G3 through R2
        call z_upr1fact_rot3throughtri(.TRUE.,D2((2*ii+1):(2*ii+4)) &
        ,C2((3*ii+1):(3*ii+6)),B2((3*ii+1):(3*ii+6)),G3)
        
        ! G3 = G3^-1
        G3(2) = -G3(2)
        G3(3) = -G3(3)
        
        ! update V
        if (VEC) then
          
          A(1,1) = cmplx(G3(1),G3(2),kind=8)
          A(2,1) = cmplx(G3(3),0d0,kind=8)
          A(1,2) = -A(2,1)
          A(2,2) = conjg(A(1,1))
          
          V(:,(ii+1):(ii+2)) = matmul(V(:,(ii+1):(ii+2)),A)
          
        end if

        ! pass G3 through R1
        call z_upr1fact_rot3throughtri(.FALSE.,D1((2*ii+1):(2*ii+4)) &
        ,C1((3*ii+1):(3*ii+6)),B1((3*ii+1):(3*ii+6)),G3)
        
      ! inverse hess
      else
      
        ! set P(ii)
        P(ii) = P(ii+1)
        
        ! set Q(ii)
        Q(3*ii-2) = G1(1)
        Q(3*ii-1) = G1(2)
        Q(3*ii) = G1(3)
      
        ! set Q(ii+1)
        Q(3*ii+1) = G3(1)
        Q(3*ii+2) = G3(2)
        Q(3*ii+3) = G3(3)  
        
        ! set G3 for turnover
        G3 = G3    
        
        ! pass G2 through R1
        call z_upr1fact_rot3throughtri(.TRUE.,D1((2*ii+1):(2*ii+4)) &
        ,C1((3*ii+1):(3*ii+6)),B1((3*ii+1):(3*ii+6)),G2)
        
        ! G2 = G2^-1
        G2(2) = -G2(2)
        G2(3) = -G2(3)

        ! update V
        if (VEC) then
          
          A(1,1) = cmplx(G2(1),G2(2),kind=8)
          A(2,1) = cmplx(G2(3),0d0,kind=8)
          A(1,2) = -A(2,1)
          A(2,2) = conjg(A(1,1))
          
          V(:,(ii+1):(ii+2)) = matmul(V(:,(ii+1):(ii+2)),A)
          
        end if

        ! pass G2 through R2
        call z_upr1fact_rot3throughtri(.FALSE.,D2((2*ii+1):(2*ii+4)) &
        ,C2((3*ii+1):(3*ii+6)),B2((3*ii+1):(3*ii+6)),G2)
        
        ! update W
        if (VEC) then
          
          A(1,1) = cmplx(G2(1),G2(2),kind=8)
          A(2,1) = cmplx(G2(3),0d0,kind=8)
          A(1,2) = -A(2,1)
          A(2,2) = conjg(A(1,1))
          
          V(:,(ii+1):(ii+2)) = matmul(V(:,(ii+1):(ii+2)),A)
          
        end if

        ! G2 = G2^-1
        G2(2) = -G2(2)
        G2(3) = -G2(3)

        ! set G1 for turnover
        G1 = G2

        ! set G2 for turnover
        G2(1) = Q(3*ii+4)
        G2(2) = Q(3*ii+5)
        G2(3) = Q(3*ii+6)
        
      end if

    end do  
  
    ! final turnover
    call z_rot3_turnover(G1,G2,G3)
      
    ! set P(N-1)
    P(N-2) = final_flag
    
    ! finish transformation based on P(N-1)
    ! hess
    if (.NOT.P(N-2)) then
    
      ! set Q(N-2)
      Q(3*(N-2)-2) = G1(1)
      Q(3*(N-2)-1) = G1(2)
      Q(3*(N-2)) = G1(3)   
      
      ! set Q(N-1)
      Q(3*(N-1)-2) = G2(1)
      Q(3*(N-1)-1) = G2(2)
      Q(3*(N-1)) = G2(3)  
      
      ! update W
      if (VEC) then
          
        A(1,1) = cmplx(G3(1),G3(2),kind=8)
        A(2,1) = cmplx(G3(3),0d0,kind=8)
        A(1,2) = -A(2,1)
        A(2,2) = conjg(A(1,1))
          
        W(:,(N-1):N) = matmul(W(:,(N-1):N),A)
         
      end if    
   
      ! G3 = G3^-1
      G3(2) = -G3(2)
      G3(3) = -G3(3)
 
      ! pass G3 through R2
      call z_upr1fact_rot3throughtri(.TRUE.,D2((2*N-3):(2*N)) &
      ,C2((3*N-5):(3*N)),B2((3*N-5):(3*N)),G3)
        
      ! G3 = G3^-1
      G3(2) = -G3(2)
      G3(3) = -G3(3)
 
      ! update V
      if (VEC) then
          
        A(1,1) = cmplx(G3(1),G3(2),kind=8)
        A(2,1) = cmplx(G3(3),0d0,kind=8)
        A(1,2) = -A(2,1)
        A(2,2) = conjg(A(1,1))
          
        V(:,(N-1):N) = matmul(V(:,(N-1):N),A)
         
      end if    
    
      ! pass G3 through R1
      call z_upr1fact_rot3throughtri(.FALSE.,D1((2*N-3):(2*N)) &
      ,C1((3*N-5):(3*N)),B1((3*N-5):(3*N)),G3)
        
      ! merge bulge 
      call z_upr1fact_mergebulge(.FALSE.,N,P,Q,D1(1:(2*N)),G3)
      
    ! inverse hess
    else
    
      ! set Q(N-2)
      Q(3*(N-2)-2) = G1(1)
      Q(3*(N-2)-1) = G1(2)
      Q(3*(N-2)) = G1(3)   
      
      ! set Q(N-1)
      Q(3*(N-1)-2) = G3(1)
      Q(3*(N-1)-1) = G3(2)
      Q(3*(N-1)) = G3(3)  
      
      ! pass G2 through R1
      call z_upr1fact_rot3throughtri(.TRUE.,D1((2*N-3):(2*N)) &
      ,C1((3*N-5):(3*N)),B1((3*N-5):(3*N)),G2)
        
      ! G2 = G2^-1
      G2(2) = -G2(2)
      G2(3) = -G2(3)
 
      ! update V
      if (VEC) then
          
        A(1,1) = cmplx(G2(1),G2(2),kind=8)
        A(2,1) = cmplx(G2(3),0d0,kind=8)
        A(1,2) = -A(2,1)
        A(2,2) = conjg(A(1,1))
          
        V(:,N:(N+1)) = matmul(V(:,N:(N+1)),A)
         
      end if  
      
      ! pass G2 through R2
      call z_upr1fact_rot3throughtri(.FALSE.,D2((2*N-3):(2*N)) &
      ,C2((3*N-5):(3*N)),B2((3*N-5):(3*N)),G2)
        
      ! update W
      if (VEC) then
          
        A(1,1) = cmplx(G2(1),G2(2),kind=8)
        A(2,1) = cmplx(G2(3),0d0,kind=8)
        A(1,2) = -A(2,1)
        A(2,2) = conjg(A(1,1))
          
        W(:,N:(N+1)) = matmul(W(:,N:(N+1)),A)
         
      end if  
      
      ! G2 = G2^-1
      G2(2) = -G2(2)
      G2(3) = -G2(3)
 
      ! merge bulge 
      call z_upr1fact_mergebulge(.FALSE.,N,P,Q,D1(1:(2*N)),G2)
      
    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! iteration for QR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  else
  
    ! update V
    if (VEC) then
      
      A(1,1) = cmplx(G2(1),G2(2),kind=8)
      A(2,1) = cmplx(G2(3),0d0,kind=8)
      A(1,2) = -A(2,1)
      A(2,2) = conjg(A(1,1))
      
      V(:,1:2) = matmul(V(:,1:2),A)
      
    end if
    
    ! initialize turnover 
    ! hess
    if (.NOT.P(1)) then
    
      ! set G1 as G2^-1 
      G1(1) = G2(1)
      G1(2) = -G2(2)
      G1(3) = -G2(3)
    
      ! merge from left
      call z_upr1fact_mergebulge(.TRUE.,N,P,Q,D1(1:(2*N)),G1)
      
      ! set G1 for turnover
      G1 = Q(1:3)

      ! pass G2 through triangular part
      call z_upr1fact_rot3throughtri(.FALSE.,D1(1:4),C1(1:6),B1(1:6),G2)
    
      ! set G3 for turnover
      G3 = G2

      ! set G2 for turnover
      G2 = Q(4:6)

    ! inverse hess
    else
    
      ! set G1 as G2^-1 
      G1(1) = G2(1)
      G1(2) = -G2(2)
      G1(3) = -G2(3)
    
      ! pass G2 through triangular part
      call z_upr1fact_rot3throughtri(.FALSE.,D1(1:4),C1(1:6),B1(1:6),G2)
    
      ! merge from right
      call z_upr1fact_mergebulge(.TRUE.,N,P,Q,D1(1:(2*N)),G2)
      
      ! set G3 for turnover
      G3 = Q(1:3)

      ! set G2 for turnover
      G2 = Q(4:6)

    end if
    
    ! chase bulge
    do ii=1,(N-3)
    
      ! execute turnover of G1*G2*G3
      call z_rot3_turnover(G1,G2,G3)
      
      ! prepare for next turnover based on P(ii+1)
      ! hess
      if (.NOT.P(ii+1)) then
      
        ! set P(ii)
        P(ii) = P(ii+1)
        
        ! set Q(ii)
        Q(3*ii-2) = G1(1)
        Q(3*ii-1) = G1(2)
        Q(3*ii) = G1(3)
      
        ! set Q(ii+1)
        Q(3*ii+1) = G2(1)
        Q(3*ii+2) = G2(2)
        Q(3*ii+3) = G2(3)        
 
        ! set G1 for turnover
        G1 = G2     
        
        ! set G2 for turnover
        G2(1) = Q(3*ii+4)
        G2(2) = Q(3*ii+5)
        G2(3) = Q(3*ii+6)
        
        ! update V
        if (VEC) then
          
          A(1,1) = cmplx(G3(1),G3(2),kind=8)
          A(2,1) = cmplx(G3(3),0d0,kind=8)
          A(1,2) = -A(2,1)
          A(2,2) = conjg(A(1,1))
          
          V(:,(ii+1):(ii+2)) = matmul(V(:,(ii+1):(ii+2)),A)
          
        end if
        
        ! pass G3 through upper triangular part
        call z_upr1fact_rot3throughtri(.FALSE.,D1((2*ii+1):(2*ii+4)) &
        ,C1((3*ii+1):(3*ii+6)),B1((3*ii+1):(3*ii+6)),G3)
        
      ! inverse hess
      else
      
        ! set P(ii)
        P(ii) = P(ii+1)
        
        ! set Q(ii)
        Q(3*ii-2) = G1(1)
        Q(3*ii-1) = G1(2)
        Q(3*ii) = G1(3)
      
        ! set Q(ii+1)
        Q(3*ii+1) = G3(1)
        Q(3*ii+2) = G3(2)
        Q(3*ii+3) = G3(3)  
        
        ! set G3 for turnover
        G3 = G3    
        
        ! pass G2 through upper triangular part
        call z_upr1fact_rot3throughtri(.TRUE.,D1((2*ii+1):(2*ii+4)) &
        ,C1((3*ii+1):(3*ii+6)),B1((3*ii+1):(3*ii+6)),G2)
        
        ! update V
        if (VEC) then
          
          A(1,1) = cmplx(G2(1),-G2(2),kind=8)
          A(2,1) = cmplx(-G2(3),0d0,kind=8)
          A(1,2) = -A(2,1)
          A(2,2) = conjg(A(1,1))
          
          V(:,(ii+1):(ii+2)) = matmul(V(:,(ii+1):(ii+2)),A)
          
        end if

        ! set G1 for turnover
        G1 = G2

        ! set G2 for turnover
        G2(1) = Q(3*ii+4)
        G2(2) = Q(3*ii+5)
        G2(3) = Q(3*ii+6)
        
      end if

    end do
    
    ! final turnover
    call z_rot3_turnover(G1,G2,G3)
      
    ! set P(N-1)
    P(N-2) = final_flag
    
    ! finish transformation based on P(N-1)
    ! hess
    if (.NOT.P(N-2)) then
    
      ! set Q(N-2)
      Q(3*(N-2)-2) = G1(1)
      Q(3*(N-2)-1) = G1(2)
      Q(3*(N-2)) = G1(3)   
      
      ! set Q(N-1)
      Q(3*(N-1)-2) = G2(1)
      Q(3*(N-1)-1) = G2(2)
      Q(3*(N-1)) = G2(3)  
      
      ! update V
      if (VEC) then
          
        A(1,1) = cmplx(G3(1),G3(2),kind=8)
        A(2,1) = cmplx(G3(3),0d0,kind=8)
        A(1,2) = -A(2,1)
        A(2,2) = conjg(A(1,1))
          
        V(:,(N-1):N) = matmul(V(:,(N-1):N),A)
         
      end if    
    
      ! pass G3 through upper triangular part
      call z_upr1fact_rot3throughtri(.FALSE.,D1((2*N-3):(2*N)) &
      ,C1((3*N-5):(3*N)),B1((3*N-5):(3*N)),G3)
        
      ! merge bulge 
      call z_upr1fact_mergebulge(.FALSE.,N,P,Q,D1(1:(2*N)),G3)
      
    ! inverse hess
    else
    
      ! set Q(N-2)
      Q(3*(N-2)-2) = G1(1)
      Q(3*(N-2)-1) = G1(2)
      Q(3*(N-2)) = G1(3)   
      
      ! set Q(N-1)
      Q(3*(N-1)-2) = G3(1)
      Q(3*(N-1)-1) = G3(2)
      Q(3*(N-1)) = G3(3)  
      
      ! pass G2 through upper triangular part
      call z_upr1fact_rot3throughtri(.TRUE.,D1((2*N-3):(2*N)) &
      ,C1((3*N-5):(3*N)),B1((3*N-5):(3*N)),G2)
        
      ! update V
      if (VEC) then
          
        A(1,1) = cmplx(G2(1),-G2(2),kind=8)
        A(2,1) = cmplx(-G2(3),0d0,kind=8)
        A(1,2) = -A(2,1)
        A(2,2) = conjg(A(1,1))
          
        V(:,N:(N+1)) = matmul(V(:,N:(N+1)),A)
         
      end if  
      
      ! merge bulge 
      call z_upr1fact_mergebulge(.FALSE.,N,P,Q,D1(1:(2*N)),G2)
      
    end if
    
  end if
  
end subroutine z_upr1fact_singlestep
