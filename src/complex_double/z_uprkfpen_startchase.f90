#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkfpen_startchase 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine initializes one iteration of Francis' singleshift 
! algorithm for a factored unitary plus rank one (uprkfpen) matrix pencil.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvector
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
!  V,W             COMPLEX(8) array of dimension (M,N)
!                    right and left schurvectors 
!
! OUTPUT VARIABLES:
!
!  G               REAL(8) array of dimension 3
!                    generators for bulge core transformation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkfpen_startchase(VEC,N,K,STR,STP,P,Q,D1,C1,B1,D2,C2,B2,M,V,W,ITCNT,G)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: M, N, K, STR, STP, ITCNT
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*N*K), C1(3*N*K), B1(3*N*K)
  real(8), intent(inout) :: D2(2*N*K), C2(3*N*K), B2(3*N*K), G(3)
  complex(8), intent(inout) :: V(M,2), W(M,2)
  
  ! compute variables
  integer :: ii, ir1, ir2, id1, id2
  logical :: final_flag, tp(2)
  real(8) :: nrm, Ginv(3)
  real(8) :: tq(6), td1(6), tc1(9), tb1(9)
  real(8) :: td2(6), tc2(9), tb2(9)
  complex(8) :: shift, A(2,2)
  
  ! compute shift
  ! random shift
  if ((mod(ITCNT,15).EQ.0).AND.(ITCNT.GT.0)) then
     call random_number(G(1))
     call random_number(G(2))
     shift = cmplx(G(1),G(2),kind=8)
     print*, "Random shift", shift
          
  ! wilkinson shift
  else
  
    ! compute wilkinson shift
    call z_uprkfpen_singleshift(N,K,STR,STP,P,Q,D1,C1,B1,D2,C2,B2,shift)

  end if

  ! build bulge
  call z_uprkfpen_buildbulge(N,K,STR,STP,P,Q,D1,C1,B1,D2,C2,B2,shift,G)


  ! set Ginv
  Ginv(1) = G(1)
  Ginv(2) = -G(2)
  Ginv(3) = -G(3)
  
  ! update left schurvectors with G
  if (VEC) then
    
    A(1,1) = cmplx(G(1),G(2),kind=8)
    A(2,1) = cmplx(G(3),0d0,kind=8)
    A(1,2) = -A(2,1)
    A(2,2) = conjg(A(1,1))
    
    W = matmul(W,A)
    
  end if

  ! copy Ginv into G
  G = Ginv

  ! initialize turnover 
  if (STP.LE.STR) then

    ! fuse Ginv and Q, Ginv is now a diagonal rotation
    call z_rot3_fusion(.FALSE.,Ginv,Q(3*STR-2:3*STR))

    ! pass G through R2
    call z_uprkutri_rot3swap(.TRUE.,N,K,1,K,D2,C2,B2,G,STR)
    
    ! update left schurvectors diagonal Ginv
    if (VEC) then
      
      W(:,1) = W(:,1)*cmplx(Ginv(1),Ginv(2),kind=8)
      W(:,2) = W(:,2)*cmplx(Ginv(1),-Ginv(2),kind=8)
      
    end if

    ! Ginv scales the rows of R2
    call z_upr1utri_unimodscale(.TRUE.,D2(2*STR-1:2*STR),C2(3*STR-2:3*STR),B2(3*STR-2:3*STR), &
                                cmplx(Ginv(1),-Ginv(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D2(2*STR+1:2*STR+2),C2(3*STR+1:3*STR+3),B2(3*STR+1:3*STR+3), &
                                cmplx(Ginv(1),Ginv(2),kind=8))
    ! true -> does not touch C2 and B2


    ! invert G
    G(2) = -G(2)
    G(3) = -G(3)

    ! update right schurvectors with G
    if (VEC) then
      
      A(1,1) = cmplx(G(1),G(2),kind=8)
      A(2,1) = cmplx(G(3),0d0,kind=8)
      A(1,2) = -A(2,1)
      A(2,2) = conjg(A(1,1))
      
      V = matmul(V,A)
      
    end if

    ! pass G through R1
    call z_uprkutri_rot3swap(.FALSE.,N,K,1,K,D1,C1,B1,G,STR)

  ! hess
  elseif (.NOT.P(STR)) then
  
    ! fuse Ginv and Q, Ginv is now a diagonal rotation
    call z_rot3_fusion(.FALSE.,Ginv,Q(3*STR-2:3*STR))

    ! pass G through R2
    call z_uprkutri_rot3swap(.TRUE.,N,K,1,K,D2,C2,B2,G,STR)
  
    ! update left schurvectors diagonal Ginv
    if (VEC) then
      
      W(:,1) = W(:,1)*cmplx(Ginv(1),Ginv(2),kind=8)
      W(:,2) = W(:,2)*cmplx(Ginv(1),-Ginv(2),kind=8)
      
    end if

    ! Ginv scales the rows of R2
    call z_upr1utri_unimodscale(.TRUE.,D2(2*STR-1:2*STR),C2(3*STR-2:3*STR),&
         &B2(3*STR-2:3*STR), cmplx(Ginv(1),-Ginv(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D2(2*STR+1:2*STR+2),C2(3*STR+1:3*STR+3),&
         &B2(3*STR+1:3*STR+3), cmplx(Ginv(1),Ginv(2),kind=8))
    ! true -> does not touch C2 and B2

    ! invert G
    G(2) = -G(2)
    G(3) = -G(3)

    ! update right schurvectors with G
    if (VEC) then
      
      A(1,1) = cmplx(G(1),G(2),kind=8)
      A(2,1) = cmplx(G(3),0d0,kind=8)
      A(1,2) = -A(2,1)
      A(2,2) = conjg(A(1,1))
      
      V = matmul(V,A)
      
    end if

    ! pass G through R1
    call z_uprkutri_rot3swap(.FALSE.,N,K,1,K,D1,C1,B1,G,STR)

  ! inverse hess
  else
  
    ! pass Ginv through R2
    call z_uprkutri_rot3swap(.TRUE.,N,K,1,K,D2,C2,B2,Ginv,STR)
  
    ! invert Ginv
    Ginv(2) = -Ginv(2)
    Ginv(3) = -Ginv(3)

    ! update right schurvectors with Ginv
    if (VEC) then
      
      A(1,1) = cmplx(Ginv(1),Ginv(2),kind=8)
      A(2,1) = cmplx(Ginv(3),0d0,kind=8)
      A(1,2) = -A(2,1)
      A(2,2) = conjg(A(1,1))
      
      V = matmul(V,A)
      
    end if

    ! pass Ginv through R1
    call z_uprkutri_rot3swap(.FALSE.,N,K,1,K,D1,C1,B1,Ginv,STR)

    ! fuse Ginv and Q, Ginv is now a diagonal rotation
    call z_rot3_fusion(.TRUE.,Q(3*STR-2:3*STR),Ginv)

    ! Ginv scales the rows of R1
    call z_upr1utri_unimodscale(.TRUE.,D1(1:2),C1(1:3),B1(1:3), &
                                cmplx(Ginv(1),Ginv(2),kind=8))
    call z_upr1utri_unimodscale(.TRUE.,D1(3:4),C1(4:6),B1(4:6), &
                                cmplx(Ginv(1),-Ginv(2),kind=8))
    ! C1, B1 not used

  end if
  
end subroutine z_uprkfpen_startchase
