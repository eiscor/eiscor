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
!  D               REAL(8) arrays of dimension (2*(N+1))
!                    array of generators for complex diagonal matrices
!
!  C,B             REAL(8) arrays of dimension (3*N)
!                    array of generators for upper-triangular parts of the pencil
!
!  M               INTEGER
!                    leading dimesnion of V and W
!
!  V               COMPLEX(8) array of dimension (M,N)
!                    right schur vectors
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. update V to store right schurvectors 
!
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_singlestep(VEC,FUN,N,P,Q,D,C,B,M,V,ITCNT)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: M, N
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N), C(3*N), B(3*N)
  complex(8), intent(inout) :: V(M,N)
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
  real(8) :: MISFIT(3), G1(3), G2(3), G3(3)
  complex(8) :: shift, A(2,2)
  
  ! compute final_flag
  if (N.LT.3) then 
    final_flag = .FALSE.
  else
    final_flag = FUN(N,P)
  end if  

  ! initialize core chasing
  call z_upr1fact_startchase(VEC,N,P,Q,D,C,B,M,V,ITCNT,MISFIT)
  
  ! core chasing loop
!  do ii=1,(N-3)
!  
!    ! execute turnover of G1*G2*G3
!    call z_rot3_turnover(G1,G2,G3)
!    
!    ! prepare for next turnover based on P(ii+1)
!    ! hess
!    if (.NOT.P(ii+1)) then
!    
!      ! set P(ii)
!      P(ii) = P(ii+1)
!      
!      ! set Q(ii)
!      Q(3*ii-2) = G1(1)
!      Q(3*ii-1) = G1(2)
!      Q(3*ii) = G1(3)
!    
!      ! set Q(ii+1)
!      Q(3*ii+1) = G2(1)
!      Q(3*ii+2) = G2(2)
!      Q(3*ii+3) = G2(3)        
!
!      ! set G1 for turnover
!      G1 = G2     
!      
!      ! set G2 for turnover
!      G2(1) = Q(3*ii+4)
!      G2(2) = Q(3*ii+5)
!      G2(3) = Q(3*ii+6)
!      
!      ! update V
!      if (VEC) then
!        
!        A(1,1) = cmplx(G3(1),G3(2),kind=8)
!        A(2,1) = cmplx(G3(3),0d0,kind=8)
!        A(1,2) = -A(2,1)
!        A(2,2) = conjg(A(1,1))
!        
!        V(:,(ii+1):(ii+2)) = matmul(V(:,(ii+1):(ii+2)),A)
!        
!      end if
!      
!      ! pass G3 through upper triangular part
!      call z_upr1fact_rot3throughtri(.FALSE.,D1((2*ii+1):(2*ii+4)) &
!      ,C1((3*ii+1):(3*ii+6)),B1((3*ii+1):(3*ii+6)),G3)
!      
!    ! inverse hess
!    else
!    
!      ! set P(ii)
!      P(ii) = P(ii+1)
!      
!      ! set Q(ii)
!      Q(3*ii-2) = G1(1)
!      Q(3*ii-1) = G1(2)
!      Q(3*ii) = G1(3)
!    
!      ! set Q(ii+1)
!      Q(3*ii+1) = G3(1)
!      Q(3*ii+2) = G3(2)
!      Q(3*ii+3) = G3(3)  
!      
!      ! set G3 for turnover
!      G3 = G3    
!      
!      ! pass G2 through upper triangular part
!      call z_upr1fact_rot3throughtri(.TRUE.,D1((2*ii+1):(2*ii+4)) &
!      ,C1((3*ii+1):(3*ii+6)),B1((3*ii+1):(3*ii+6)),G2)
!      
!      ! update V
!      if (VEC) then
!        
!        A(1,1) = cmplx(G2(1),-G2(2),kind=8)
!        A(2,1) = cmplx(-G2(3),0d0,kind=8)
!        A(1,2) = -A(2,1)
!        A(2,2) = conjg(A(1,1))
!        
!        V(:,(ii+1):(ii+2)) = matmul(V(:,(ii+1):(ii+2)),A)
!        
!      end if
!
!      ! set G1 for turnover
!      G1 = G2
!
!      ! set G2 for turnover
!      G2(1) = Q(3*ii+4)
!      G2(2) = Q(3*ii+5)
!      G2(3) = Q(3*ii+6)
!      
!    end if
!
!  end do

  ! finish core chasing
  call z_upr1fact_endchase(VEC,N,P,Q,D,C,B,M,V,MISFIT,final_flag)
  
end subroutine z_upr1fact_singlestep
