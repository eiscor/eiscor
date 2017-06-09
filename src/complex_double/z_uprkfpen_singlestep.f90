#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkfpen_singlestep 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift 
! algorithm on a factored unitary plus rank k (uprkfpen) matrix.. 
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
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkfpen_singlestep(VECR,VECL,FUN,N,K,STR,STP,P,Q,D1,C1,B1,D2,C2,B2,M,V,W,ITCNT)

  implicit none
  
  ! input variables
  logical, intent(in) :: VECR, VECL
  integer, intent(in) :: M, N, K, STR, STP
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*N*K), C1(3*N*K), B1(3*N*K)
  real(8), intent(inout) :: D2(2*N*K), C2(3*N*K), B2(3*N*K)
  complex(8), intent(inout) :: V(M,N),W(M,N)
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
  real(8) :: MISFIT(3)
  
  ! compute final_flag
  if (N.LT.3) then 
    final_flag = .FALSE.
  else
    final_flag = FUN(N,P)
  end if  

  ! initialize core chasing
  call z_uprkfpen_startchase(VECR,VECL,N,K,STR,STP,P,Q,D1,C1,B1,&
       &D2,C2,B2,M,V(:,STR:STR+1),W(:,STR:STR+1),ITCNT,MISFIT)
  
  ! core chasing loop
  do ii=STR,(STP-2)

    ! move misfit down one row
    call z_uprkfpen_chasedown(VECR,VECL,N,K,ii,P,Q,D1,C1,B1,D2,C2,B2,&
         &M,V(:,ii+1:ii+2),W(:,ii+1:ii+2),MISFIT)

    !call z_upr1fpen_chasedown(VEC,P(ii:ii+1),Q((ir1-3):(ir2-3)),D1(id1:id2), & 
    !                          C1(ir1:ir2),B1(ir1:ir2),D2(id1:id2),C2(ir1:ir2), &
    !                          B2(ir1:ir2),M,V(:,ii+1:ii+2),W(:,ii+1:ii+2),MISFIT)

  end do

  ! finish core chasing
  call z_uprkfpen_endchase(VECR,VECL,N,K,STR,STP,P,Q,D1,C1,B1,D2,C2,B2,M,V(:,STP:STP+1),W(:,STP:STP+1),MISFIT,final_flag)
  
end subroutine z_uprkfpen_singlestep
