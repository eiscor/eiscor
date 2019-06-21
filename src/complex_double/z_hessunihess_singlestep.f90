#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_hessunihess_singlestep 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift 
! algorithm on a factored unitary plus rank one (upr1fpen) matrix.. 
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
!  A               COMPLEX(8) array of dimension (N,N)
!                    array of the dense Hessenberg matrix
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
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
subroutine z_hessunihess_singlestep(VEC,N,A,Q,M,V,W,ITCNT)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: M, N
  real(8), intent(inout) :: Q(3*(N-1))
  complex(8), intent(inout) :: A(N,N),V(M,N),W(M,N)
  integer, intent(in) :: ITCNT
  
  ! compute variables
  integer :: ii, ir1, ir2, id1, id2
  real(8) :: MISFIT(3)
  
  ! initialize core chasing
  !call z_hessunihess_startchase(VEC,N,A,Q,M,V(:,1:2),W(:,1:2),ITCNT,MISFIT)
  
  ! core chasing loop
  do ii=1,(N-3)

    ! compute indices
    ir1 = 3*(ii)+1
    ir2 = 3*(ii+2)
    id1 = 2*(ii)+1
    id2 = 2*(ii+2)

    ! move misfit down one row
    !call z_hessunihess_chasedown(VEC,N,A,Q((ir1-3):(ir2-3)), & 
    !                          M,V(:,ii+1:ii+2),W(:,ii+1:ii+2),MISFIT)

  end do

  ! finish core chasing
  !call z_hessunihess_endchase(VEC,N,A,Q,M,V(:,N-1:N),W(:,N-1:N),MISFIT,final_flag)
  
end subroutine z_hessunihess_singlestep
