#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_urffact_singlestep_shift2 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift algorithm on a
! unitary upper hessenberg matrix that is stored as a product of givens
! rotations. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER 
!                    dimension of matrix, must be >= 2
!
!  U               REAL(8) array of dimension (3*N)
!                    array of generators for givens rotations
!
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_urffact_singlestep_shift2(N,U,V,ITCNT,SHIFT)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: ITCNT
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: V(N)
  complex(8), intent(in) :: SHIFT
  
  ! compute variables
  integer :: ii
  real(8) :: nrm, cc, ss
  complex(8) :: rho, swap, w
  real(8) :: G1(3), G2(3), G3(3)

  rho = SHIFT

  ! initialize
  call z_rot3_vec3gen(dble(U(1)-rho),aimag(U(1)-rho),V(1),G1(1),G1(2),G1(3),nrm)
  G3 = G1

!print*,""
  ! main chasing loop
  do ii=1,(N-1)

    ! set G2
    G2 = (/ dble(U(ii+1)), aimag(U(ii+1)), V(ii+1) /)

w = cmplx(G1(1),G1(2),kind=8)/cmplx(G3(1),G3(2),kind=8)
cc = abs((cmplx(G1(1),G1(2),kind=8)*cmplx(G3(1),G3(2),kind=8))/w)
ss = G1(3)*G3(3)
!print*,w,cc,ss,cc+ss

    ! turnover
    call z_rot3_turnover(G1,G2,G3)

    ! update U and V
    U(ii) = -conjg(rho)*cmplx(G1(1),G1(2),kind=8) 
    V(ii) = G1(3)

    ! set G1
    G1 = G2

    ! update phase of G1
    swap = -rho*cmplx(G1(1),G1(2),kind=8)
    G1(1) = dble(swap)
    G1(2) = aimag(swap)
  end do
  
end subroutine z_urffact_singlestep_shift2
