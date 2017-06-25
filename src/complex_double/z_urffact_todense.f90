#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_urffact_todense
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine diagonalizes a unitary upper hessenberg matrix that is stored as
! a product of N Givens rotations, without computing square roots.
!
! | u1       -v1 |
! | v1  conj(u1) | | u2       -v2 | 
!                  | v2  conj(u2) | | u3       -v3 | | 1   0 |
!                                   | v3  conj(u3) | | 0  u4 | 
!                                                                     
! The square root free algorithm only requires the storage of the vi^2,
! so the arrays U and VV contain the following:
!
!  U(i) = ui
! VV(i) = vi^2
!
! The input must satisfy the following:
!
!  |U(i)|^2 + VV(i) = 1
!             VV(N) = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  U               COMPLEX(8) array of dimension N
!                    array of complex generators for Givens rotations
!                    on output contains eigenvalues
!
!  VV              REAL(8) array of dimension N
!                    array of real generators (squared) for Givens rotations
!
! OUTPUT VARIABLES:
!
!  ITS             INTEGER array of dimension N-1
!                    contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO = 1 implies no convergence
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies N is invalid
!                    INFO = -2 implies U or VV is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_urffact_todense(N,U,VV,M,H)

  implicit none
  
  ! input variables
  integer, intent(in) :: N,M
  complex(8), intent(in) :: U(N)
  real(8), intent(in) :: VV(N)
  complex(8), intent(inout) :: H(M,N)
  
  ! compute variables
  integer :: ii, low
  complex(8) :: A(2,2)
  
  ! initialize H
  H = cmplx(0d0,0d0,kind=8)
  do ii = 1,min(M,N)
    H(ii,ii) = cmplx(1d0,0d0,kind=8)
  end do
  
  ! fill H
  do ii = 1,N-1

    ! fill A
    A(1,1) = U(ii)
    A(2,1) = sqrt(VV(ii))
    A(2,2) = conjg(A(1,1))
    A(1,2) = -A(2,1)

    ! set low 
    low = min(ii+1,M)

    ! update H
    H(1:low,ii:ii+1) = matmul(H(1:low,ii:ii+1),A)
    
  end do
  
  ! last column
  H(:,N) = H(:,N)*U(N)

end subroutine z_urffact_todense
