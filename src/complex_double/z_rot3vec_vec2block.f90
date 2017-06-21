#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3vec_vec2block(ID,N,Q,M,Z)

  implicit none

  ! input variables
  logical, intent(in) :: ID
  integer, intent(in) :: N,M
  real(8), intent(in) :: Q(3*N*N)
  complex(8), intent(inout) :: Z(M,2*N)

  ! compute variables
  integer :: ii, jj, zind, qind
  complex(8) :: A(2,2)

  ! check ID
  if ( ID ) then 
    Z = cmplx(0d0,0d0,kind=8)
    do ii = 1,min(M,2*N)
      Z(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if

  ! fill up Z
  do ii = 1,N

    ! apply rotations
    do jj = 1,N
    
      ! current index in rotation vector
      qind = N*(ii-1) + jj    
 
      ! current column in Z
      zind = N + (ii-1) - (jj-1)

      ! build 2x2 matrix
      A(1,1) = cmplx(Q(3*qind-2),Q(3*qind-1),kind=8)
      A(2,2) = cmplx(Q(3*qind-2),-Q(3*qind-1),kind=8)
      A(2,1) = cmplx(Q(3*qind),0d0,kind=8)
      A(1,2) = cmplx(-Q(3*qind),0d0,kind=8)

      ! multiply into Z
      Z(:,zind:zind+1) = matmul(Z(:,zind:zind+1),A)

    end do

  end do

  
end subroutine z_rot3vec_vec2block
