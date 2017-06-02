#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_urffact_singlestep 
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
subroutine z_urffact_singlestep(N,U,VV,ITCNT)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: ITCNT
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: VV(N)
  
  ! compute variables
  integer :: ii
  real(8) :: nn, zz, cc, ss, xx
  complex(8) :: z, w, rho
  complex(8) :: block(2,2), t1(2,2), t2(2,2)

  ! get 2x2 block
  block(1,1) =  U(N-1)
  block(2,2) =  conjg(U(N-1))
  block(1,2) = -sqrt(VV(N-1))
  block(2,1) =  sqrt(VV(N-1))
  block(:,2) =  block(:,2)*U(N)
  if (N > 2) then
    block(1,:) = conjg(U(N-2))*block(1,:)
  end if
    
  ! compute eigenvalues and eigenvectors
  t1 = block
  call z_2x2array_eig(.FALSE.,t1,t1,t2,t2)
    
  ! choose wikinson shift
  ! complex abs does not matter here
  if(abs(block(2,2)-t1(1,1)) < abs(block(2,2)-t1(2,2)))then
    rho = t1(1,1)
  else
    rho = t1(2,2)
  end if

  ! compute a nonzero shift
  ! random shift
  xx = abs(rho)
  if (xx == 0) then
    call random_number(xx)
    rho = cmplx(cos(xx),sin(xx),kind=8)
  ! wilkinson shift
  else
    rho = rho/xx
  end if
!print*,""
!print*,"rho"
!print*,rho

  ! initialize
  w = -rho
  z = U(1) + w
  zz = dble(z)**2 + aimag(z)**2
  nn = VV(1) + zz
  cc = 1d0

  ! main chasing loop
  do ii=1,(N-1)
    cc = cc*zz/nn
    ss = VV(ii)/nn
    w = -rho*conjg(w)*(z*z)/zz
    z = U(ii+1) + w
    zz = dble(z)**2 + aimag(z)**2
    nn = VV(ii+1) + cc*zz
    U(ii) = conjg(rho)*(U(ii+1) - cc*z)
    VV(ii)   = ss*nn
    xx = dble(U(ii))**2 + aimag(U(ii))**2 + VV(ii)
    U(ii) = 5d-1*U(ii)*(3d0 - xx)
    VV(ii) = VV(ii)*(2d0 - xx)
  end do
  
end subroutine z_urffact_singlestep
