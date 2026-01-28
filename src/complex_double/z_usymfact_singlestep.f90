#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_usymfact_singlestep 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift algorithm on a
! unitary upper hessenberg matrix that is stored as a product of givens
! rotations. 
!                                                                               
! | nu  0 | | u1       -v1 |
! |  0  1 | | v1  conj(u1) | | u2       -v2 | 
!                            | v2  conj(u2) | | u3       -v3 | | 1   0 |
!                                             | v3  conj(u3) | | 0  u4 | 
!
! The square root free algorithm only requires the storage of the vi^2,
! so the arrays U and V contain the following:
!
!  U(i) = ui
!  V(i) = vi
!
! The input must satisfy the following:
!
!  |U(i)|^2 + V(i)^2 = 1
!               V(N) = 0
!               |NU| = 1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER 
!                    dimension of matrix, must be >= 2
!
!  U               COMPLEX(8) array of dimension N
!                    array of complex generators for Givens rotations
!
!  V               REAL(8) array of dimension N
!                    array of real generators for Givens rotations
!
!  NU              COMPLEX(8)
!                    unimodular phase 
!
!  ITCNT           INTEGER
!                    Contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_usymfact_singlestep(N,U,V,NU,ITCNT)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  integer, intent(inout) :: ITCNT
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: V(N)
  complex(8), intent(in) :: NU
  
  ! compute variables
  integer :: ii
  real(8) :: vt, c, s, xx
  complex(8) :: ut, w, rho
  complex(8) :: block(2,2), t1(2,2), t2(2,2)

  ! get 2x2 block
  if (N.EQ.2) then
    block(1,1) =  NU*U(N-1)
    block(2,2) =  conjg(NU*U(N-1))
  else 
    block(1,1) =  U(N-1)
    block(2,2) =  conjg(U(N-1))
  end if
  block(1,2) = -V(N-1)
  block(2,1) =  V(N-1)
  block(:,2) =  block(:,2)*U(N)
  if (N > 2) then
    xx = abs(U(N-2))
    if (xx.GT.0) then
      block(1,:) = conjg(U(N-2))*block(1,:)/xx
    end if
  end if
    
  ! compute eigenvalues and eigenvectors
  t1 = block
  call z_2x2array_eig(.FALSE.,t1,t1,t2,t2)
    
  ! choose Wilkinson shift
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
  ! Wilkinson shift
  else
    rho = rho/xx
  end if

  ! initialize
  w = -rho
  c = 1d0
  s = 0d0

  ! main chasing loop
  do ii=0,(N-1)

    ! set ut and vt
    ut = NU*U(ii+1)
    vt = V(ii+1)

    ! turnover
    call z_rot3_symturnover(w,c,s,ut,vt,rho)

    ! store ut and vt
    if ( ii > 0 ) then
      U(ii) = conjg(NU)*ut
      V(ii) = vt
    end if
    
  end do

end subroutine z_usymfact_singlestep
