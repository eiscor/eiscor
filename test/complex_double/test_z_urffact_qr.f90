#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_urffact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_urffact_qr. The following tests are run:
!
! 1) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_urffact_qr

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2**13
  real(8), parameter :: tol = (EISCOR_DBL_EPS)
  integer :: ii, jj, INFO, ITCNT(N-1)
  complex(8) :: U(N)
  real(8) :: VV(N), A(N)
  integer :: id
  real(8) :: small, twopi
  complex(8) :: swap, E(N)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  
 
  ! Check 1)
  ! initialize U and VV
  U = cmplx(0d0,0d0,kind=8)
  U(N) = cmplx(sign(1d0,(-1d0)**(N-1)),0d0,kind=8)
  VV = 1d0
  VV(N) = 0d0
    
  ! call singlestep
  call z_urffact_rfqr(N,U,VV,ITCNT,INFO)

  ! compute by argument
  do ii = 1,N
    call z_scalar_argument(dble(U(ii)),aimag(U(ii)),A(ii),INFO)
  end do

  ! sort by argument
  do ii = 1,N
    small = 10d0
    id = ii
    do jj = ii,N
      if ( A(jj) < small ) then
        id = jj
        small = A(id)
      end if
    end do
    A(id) = A(ii)
    A(ii) = small
    swap = U(id)
    U(id) = U(ii)
    U(ii) = swap
  end do
    
  ! true eigenvalues
  if ( abs(U(1)-cmplx(1d0,0d0,kind=8)) < dble(N)*(EISCOR_DBL_EPS) ) then
    id = 1
  else
    id = 0
  end if
  twopi = 2d0*(EISCOR_DBL_PI)
  do ii = 1,N
    small = twopi*dble(ii-id)/dble(N)
    E(ii) = cmplx(cos(small),sin(small),kind=8)
  end do

  ! compute maximum forward error
  small = 0d0
  do ii = 1,N
    if (abs(U(ii)-E(ii)) > small) then
      small = abs(U(ii)-E(ii))
    end if
  end do

  ! check maximum error
  if (small >= dble(N)*tol) then
    call u_test_failed(__LINE__)
  end if  

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_urffact_qr
