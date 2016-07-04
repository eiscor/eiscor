#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1utri_decompress
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1utri_decompress. 
! The following tests are run (both for row and column scaling):
!
! 1)  D = 1, C^-1 = B = [0+0i,1] 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1utri_decompress

  implicit none

  ! compute variables
  integer, parameter :: N = 2**2
  real(8), parameter :: tol = EISCOR_DBL_EPS
  real(8) :: D(2*N), C(3*N), B(3*N)
  complex(8) :: T(N,N), Ttrue(N,N)
  integer :: ii

  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)




  !!!!!!!!!!!!!!!!!!!!
  ! check 1)
 
  ! set variables
  D = 0d0
  do ii = 1,N
    D(2*ii-1) = 1d0
  end do

  B = 0d0 
  C = 0d0
  do ii = 1,N
    C(3*ii) = -1d0
    B(3*ii) = 1d0
  end do
  
  ! set Ttrue
  Ttrue = cmplx(0d0,0d0,kind=8)
  do ii = 1,N
    Ttrue(ii,1) = cmplx(1d0,0d0,kind=8)
  end do
  
  ! call decompress
  call z_upr1utri_decompress(.TRUE.,N,D,C,B,T)
  
  ! check T
  if ( maxval(dble(abs(T(:,1)-Ttrue(:,1)))) >= dble(N)*tol ) then
    call u_test_failed(__LINE__)
  end if

  ! end check 1





  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
 
  ! set variables
  D = 0d0
  do ii = 1,N
    D(2*ii-1) = 1d0
  end do

  B = 0d0 
  C = 0d0
  do ii = 1,N
    C(3*ii) = -1d0
    B(3*ii) = 1d0
  end do
  
  ! set Ttrue
  Ttrue = cmplx(0d0,0d0,kind=8)
  do ii = 1,N
    Ttrue(ii,ii) = cmplx(1d0,0d0,kind=8)
  end do
  
  ! call decompress
  call z_upr1utri_decompress(.FALSE.,N,D,C,B,T)
  
  ! check T
  if ( maxval(dble(abs(T-Ttrue))) >= dble(N)*tol ) then
    call u_test_failed(__LINE__)
  end if

  ! end check 2





  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_upr1utri_decompress
