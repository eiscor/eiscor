#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_orthfact_real2complex
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine d_orthfact_real2complex. The 
! following tests are run:
!
! 1) Convert real Schur form for Nth roots of unity to complex Schur
!    form.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_orthfact_real2complex

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2*8
  real(8), parameter :: PI = EISCOR_DBL_PI 
  real(8), parameter :: tol = 1d1*EISCOR_DBL_EPS 
  integer :: ii, INFO
  real(8) :: Q(2*(N-1)), D(N), Z(N,N) 
  real(8) :: theta
  complex(8) :: b(2,2), t1(2,2), t2(2,2), E(N), V(N,N)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)

  ! initialize INFO
  INFO = 0

  ! initialize b
  b(1,1) = cmplx(0d0,1d0,kind=8)
  b(1,2) = cmplx(-1d0,0d0,kind=8)
  b(2,1) = cmplx(1d0,0d0,kind=8)
  b(2,2) = cmplx(0d0,-1d0,kind=8)
  b = b/sqrt(2d0)

  ! initialize Q
  Q = 0d0
  do ii=1,(N/2)-1
    theta = 2d0*PI*dble(ii-5d-1)/dble(N)
    Q(4*(ii-1)+1) = cos(theta)
    Q(4*(ii-1)+2) = sin(theta)
    Q(4*(ii-1)+3) = 1d0
  end do

  ! initialize D
  D = 1d0

  ! call d_orthfact_real2complex
  call d_orthfact_real2complex(.FALSE.,N,Q,D,N,Z,E,V,INFO)

  ! check INFO
  if (INFO.NE.0) then
    call u_test_failed(__LINE__)
  end if
  
  ! check residuals
  do ii=1,(N/2)
    t1(1,1) = cmplx(Q(4*(ii-1)+1),0d0,kind=8)
    t1(2,1) = cmplx(Q(4*(ii-1)+2),0d0,kind=8)
    t1(1,2) = cmplx(-Q(4*(ii-1)+2),0d0,kind=8)
    t1(2,2) = cmplx(Q(4*(ii-1)+1),0d0,kind=8)
    t1 = matmul(t1,b)

    t2 = b
    t2(:,1) = t2(:,1)*E(2*ii-1)
    t2(:,2) = t2(:,2)*E(2*ii)

    if (maxval(abs(t1-t2)) > tol) then
      call u_test_failed(__LINE__)
    end if

  end do

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_d_orthfact_real2complex
