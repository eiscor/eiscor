#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_orthfact_2x2deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_d_orthfact_2x2deflation.
! The following tests are run:
!
! 1) Q = [0.8, 0.6], D = diag([1, -1])
! 2) Q = [0.8, 0.6], D = diag([-1, -1])
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_orthfact_2x2deflation

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*EISCOR_DBL_EPS ! accuracy (tolerance)
  
  ! compute variables
  integer, parameter :: M = 2
  real(8) :: Q(2), D(2), Z(M,2)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 1)
  Q(1) = 8d-1
  Q(2) = 6d-1
  
  D(1) = 1d0
  D(2) = -1d0
  
  call d_orthfact_2x2deflation(.FALSE.,Q,D,M,Z)

  ! check Q
  if (abs(Q(1)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(2))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check D
  if (abs(D(1)+1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(2)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  Q(1) = 8d-1
  Q(2) = 6d-1
  
  D(1) = -1d0
  D(2) = -1d0
  
  call d_orthfact_2x2deflation(.FALSE.,Q,D,M,Z)

  ! check Q
  if (abs(Q(1)+8d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)+6d-1)>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check D
  if (abs(D(1)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(2)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_d_orthfact_2x2deflation
