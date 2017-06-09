#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_orthfact_2x2diagblock
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_d_orthfact_2x2diagblock.
! The following tests are run:
!
! 1) Q1=[0.8, 0.6], Q2=[0.5,sqrt(3)/2], Q3=[-1/sqrt(2),1/sqrt(2)],
!    D = diag([1, -1, -1, 1])
!    The blocks (1:2,1:2), (2:3,2:3), and (3:4,3:4) are checked.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_orthfact_2x2diagblock

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*EISCOR_DBL_EPS ! accuracy (tolerance)
  
  ! compute variables
  real(8) :: Q(6) ! N = 4
  real(8) :: D(4)
  real(8) :: H(2,2)
  
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
  Q(3) = 5d-1
  Q(4) = -sqrt(3d0)/2d0
  Q(5) = -1d0/sqrt(2d0)
  Q(6) = 1d0/sqrt(2d0)
  
  D(1) = 1d0
  D(2) = -1d0
  D(3) = -1d0
  D(4) = 1d0

  call d_orthfact_2x2diagblock(.TRUE.,Q(1:4),D(1:2),H)

  ! check H
  if (abs(H(1,1)-8d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-6d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)-3d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)+4d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  
  call d_orthfact_2x2diagblock(.FALSE.,Q(3:6),D(3:4),H)

  ! check H
  if (abs(H(1,1)-5d-1/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)+1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+5d-1/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)+1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_d_orthfact_2x2diagblock
