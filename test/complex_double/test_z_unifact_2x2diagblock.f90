#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_unifact_2x2diagblock
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_z_unifact_2x2diagblock.
! The following tests are run:
!
! 1) Q1=[0.8, 0.6, 0.0], Q2=[0.3,0.4,-sqrt(3)/2], 
!    Q3=[-1/sqrt(3),1/sqrt(3), 1/sqrt(3)],
!    D = diag([1, -1, -1, 1])
!    The top and bottom blocks are checked.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_unifact_2x2diagblock

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*EISCOR_DBL_EPS ! accuracy (tolerance)
  
  ! compute variables
  integer :: info
  real(8) :: Q(9) ! N = 4
  real(8) :: D(8)
  complex(8) :: H(2,2)
  
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
  Q(3) = 0d0
  Q(4) = 3d-1
  Q(5) = 4d-1
  Q(6) = -sqrt(3d0)/2d0
  Q(7) = -1d0/sqrt(3d0)
  Q(8) = 1d0/sqrt(3d0)
  Q(9) = 1d0/sqrt(3d0)
  
  D(1) = 1d0
  D(2) = 0d0
  D(3) = -1d0
  D(4) = 0d0
  D(5) = -1d0
  D(6) = 0d0
  D(7) = 1d0
  D(8) = 0d0

  call z_unifact_2x2diagblock(.TRUE.,Q(1:6),D(1:4),H)

  ! check H
  if (abs(dble(H(1,1))-8d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(1,1))-6d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(dble(H(2,1)))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(2,1)))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(dble(H(1,2)))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(1,2)))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(dble(H(2,2))+0.48d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(2,2))+0.14d0)>tol) then
     call u_test_failed(__LINE__)
  end if

  call z_unifact_2x2diagblock(.FALSE.,Q(4:9),D(5:8),H)

  ! check H
  if (abs(dble(H(1,1))+0.1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(1,1))+0.7d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(dble(H(2,1))+1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(2,1)))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(dble(H(1,2))+0.3d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(1,2))-0.4d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(dble(H(2,2))+1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(2,2))+1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_unifact_2x2diagblock
