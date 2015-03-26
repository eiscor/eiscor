#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_unifact_buildbulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_z_unifact_buildbulge.
! The following tests are run:
!
! Matrix:
! Q1=[0.8, 0.6, 0.0], Q2=[0.3,0.4,-sqrt(3)/2], 
! Q3=[-1/sqrt(3),1/sqrt(3), 1/sqrt(3)]
! D = diag([1, 0.8+i0,6, 1/sqrt(2)+i/sqrt(2), 1])
!
! 1) shift = 0, shift = 0.8-1/sqrt(2) + i(0.6-1/sqrt(2))
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_unifact_buildbulge

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*EISCOR_DBL_EPS ! accuracy (tolerance)
  
  ! compute variables
  real(8) :: Q(9), D(8), B(3), B2(3), nul
  complex(8) :: H(2,2), shift
  
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
  D(3) = 0.8d0
  D(4) = 0.6d0
  D(5) = 1d0/sqrt(2d0)
  D(6) = 1d0/sqrt(2d0)
  D(7) = 1d0
  D(8) = 0d0
  
  SHIFT = cmplx(0d0,0d0,kind=8)

  call z_unifact_buildbulge(Q(1:6),D(1:4),SHIFT,B)

  if (abs(B(1)-0.8d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(2)-0.6d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(3))>tol) then
    call u_test_failed(__LINE__)
  end if

  SHIFT = cmplx(+0.8d0-1d0/sqrt(2d0),+0.6d0-1d0/sqrt(2d0),kind=8)

  call z_unifact_buildbulge(Q(1:6),D(1:4),SHIFT,B)

  if (abs(B(1)-1d0/sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(2)-1d0/sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(3))>tol) then
    call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_unifact_buildbulge
