#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_orthfact_deflationcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_d_orthfact_deflationcheck.
! The following tests are run:
!
! 1) no deflation in 
!    Q1=[0.8, 0.6], Q2=[0.5,sqrt(3)/2], Q3=[-1/sqrt(2),1/sqrt(2)],
!    D = diag([1, -1, -1, 1])
!
! 2) deflation in position 3
!    Q1=[0.8, 0.6], Q2=[0.5,sqrt(3)/2], Q3=[-1,eps/2], 
!    Q4=[-1/sqrt(2),1/sqrt(2)],
!    D = diag([1, -1, -1, 1, 1])
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_orthfact_deflationcheck

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*EISCOR_DBL_EPS ! accuracy (tolerance)
  
  ! compute variables
  integer :: zero
  real(8) :: Q(8) ! N = 5
  real(8) :: D(5)
  
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
  
  call d_orthfact_deflationcheck(4,Q,D,zero)

  ! check Q
  if (abs(Q(1)-8d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-6d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-5d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(4)+sqrt(3d0)/2d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(5)+1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(6)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check D
  if (abs(D(1)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(2)+1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(3)+1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(4)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check integers
  if (zero.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  Q(1) = 8d-1
  Q(2) = 6d-1
  Q(3) = 5d-1
  Q(4) = -sqrt(3d0)/2d0
  Q(5) = -1d0
  Q(6) = tol/2d1
  Q(7) = -1d0/sqrt(2d0)
  Q(8) = 1d0/sqrt(2d0)
  
  D(1) = 1d0
  D(2) = -1d0
  D(3) = -1d0
  D(4) = 1d0
  D(5) = 1d0
  
  call d_orthfact_deflationcheck(5,Q,D,zero)

  ! check Q
  if (abs(Q(1)-8d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-6d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-5d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(4)+sqrt(3d0)/2d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(5)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(6))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(7)+1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(8)+1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check D
  if (abs(D(1)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(2)+1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(3)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(4)+1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(5)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check integers
  if (zero.NE.3) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_d_orthfact_deflationcheck
