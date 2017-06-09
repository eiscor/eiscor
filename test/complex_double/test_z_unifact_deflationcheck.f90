#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_unifact_deflationcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_z_unifact_deflationcheck.
! The following tests are run:
!
! 1) Q1=[0.8, 0.6, 0.0], Q2=[0.3,0.4,-sqrt(3)/2], 
!    Q3=[-1/sqrt(3),1/sqrt(3), 1/sqrt(3)]
!    D = diag([1, 0.8+i0,6, 1/sqrt(2)+i/sqrt(2), 1])
!
! 1B) set STR = 2
!
! 2) set Q2 = [5/sqrt(60), -7/sqrt(60), 0], STR = 2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_unifact_deflationcheck

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*EISCOR_DBL_EPS ! accuracy (tolerance)
  
  ! compute variables
  integer :: ZERO
  real(8) :: Q(9) ! N = 4
  real(8) :: D(8)
  
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
  
  ZERO = 0

  call z_unifact_deflationcheck(4,Q,D,ZERO)

  ! check Q
  if (abs(Q(1)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(2))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(4)-48d-2)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(5)-14d-2)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(6)+sqrt(3d0)/2d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(7)+0.2d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(8)-1.4d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(9)-1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check D
  if (abs(D(1)-0.8d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(2)-0.6d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(3)-0.8d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(4)-0.6d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(5)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(6)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(7)-0.8d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(8)+0.6d0)>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check integers
  if (zero.NE.1) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 1B)
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
  
  ZERO = 0

  call z_unifact_deflationcheck(3,Q(4:9),D(3:8),ZERO)

  ! check Q
  if (abs(Q(1)-8d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-6d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(4)-3d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(5)-4d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(6)+sqrt(3d0)/2d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(7)+1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(8)-1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(9)-1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check D
  if (abs(D(1)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(2))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(3)-0.8d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(4)-0.6d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(5)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(6)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(7)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(8))>tol) then
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
  Q(3) = 0d0
  Q(4) = 5d0/sqrt(61d0)
  Q(5) = -6d0/sqrt(61d0)
  Q(6) = 0d0
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
  
  ZERO = 0

  call z_unifact_deflationcheck(3,Q(4:9),D(3:8),ZERO)

  ! check Q
  if (abs(Q(1)-8d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-6d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(4)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(5))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(6))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(7)+11d0/sqrt(183d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(8)+1d0/sqrt(183d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(9)-1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check D
  if (abs(D(1)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(2))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(3)-7.6d0/sqrt(61d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(4)+1.8d0/sqrt(61d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(5)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(6)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(7)-5d0/sqrt(61d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(8)-6d0/sqrt(61d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check integers
  if (zero.NE.1) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_unifact_deflationcheck
