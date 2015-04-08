#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_rot2_swapdiag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine d_rot2_swapdiag (swap rotation and diagonal). 
!
! 1)  D = diag([ 1  1]), B = [2/sign(5), 1/sign(5)]!
! 2)  D = diag([-1  1]), B = [1/sign(5), -2/sign(5)]
! 3)  D = diag([ 1 -1]), B = [-3/5, 4/5]
! 4)  D = diag([-1 -1]), B = [0,1]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_rot2_swapdiag

  implicit none

  ! parameter
  real(8) :: tol = 2d0*EISCOR_DBL_EPS ! accuracy (tolerance)

  ! compute variables
  real(8) :: H(2,2), D(2), B(2), a, c, nrm

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
  D(1) = 1d0
  D(2) = 1d0
  a = 2d0
  c = 1d0
  call d_rot2_vec2gen(a,c,B(1),B(2),nrm)

  ! compute H
  H(1,1) = B(1)*D(1)
  H(2,1) = B(2)*D(1)
  H(1,2) = -B(2)*D(2)
  H(2,2) = B(1)*D(2)

  call d_rot2_swapdiag(D,B)

  if (abs(H(1,1)-B(1)*D(1))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(2)*D(2))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(2)*D(1))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-B(1)*D(2))>tol) then
     call u_test_failed(__LINE__)
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  D(1) = -1d0
  D(2) = 1d0
  a = 1d0
  c = -2d0
  call d_rot2_vec2gen(a,c,B(1),B(2),nrm)

  ! compute H
  H(1,1) = B(1)*D(1)
  H(2,1) = B(2)*D(1)
  H(1,2) = -B(2)*D(2)
  H(2,2) = B(1)*D(2)

  call d_rot2_swapdiag(D,B)

  if (abs(H(1,1)-B(1)*D(1))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(2)*D(2))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(2)*D(1))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-B(1)*D(2))>tol) then
     call u_test_failed(__LINE__)
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  D(1) = 1d0
  D(2) = -1d0
  a = -3d0
  c = 4d0
  call d_rot2_vec2gen(a,c,B(1),B(2),nrm)

  ! compute H
  H(1,1) = B(1)*D(1)
  H(2,1) = B(2)*D(1)
  H(1,2) = -B(2)*D(2)
  H(2,2) = B(1)*D(2)

  call d_rot2_swapdiag(D,B)

  if (abs(H(1,1)-B(1)*D(1))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(2)*D(2))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(2)*D(1))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-B(1)*D(2))>tol) then
     call u_test_failed(__LINE__)
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  D(1) = -1d0
  D(2) = -1d0
  a = 0d0
  c = 1d0
  call d_rot2_vec2gen(a,c,B(1),B(2),nrm)

  ! compute H
  H(1,1) = B(1)*D(1)
  H(2,1) = B(2)*D(1)
  H(1,2) = -B(2)*D(2)
  H(2,2) = B(1)*D(2)

  call d_rot2_swapdiag(D,B)

  if (abs(H(1,1)-B(1)*D(1))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(2)*D(2))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(2)*D(1))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-B(1)*D(2))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_d_rot2_swapdiag
