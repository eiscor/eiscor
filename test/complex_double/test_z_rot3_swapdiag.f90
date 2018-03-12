#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_rot3_swapdiag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_rot3_swapdiag (swap rotation and diagonal). 
! The following tests are run (once with .TRUE. and once with .FALSE.):
!
! 1)  D = diag([ 1  1]), B = normalized([2,1,3])=[2/sign(5), 1/sign(5)]
! 2)  D = diag([-1  1]), B = normalized([1,-2,-1])=[1/sign(5), -2/sign(5)]
! 3)  D = diag([ 1 -1]), B = normalized([-3,4,0])=[-3/5, 4/5, 0]
! 4)  D = diag([-1 -1]), B = [0,0,1]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rot3_swapdiag

  implicit none

  ! parameter
  real(8) :: tol = 2d0*EISCOR_DBL_EPS ! accuracy (tolerance)

  ! compute variables
  real(8) :: D(4), B(3), a, c, e, nrm
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

  ! set variables
  D(1) = 1d0
  D(2) = 0d0
  D(3) = 1d0
  D(4) = 0d0
  a = 2d0
  c = 1d0
  e = 3d0
  call z_rot3_vec3gen(a,c,e,B(1),B(2),B(3),nrm)

  ! compute H
  H(1,1) = cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(3)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(3)*cmplx(D(3),D(4),kind=8)
  H(2,2) = cmplx(B(1),B(2),kind=8)*cmplx(D(3),D(4),kind=8)

  call z_rot3_swapdiag(D,B)

  if (abs(H(1,1)-cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(3)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(3)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-cmplx(B(1),B(2),kind=8)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! set variables
  D(1) = 1d0
  D(2) = 0d0
  D(3) = 1d0
  D(4) = 0d0
  a = 2d0
  c = 1d0
  call z_rot3_vec3gen(a,c,e,B(1),B(2),B(3),nrm)

  ! compute H
  H(1,1) = cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(3)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(3)*cmplx(D(3),D(4),kind=8)
  H(2,2) = cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8)

  call z_rot3_swapdiag(D,B)

  if (abs(H(1,1)-cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(3)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(3)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  D(1) = -1d0
  D(2) = 0d0
  D(3) = 1d0
  D(4) = 0d0
  a = 1d0
  c = -2d0
  e = -1d0
  call z_rot3_vec3gen(a,c,e,B(1),B(2),B(3),nrm)

  ! compute H
  H(1,1) = cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(3)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(3)*cmplx(D(3),D(4),kind=8)
  H(2,2) = cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8)

  call z_rot3_swapdiag(D,B)

  if (abs(H(1,1)-cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(3)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(3)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! set variables
  D(1) = -1d0
  D(2) = 0d0
  D(3) = 1d0
  D(4) = 0d0
  a = 1d0
  c = -2d0
  e = -1d0
  call z_rot3_vec3gen(a,c,e,B(1),B(2),B(3),nrm)

  ! compute H
  H(1,1) = cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(3)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(3)*cmplx(D(3),D(4),kind=8)
  H(2,2) = cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8)

  call z_rot3_swapdiag(D,B)

  if (abs(H(1,1)-cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(3)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(3)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  D(1) = 1d0
  D(2) = 0d0
  D(3) = -1d0
  D(4) = 0d0
  a = -3d0
  c = 4d0
  e = 0d0
  call z_rot3_vec3gen(a,c,e,B(1),B(2),B(3),nrm)

  ! compute H
  H(1,1) = cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(3)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(3)*cmplx(D(3),D(4),kind=8)
  H(2,2) = cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8)

  call z_rot3_swapdiag(D,B)

  if (abs(H(1,1)-cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(3)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(3)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! set variables
  D(1) = 1d0
  D(2) = 0d0
  D(3) = -1d0
  D(4) = 0d0
  a = -3d0
  c = 4d0
  e = 0d0
  call z_rot3_vec3gen(a,c,e,B(1),B(2),B(3),nrm)

  ! compute H
  H(1,1) = cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(3)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(3)*cmplx(D(3),D(4),kind=8)
  H(2,2) = cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8)

  call z_rot3_swapdiag(D,B)

  if (abs(H(1,1)-cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(3)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(3)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  D(1) = -1d0
  D(2) = 0d0
  D(3) = -1d0
  D(4) = 0d0
  a = 0d0
  c = 0d0
  e = 1d0
  call z_rot3_vec3gen(a,c,e,B(1),B(2),B(3),nrm)

  ! compute H
  H(1,1) = cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(3)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(3)*cmplx(D(3),D(4),kind=8)
  H(2,2) = cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8)

  call z_rot3_swapdiag(D,B)

  if (abs(H(1,1)-cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(3)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(3)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! set variables
  D(1) = -1d0
  D(2) = 0d0
  D(3) = -1d0
  D(4) = 0d0
  a = 0d0
  c = 0d0
  e = 1d0
  call z_rot3_vec3gen(a,c,e,B(1),B(2),B(3),nrm)

  ! compute H
  H(1,1) = cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(3)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(3)*cmplx(D(3),D(4),kind=8)
  H(2,2) = cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8)

  call z_rot3_swapdiag(D,B)

  if (abs(H(1,1)-cmplx(B(1),B(2),kind=8)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(3)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(3)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-cmplx(B(1),-B(2),kind=8)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rot3_swapdiag
