#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_rot3_fusion
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_rot3_fusion (swap rotation and diagonal). 
! The following tests are run (once with .TRUE. and once with .FALSE.):
!
! 1)  G1 = normalized([1, 2, 3])
!     G2 = normalized([3, 2, 1])
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rot3_fusion

  implicit none

  ! parameter
  real(8) :: tol = 2d0*EISCOR_DBL_EPS ! accuracy (tolerance)

  ! compute variables
  real(8) :: G1(3), G2(3), nrm

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
  call z_rot3_vec3gen(1d0,2d0,3d0,G1(1),G1(2),G1(3),nrm)
  call z_rot3_vec3gen(3d0,2d0,1d0,G2(1),G2(2),G2(3),nrm)

  ! call fusion
  call z_rot3_fusion(.TRUE.,G1,G2)

  ! check G1
  if (abs(G1(1)+4d0/sqrt(5684d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(G1(2)-48d0/sqrt(5684d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(G1(3)-58d0/sqrt(5684d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check G2
  if (abs(G2(1)-10d0/sqrt(116d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(G2(2)-4d0/sqrt(116d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (G2(3).NE.0d0) then
     call u_test_failed(__LINE__)
  end if

  ! set variables
  call z_rot3_vec3gen(1d0,2d0,3d0,G1(1),G1(2),G1(3),nrm)
  call z_rot3_vec3gen(3d0,2d0,1d0,G2(1),G2(2),G2(3),nrm)

  ! call fusion
  call z_rot3_fusion(.FALSE.,G1,G2)

  ! check G1
  if (abs(G1(1)-10d0/sqrt(116d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(G1(2)+4d0/sqrt(116d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (G1(3).NE.0d0) then
     call u_test_failed(__LINE__)
  end if

  ! check G2
  if (abs(G2(1)+36d0/sqrt(5684d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(G2(2)-32d0/sqrt(5684d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(G2(3)-58d0/sqrt(5684d0))>tol) then
     call u_test_failed(__LINE__)
  end if


  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rot3_fusion
