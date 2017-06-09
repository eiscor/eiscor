#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_rot2_fuse
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine d_rot2_fuse (fuse rotation). 
! The following tests are run (once with .TRUE. and once with .FALSE.):
!
! 1)  fuse([2/sqrt(5),1/sqrt(5)], [2/sqrt(5),-1/sqrt(5)])
! 2)  fuse([-2/sqrt(5),1/sqrt(5)], [2/sqrt(5),1/sqrt(5)])
! 3)  fuse([0.6, 0.8], [0.8,0.6])
! 4)  fuse([0.6, 0.8], [-0.8,-0.6])
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_rot2_fuse

  implicit none

  ! parameter
  real(8) :: tol = 2d0*EISCOR_DBL_EPS ! accuracy (tolerance)

  ! compute variables
  logical :: DIR
  real(8) :: A(2), B(2), C(2), D(2), e, f, nrm

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
  e = 2d0
  f = 1d0
  call d_rot2_vec2gen(e,f,B(1),B(2),nrm)
  e = 2d0
  f = -1d0
  call d_rot2_vec2gen(e,f,A(1),A(2),nrm)
  C = A
  D = B

  DIR = .FALSE.
  call d_rot2_fuse(DIR,A,B)

  DIR = .TRUE.
  call d_rot2_fuse(DIR,C,D)

  ! check result
  if (maxval(abs(B-C))>tol) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  e = -2d0
  f = 1d0
  call d_rot2_vec2gen(e,f,B(1),B(2),nrm)
  e = 2d0
  f = 1d0
  call d_rot2_vec2gen(e,f,A(1),A(2),nrm)
  C = A
  D = B

  DIR = .FALSE.
  call d_rot2_fuse(DIR,A,B)

  DIR = .TRUE.
  call d_rot2_fuse(DIR,C,D)

  ! check result
  if (maxval(abs(B-C))>tol) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  e = 3d0
  f = 4d0
  call d_rot2_vec2gen(e,f,B(1),B(2),nrm)
  e = 4d0
  f = 3d0
  call d_rot2_vec2gen(e,f,A(1),A(2),nrm)
  C = A
  D = B

  DIR = .FALSE.
  call d_rot2_fuse(DIR,A,B)

  DIR = .TRUE.
  call d_rot2_fuse(DIR,C,D)

  ! check result
  if (maxval(abs(B-C))>tol) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  e = 3d0
  f = 4d0
  call d_rot2_vec2gen(e,f,B(1),B(2),nrm)
  e = -4d0
  f = -3d0
  call d_rot2_vec2gen(e,f,A(1),A(2),nrm)
  C = A
  D = B

  DIR = .FALSE.
  call d_rot2_fuse(DIR,A,B)

  DIR = .TRUE.
  call d_rot2_fuse(DIR,C,D)

  ! check result
  if (maxval(abs(B-C))>tol) then
     call u_test_failed(__LINE__)
  end if
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
end program test_d_rot2_fuse
