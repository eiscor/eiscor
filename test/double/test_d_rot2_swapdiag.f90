#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_rot2_swapdiag
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine d_rot2_swapdiag (swap rotation and diagonal). 
! The following tests are run (once with 'L' and once with 'R'):
!
! 1)  D = diag([ 1  1]), B = [2/sign(5), 1/sign(5)]!
! 2)  D = diag([-1  1]), B = [1/sign(5), -2/sign(5)]
! 3)  D = diag([ 1 -1]), B = [-3/5, 4/5]
! 4)  D = diag([-1 -1]), B = [0,1]
!
! in DEBUG mode additionally the check of the input data is checked
! A) job='K'
! B) diagonal (D) not real, diagonal not of absolute value 1
! C) INF and NAN in rotation (B)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_rot2_swapdiag

  implicit none

  ! parameter
  real(8) :: tol = 2d0*epsilon(1d0) ! accuracy (tolerance)

  ! compute variables
  character :: job
  real(8) :: D(4), B(2), a, c, nrm
  complex(8) :: H(2,2)
  integer :: info

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
  job = 'L'
  D(1) = 1d0
  D(2) = 0d0
  D(3) = 1d0
  D(4) = 0d0
  a = 2d0
  c = 1d0
  call d_rot2_vec2gen(a,c,B(1),B(2),nrm,info)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! compute H
  H(1,1) = B(1)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(2)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(2)*cmplx(D(3),D(4),kind=8)
  H(2,2) = B(1)*cmplx(D(3),D(4),kind=8)

  call d_rot2_swapdiag(JOB,D,B,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  if (abs(H(1,1)-B(1)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(2)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(2)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-B(1)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! set variables
  job = 'R'
  D(1) = 1d0
  D(2) = 0d0
  D(3) = 1d0
  D(4) = 0d0
  a = 2d0
  c = 1d0
  call d_rot2_vec2gen(a,c,B(1),B(2),nrm,info)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! compute H
  H(1,1) = B(1)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(2)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(2)*cmplx(D(3),D(4),kind=8)
  H(2,2) = B(1)*cmplx(D(3),D(4),kind=8)

  call d_rot2_swapdiag(JOB,D,B,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  if (abs(H(1,1)-B(1)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(2)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(2)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-B(1)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  job = 'L'
  D(1) = -1d0
  D(2) = 0d0
  D(3) = 1d0
  D(4) = 0d0
  a = 1d0
  c = -2d0
  call d_rot2_vec2gen(a,c,B(1),B(2),nrm,info)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! compute H
  H(1,1) = B(1)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(2)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(2)*cmplx(D(3),D(4),kind=8)
  H(2,2) = B(1)*cmplx(D(3),D(4),kind=8)

  call d_rot2_swapdiag(JOB,D,B,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  if (abs(H(1,1)-B(1)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(2)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(2)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-B(1)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! set variables
  job = 'R'
  D(1) = -1d0
  D(2) = 0d0
  D(3) = 1d0
  D(4) = 0d0
  a = 1d0
  c = -2d0
  call d_rot2_vec2gen(a,c,B(1),B(2),nrm,info)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! compute H
  H(1,1) = B(1)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(2)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(2)*cmplx(D(3),D(4),kind=8)
  H(2,2) = B(1)*cmplx(D(3),D(4),kind=8)

  call d_rot2_swapdiag(JOB,D,B,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  if (abs(H(1,1)-B(1)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(2)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(2)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-B(1)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  job = 'L'
  D(1) = 1d0
  D(2) = 0d0
  D(3) = -1d0
  D(4) = 0d0
  a = -3d0
  c = 4d0
  call d_rot2_vec2gen(a,c,B(1),B(2),nrm,info)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! compute H
  H(1,1) = B(1)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(2)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(2)*cmplx(D(3),D(4),kind=8)
  H(2,2) = B(1)*cmplx(D(3),D(4),kind=8)

  call d_rot2_swapdiag(JOB,D,B,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  if (abs(H(1,1)-B(1)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(2)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(2)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-B(1)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! set variables
  job = 'R'
  D(1) = 1d0
  D(2) = 0d0
  D(3) = -1d0
  D(4) = 0d0
  a = -3d0
  c = 4d0
  call d_rot2_vec2gen(a,c,B(1),B(2),nrm,info)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! compute H
  H(1,1) = B(1)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(2)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(2)*cmplx(D(3),D(4),kind=8)
  H(2,2) = B(1)*cmplx(D(3),D(4),kind=8)

  call d_rot2_swapdiag(JOB,D,B,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  if (abs(H(1,1)-B(1)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(2)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(2)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-B(1)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  job = 'L'
  D(1) = -1d0
  D(2) = 0d0
  D(3) = -1d0
  D(4) = 0d0
  a = 0d0
  c = 1d0
  call d_rot2_vec2gen(a,c,B(1),B(2),nrm,info)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! compute H
  H(1,1) = B(1)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(2)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(2)*cmplx(D(3),D(4),kind=8)
  H(2,2) = B(1)*cmplx(D(3),D(4),kind=8)

  call d_rot2_swapdiag(JOB,D,B,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  if (abs(H(1,1)-B(1)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(2)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(2)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-B(1)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! set variables
  job = 'R'
  D(1) = -1d0
  D(2) = 0d0
  D(3) = -1d0
  D(4) = 0d0
  a = 0d0
  c = 1d0
  call d_rot2_vec2gen(a,c,B(1),B(2),nrm,info)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! compute H
  H(1,1) = B(1)*cmplx(D(1),D(2),kind=8)
  H(2,1) = B(2)*cmplx(D(1),D(2),kind=8)
  H(1,2) = -B(2)*cmplx(D(3),D(4),kind=8)
  H(2,2) = B(1)*cmplx(D(3),D(4),kind=8)

  call d_rot2_swapdiag(JOB,D,B,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  if (abs(H(1,1)-B(1)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,1)-B(2)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(1,2)+B(2)*cmplx(D(1),D(2),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(H(2,2)-B(1)*cmplx(D(3),D(4),kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  
  if (DEBUG) then
     !!!!!!!!!!!!!!!!!!!!
     ! check A)
     ! wrong JOB
     ! set variables
     job = 'K'
     D(1) = 1d0
     D(2) = 0d0
     D(3) = 1d0
     D(4) = 0d0
     a = 1d0
     c = 1d0
     call d_rot2_vec2gen(a,c,B(1),B(2),nrm,info)
     ! check info
     if (info.NE.0) then
        call u_test_failed(__LINE__)
     end if
     
     call d_rot2_swapdiag(JOB,D,B,INFO)
     ! check info
     if (info.NE.-1) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check B)
     ! wrong D
     ! set variables
     job = 'L'
     a = 1d0
     c = 1d0
     call d_rot2_vec2gen(a,c,B(1),B(2),nrm,info)
     ! check info
     if (info.NE.0) then
        call u_test_failed(__LINE__)
     end if
     
     D(1) = 2d0
     D(2) = 0d0
     D(3) = 1d0
     D(4) = 0d0
     call d_rot2_swapdiag(JOB,D,B,INFO)
     ! check info
     if (info.NE.-2) then
        print*, info, D
        call u_test_failed(__LINE__)
     end if

     D(1) = 1d0
     D(2) = -1d0
     D(3) = 1d0
     D(4) = 0d0
     call d_rot2_swapdiag(JOB,D,B,INFO)
     ! check info
     if (info.NE.-2) then
        call u_test_failed(__LINE__)
     end if

     D(1) = 1d0
     D(2) = -1d0
     D(3) = 1d0
     D(4) = 0d0
     call d_rot2_swapdiag(JOB,D,B,INFO)
     ! check info
     if (info.NE.-2) then
        call u_test_failed(__LINE__)
     end if
     
     D(1) = 1d0
     D(2) = 0d0
     D(3) = 2d-1
     D(4) = 0d0
     call d_rot2_swapdiag(JOB,D,B,INFO)
     ! check info
     if (info.NE.-2) then
        call u_test_failed(__LINE__)
     end if

     D(1) = 1d0
     D(2) = 0d0
     D(3) = 1d0
     D(4) = 2d0
     call d_rot2_swapdiag(JOB,D,B,INFO)
     ! check info
     if (info.NE.-2) then
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check C)
     ! wrong B
     ! set variables
     job = 'L'
     D(1) = 1d0
     D(2) = 0d0
     D(3) = 1d0
     D(4) = 0d0
     a = huge(1d0)
     c = 0d0
     B(1) = a/c
     B(2) = 0d0
     call d_rot2_swapdiag(JOB,D,B,INFO)
     ! check info
     if (info.NE.-3) then
        call u_test_failed(__LINE__)
     end if
     
     a = huge(1d0)
     c = 0d0
     B(1) = 1d0
     B(2) = -a/c
     call d_rot2_swapdiag(JOB,D,B,INFO)
     ! check info
     if (info.NE.-3) then
        call u_test_failed(__LINE__)
     end if

     a = huge(1d0)
     c = 0d0
     B(1) = c/c
     B(2) = -1d0
     call d_rot2_swapdiag(JOB,D,B,INFO)
     ! check info
     if (info.NE.-3) then
        call u_test_failed(__LINE__)
     end if

     a = huge(1d0)
     c = 0d0
     B(1) = 2d0
     B(2) = c/c
     call d_rot2_swapdiag(JOB,D,B,INFO)
     ! check info
     if (info.NE.-3) then
        call u_test_failed(__LINE__)
     end if

     B(1) = 1.2d0
     B(2) = 0d0
     call d_rot2_swapdiag(JOB,D,B,INFO)
     ! check info
     if (info.NE.-3) then
        call u_test_failed(__LINE__)
     end if
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_d_rot2_swapdiag
