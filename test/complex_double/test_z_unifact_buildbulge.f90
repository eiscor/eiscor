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
! 1) K = 1, shift = 0, shift = 0.8-1/sqrt(2) + i(0.6-1/sqrt(2))
! 2) K = 2, shift = 0, shift = 0.7 + i0.7
!
! in DEBUG mode additionally the check of the input data is checked
! A) N = 1
! B) k = 0, K = 4
! C) SHIFT contains nan and inf
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_unifact_buildbulge

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*epsilon(1d0) ! accuracy (tolerance)
  
  ! compute variables
  integer :: N = 4, K, info
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
  
  K = 1

  SHIFT = cmplx(0d0,0d0,kind=8)

  call z_unifact_buildbulge(N,K,Q,D,SHIFT,B,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if

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

  call z_unifact_buildbulge(N,K,Q,D,SHIFT,B,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if

  if (abs(B(1)-1d0/sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(2)-1d0/sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(3))>tol) then
    call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  K = 2

  SHIFT = cmplx(0d0,0d0,kind=8)

  call z_unifact_buildbulge(N,K,Q,D,SHIFT,B,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if

  if (abs(B(1)+0.48d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(2)+0.14d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(3)-sqrt(3d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  SHIFT = cmplx(+0.7d0,+0.7d0,kind=8)

  call z_unifact_buildbulge(N,K,Q,D,SHIFT,B,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if

  if (abs(B(1)-0.5d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(3)-sqrt(0.75d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (DEBUG) then
     !!!!!!!!!!!!!!!!!!!!
     ! check A)
     ! N = 1
     call z_unifact_buildbulge(1,K,Q,D,SHIFT,B,INFO)
     ! check info
     if (info.NE.-1) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check B)
     K = 0
     call z_unifact_buildbulge(N,K,Q,D,SHIFT,B,INFO)
     ! check info
     if (info.NE.-2) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     K = 4
     call z_unifact_buildbulge(N,K,Q,D,SHIFT,B,INFO)
     ! check info
     if (info.NE.-2) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check C)
     K = 1
     nul = 0d0

     SHIFT = cmplx(1d0/nul,0d0,kind=8)
     call z_unifact_buildbulge(N,K,Q,D,SHIFT,B,INFO)
     ! check info
     if (info.NE.-5) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     SHIFT = cmplx(0d0,1d0/nul,kind=8)
     call z_unifact_buildbulge(N,K,Q,D,SHIFT,B,INFO)
     ! check info
     if (info.NE.-5) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     SHIFT = cmplx(nul/nul,0d0,kind=8)
     call z_unifact_buildbulge(N,K,Q,D,SHIFT,B,INFO)
     ! check info
     if (info.NE.-5) then
        print*, info
        call u_test_failed(__LINE__)
     end if


     SHIFT = cmplx(0d0,nul/nul,kind=8)
     call z_unifact_buildbulge(N,K,Q,D,SHIFT,B,INFO)
     ! check info
     if (info.NE.-5) then
        print*, info
        call u_test_failed(__LINE__)
     end if

  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_unifact_buildbulge
