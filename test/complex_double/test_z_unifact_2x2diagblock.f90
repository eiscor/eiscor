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
!    The blocks (1:2,1:2), (2:3,2:3), and (3:4,3:4) are checked.
!
! 2) Q1=[0.8, 0.6, 0.0], Q2=[0.3,0.4,-sqrt(3)/2], 
!    Q3=[-1/sqrt(3),1/sqrt(3), 1/sqrt(3)]
!    D = diag([1, 0.8+i0,6, 1/sqrt(2)+i/sqrt(2), 1])
!    The blocks (2:3,2:3) are checked.
!
!
! in DEBUG mode additionally the check of the input data is checked
! A) N = 1
! B) K = 0
! C) K = N
! D) Q contains INF
! E) D contains INF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_unifact_2x2diagblock

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*epsilon(1d0) ! accuracy (tolerance)
  
  ! compute variables
  integer :: info
  real(8) :: Q(9),Qs(9) ! N = 4
  real(8) :: D(8), nul = 0d0
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

  call z_unifact_2x2diagblock(4,1,Q,D,H,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  
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

  call z_unifact_2x2diagblock(4,2,Q,D,H,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  
  ! check H
  if (abs(dble(H(1,1))+0.48d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(1,1))+0.14d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(dble(H(2,1))-sqrt(3d0)/2d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(2,1)))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(dble(H(1,2))-0.1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(1,2))+0.7d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(dble(H(2,2))+0.1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(2,2))+0.7d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  call z_unifact_2x2diagblock(4,3,Q,D,H,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  
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

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  Q(1) = 8d-1
  Q(2) = 6d-1
  Q(3) = 0d0
  Q(4) = 3d-1
  Q(5) = 4d-1
  Q(6) = -sqrt(3d0)/2d0
  Q(7) = -1d0/sqrt(3d0)
  Q(8) = 1d0/sqrt(3d0)
  Q(9) = 1d0/sqrt(3d0)
  
  !D = diag([1, 0.8+i0,6, 1/sqrt(2)+i/sqrt(2), 1])
  D(1) = 1d0
  D(2) = 0d0
  D(3) = 0.8d0
  D(4) = 0.6d0
  D(5) = 1d0/sqrt(2d0)
  D(6) = 1d0/sqrt(2d0)
  D(7) = 1d0
  D(8) = 0d0

  call z_unifact_2x2diagblock(4,2,Q,D,H,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! check H
  if (abs(dble(H(1,1))-0.3d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(1,1))-0.4d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(dble(H(2,1))+sqrt(3d0)*0.4d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(2,1))+sqrt(3d0)*0.3d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(dble(H(1,2))+0.8d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(1,2))-0.6d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(dble(H(2,2))+0.6d0/sqrt(6d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(aimag(H(2,2))-0.8d0/sqrt(6d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  if (DEBUG) then
     !!!!!!!!!!!!!!!!!!!!
     ! check A)
     ! N = 1
     call z_unifact_2x2diagblock(1,1,Q,D,H,INFO)
     ! check info
     if (info.NE.-1) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check B)
     call z_unifact_2x2diagblock(4,0,Q,D,H,INFO)
     ! check info
     if (info.NE.-2) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check C)
     ! K = N
     call z_unifact_2x2diagblock(4,4,Q,D,H,INFO)
     ! check info
     if (info.NE.-2) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check D)
     ! Q with INF
     Qs = Q
     Q(1) = 1d0/nul
     call z_unifact_2x2diagblock(4,1,Q,D,H,INFO)
     ! check info
     if (info.NE.-3) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check E)
     ! D with INF
     Q = Qs
     D(1) = 1d0/nul
     call z_unifact_2x2diagblock(4,1,Q,D,H,INFO)
     ! check info
     if (info.NE.-4) then
        print*, info
        call u_test_failed(__LINE__)
     end if

  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_unifact_2x2diagblock
