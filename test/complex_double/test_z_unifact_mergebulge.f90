#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_unifact_mergebulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_z_unifact_mergebulge.
! The following tests are run:
!
! Matrix:
! Q1=[0.8, 0.6, 0.0], 
! Q2=[-1/sqrt(3),1/sqrt(3), 1/sqrt(3)]
! Q3=[0.3,0.4,-sqrt(3)/2], 
! D = diag([1, 0.8+i0,6, 1/sqrt(2)+i/sqrt(2), 1])
!
! 1) merge at top  
! 2) merge at bottom
!    Note: some results are not know exact
!
! Additionally in DEBUG mode the following tests are performed
! A) incorrect jobname
! B) N too small
! C) STR out of range
! D) STP out of range
! E) B contrains inf or nan
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_unifact_mergebulge

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*epsilon(1d0) ! accuracy (tolerance)
  
  ! compute variables
  integer :: N = 4, info, str, stp
  real(8) :: Q(9), Qs(9), D(8), Ds(8), B(3), C(3), nul
  
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
  Q(4) = -1d0/sqrt(3d0)
  Q(5) = 1d0/sqrt(3d0)
  Q(6) = 1d0/sqrt(3d0)
  Q(7) = 3d-1
  Q(8) = 4d-1
  Q(9) = -sqrt(3d0)/2d0
  Qs = Q
  
  D(1) = 1d0
  D(2) = 0d0
  D(3) = 0.8d0
  D(4) = 0.6d0
  D(5) = 1d0/sqrt(2d0)
  D(6) = 1d0/sqrt(2d0)
  D(7) = 1d0
  D(8) = 0d0
  Ds = D

  nul = 0d0
 
  STR = 1
  STP = 3

  B(1) = 1/sqrt(3d0)
  B(2) = -1/sqrt(3d0)
  B(3) = -1/sqrt(3d0)

  call z_unifact_mergebulge('T',N,STR,STP,Q,D,B,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  
  if (abs(Q(1)+1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(4)-0.2d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(5)+1.4d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(6)-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(7)+0.48d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(8)+0.14d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(9)+sqrt(3d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (abs(D(1)+0.8d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(2)+0.6d0)>tol) then
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
  if (abs(D(7)+0.8d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(8)-0.6d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  Q = Qs
  D = Ds
  call z_unifact_mergebulge('B',N,STR,STP,Q,D,B,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  
  if (abs(Q(1)-0.8d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-0.6d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(4)+1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(5)-1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(6)-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(7)-0.32957819226739016d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(8)+0.25118898775644621d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(9)-0.91010016350490297d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  
  if (abs(D(1)-1.0d0)>tol) then
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
  if (abs(D(5)+0.96726920535125560d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(6)-0.25375240763222456d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(7)+0.50453256615763697d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(8)+0.86339266251595492d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (DEBUG) then
     
     !!!!!!!!!!!!!!!!!!!!
     ! check A)
     call z_unifact_mergebulge('M',N,STR,STP,Q,D,B,INFO)
     
     ! check info
     if (INFO.NE.-1) then
        call u_test_failed(__LINE__)
     end if
     
     !!!!!!!!!!!!!!!!!!!!
     ! check B)
     call z_unifact_mergebulge('T',1,STR,STP,Q,D,B,INFO)
     
     ! check info
     if (INFO.NE.-2) then
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check C)
     call z_unifact_mergebulge('T',N,0,STP,Q,D,B,INFO)
     
     ! check info
     if (INFO.NE.-3) then
        call u_test_failed(__LINE__)
     end if

     call z_unifact_mergebulge('T',N,4,STP,Q,D,B,INFO)
     
     ! check info
     if (INFO.NE.-3) then
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check D)
     call z_unifact_mergebulge('T',N,STP,0,Q,D,B,INFO)
     
     ! check info
     if (INFO.NE.-4) then
        call u_test_failed(__LINE__)
     end if

     call z_unifact_mergebulge('T',N,STR,4,Q,D,B,INFO)
     
     ! check info
     if (INFO.NE.-4) then
        call u_test_failed(__LINE__)
     end if
     
     !!!!!!!!!!!!!!!!!!!!
     ! check E)
     B(1) = 1d0/nul     
     call z_unifact_mergebulge('T',N,STR,STP,Q,D,B,INFO)

     ! check info
     if (INFO.NE.-7) then
        call u_test_failed(__LINE__)
     end if

     B(1) = -1d0/nul    
     call z_unifact_mergebulge('T',N,STR,STP,Q,D,B,INFO)

     ! check info
     if (INFO.NE.-7) then
        call u_test_failed(__LINE__)
     end if

     B(1) = nul/nul
     
     call z_unifact_mergebulge('T',N,STR,STP,Q,D,B,INFO)
     ! check info
     if (INFO.NE.-7) then
        call u_test_failed(__LINE__)
     end if

     
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_unifact_mergebulge
