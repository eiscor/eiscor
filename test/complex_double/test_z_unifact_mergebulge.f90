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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_unifact_mergebulge

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*EISCOR_DBL_EPS ! accuracy (tolerance)
  
  ! compute variables
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
 
  B(1) = 1/sqrt(3d0)
  B(2) = -1/sqrt(3d0)
  B(3) = -1/sqrt(3d0)
! this is not finished yet
  call z_unifact_mergebulge(.TRUE.,Q,D,B)

  if (abs(Q(1)+1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)+1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(1)-0.6d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(2)-0.6d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(3)-0.6d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(4)-0.6d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(1)-8d-1)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(2)-6d-1)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(3))>tol) then
    call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  Q = Qs
  D = Ds
  call z_unifact_mergebulge(.FALSE.,Q(7:9),D(5:8),B)

  if (abs(Q(7)-0.32957819226739016d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(8)+0.25118898775644621d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(9)-0.91010016350490297d0)>tol) then
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

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_unifact_mergebulge
