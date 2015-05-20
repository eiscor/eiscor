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
! Q = [-1/sqrt(3), 1/sqrt(3), 1/sqrt(3)]
! D = [0, 1, 0, 1]
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
  real(8) :: Q(3), Qs(3), D(4), Ds(4), B(3), C(3), nul
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__) 

  !!!!!!!!!!!!!!!!!!!!
  ! check 1)
  Q(1) = -1d0/sqrt(3d0)
  Q(2) = 1d0/sqrt(3d0)
  Q(3) = 1d0/sqrt(3d0)
  Qs = Q
  
  D(1) = 0d0
  D(2) = 1d0
  D(3) = 0d0
  D(4) = 1d0
  Ds = D
 
  B(1) = -1/sqrt(3d0)
  B(2) = -1/sqrt(3d0)
  B(3) = -1/sqrt(3d0)

  call z_unifact_mergebulge(.TRUE.,Q,D,B)

  if (abs(Q(1)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (abs(D(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(2)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(4)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (abs(B(1)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(B(3))>tol) then
    call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  Q = Qs
  D = Ds
  B(1) = -1/sqrt(3d0)
  B(2) = -1/sqrt(3d0)
  B(3) = -1/sqrt(3d0)

  call z_unifact_mergebulge(.FALSE.,Q,D,B)

  if (abs(Q(1)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(2))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
     call u_test_failed(__LINE__)
  end if
  
  if (abs(D(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(2)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(D(4)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_unifact_mergebulge
