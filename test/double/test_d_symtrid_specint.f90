#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_symtrid_specint
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_d_symtrid_specint.
! The following tests are run:
!
! 1) T = Toeplitz(.5,0,.5) without Newton correction 
!
! 2) T = Toeplitz(.5,0,.5) with Newton correction 
!
! 3) T = 0 matrix 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_symtrid_specint

  implicit none

  ! compute variables
  integer, parameter :: N = 1000
  real(8) :: D(N), E(N-1)
  real(8) :: A, B, rho
  logical :: flag
  
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
  D = 0d0
  E = 5d-1
  call d_symtrid_specint(.FALSE.,N,D,E,A,B,flag)

  ! check A
  if (A.NE.-1d0) then
     call u_test_failed(__LINE__)
  end if

  ! check B
  if (B.NE.1d0) then
     call u_test_failed(__LINE__)
  end if

  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  rho = cos(EISCOR_DBL_PI/(N+1d0))
  call d_symtrid_specint(.TRUE.,N,D,E,A,B,flag)

  ! check A
  if (A.NE.-rho) then
     call u_test_failed(__LINE__)
  end if

  ! check B
  if (B.NE.rho) then
     call u_test_failed(__LINE__)
  end if

  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  D = 0d0
  E = 0d0
  call d_symtrid_specint(.TRUE.,N,D,E,A,B,flag)

  ! check A
  if (A.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! check B
  if (B.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_d_symtrid_specint
