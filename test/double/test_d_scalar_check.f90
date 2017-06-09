#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_scalar_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine d_scalar_check. 
! The following tests are run:
!
! 1) test 0d0      
!
! 2) test 1.5d1             
!
! 3) test -2.42d24             
!
! 4) test INF
!
! 5) test -INF
!
! 6) test NAN
!
! 7) test huge(1d0)
!
! 8) test -huge(1d0)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_scalar_check

  implicit none

  ! compute variables
  real(8) :: num
  logical :: flag
  real(8) :: inf, nan
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)
  
  ! set inf
  inf = EISCOR_DBL_INF
  inf = 10d0*inf
  
  ! set nan
  nan = 0d0
  nan = 0d0/nan

  !!!!!!!!!!!!!!!!!!!!
  ! check 1)

  ! set variables
  num = 0d0
  call d_scalar_check(NUM,flag)
  ! check flag
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  num = 1.5d1
  call d_scalar_check(NUM,flag)
  ! check flag
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  num = -2.42d24
  call d_scalar_check(NUM,flag)
  ! check flag
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  num = inf
  call d_scalar_check(NUM,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 5)

  ! set variables
  num = -inf
  call d_scalar_check(NUM,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 6)

  ! set variables
  num = nan
  call d_scalar_check(NUM,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 7)

  ! set variables
  num = huge(1d0)
  call d_scalar_check(NUM,flag)
  ! check flag
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 8)

  ! set variables
  num = -huge(1d0)
  call d_scalar_check(NUM,flag)
  ! check flag
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_d_scalar_check
