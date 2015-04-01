#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_1Darray_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_d_1Darray_check (1d array check). 
! The following tests are run:
!
! 1) test 0d0, 1.5d1, -2.42d24             
!
! 2) test with one INF
!
! 3) test with one -INF
!
! 4) test with one NAN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_1Darray_check

  implicit none

  ! compute variables
  real(8) :: num,nul
  real(8) :: A(3)
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
  A(1) = 0d0
  A(2) = 1.5d1
  A(3) = -2.42d24
  call d_1Darray_check(3,A,flag)
  ! check flag
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  A(2) = 1.5d1
  A(3) = -2.42d24
  num = +huge(1d0)
  A(1) = num+huge(1d0)
  call d_1Darray_check(3,A,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  A(1) = 0d0
  A(3) = -2.42d24
  num = -huge(1d0)
  A(2) = num-huge(1d0)
  call d_1Darray_check(3,A,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  A(1) = 0d0
  A(2) = 1.5d1
  nul = 0d0
  A(3) = nul/nul
  call d_1Darray_check(3,A,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_d_1Darray_check
