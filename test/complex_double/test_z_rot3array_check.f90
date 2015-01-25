#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_rot3array_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_rot3array_check.
! The following tests are performed:
!
! 1) N < 1
! 
! 2) A contains an INF
!
! 3) A contains a NAN
!
! 4) A is valid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rot3array_check

  implicit none

  ! parameters
  integer, parameter :: N = 10

  ! compute variables
  logical :: FLAG
  integer :: ii 
  real(8) :: A(3*N), inf, nan
  
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
  
  ! create valid A
  A = 0d0
  do ii=1,N
    A(3*ii-2) = 1d0
  end do
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 0)

    ! N < 1
    call z_rot3array_check(0,A,FLAG)
    if (FLAG) then
      call u_test_failed(__LINE__)
    end if
    
  !!!!!!!!!!!!!!!!!!!!
  ! check 1)

    ! one INF
    A(1) = inf
    call z_rot3array_check(N,A,FLAG)
    if (FLAG) then
      call u_test_failed(__LINE__)
    end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

    ! one NAN
    A(1) = nan
    call z_rot3array_check(N,A,FLAG)
    if (FLAG) then
      call u_test_failed(__LINE__)
    end if
    
  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

    ! valid A
    A(1) = 1d0
    call z_rot3array_check(N,A,FLAG)
    if (.NOT.FLAG) then
      call u_test_failed(__LINE__)
    end if
    
  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rot3array_check
