#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_rot3_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_rot3_check. 
! The following tests are run:
!
! 0) contains an INF or NAN
!
! 1) CR = 1; CI = 0; S = 0                
!    and replace 0 by EISCOR_DBL_EPS 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rot3_check

  implicit none

  ! parameter
  real(8), parameter :: tol = EISCOR_DBL_EPS ! accuracy (tolerance)

  ! compute variables
  logical :: FLAG
  real(8) :: CR, CI, S
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
  ! check 0)

    ! one INF
    CR = 1d0; CI = 1d0; S = inf
    call z_rot3_check(CR,CI,S,FLAG)
    if (FLAG) then
      call u_test_failed(__LINE__)
    end if
    
    ! one NAN
    CR = 1d0; CI = nan; S = 0d0
    call z_rot3_check(CR,CI,S,FLAG)
    if (FLAG) then
      call u_test_failed(__LINE__)
    end if
    
  !!!!!!!!!!!!!!!!!!!!
  ! check 1)

    ! identity
    CR = 1d0; CI = 0d0; S = 0d0
    call z_rot3_check(CR,CI,S,FLAG)
    if (.NOT.FLAG) then
      call u_test_failed(__LINE__)
    end if
    
    ! nearly identity
    CR = 1d0; CI = tol; S = tol
    call z_rot3_check(CR,CI,S,FLAG)
    if (.NOT.FLAG) then
      call u_test_failed(__LINE__)
    end if
   
  
  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rot3_check
