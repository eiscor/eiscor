#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_scalar_infcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine d_scalar_infcheck (infcheck). 
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_scalar_infcheck

  implicit none

  ! compute variables
  real(8) :: num,nul
  integer :: info
  
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
  num = 0d0
  call d_scalar_infcheck(NUM,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  num = 1.5d1
  call d_scalar_infcheck(NUM,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  num = -2.42d24
  call d_scalar_infcheck(NUM,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  num = +huge(1d0)
  num = num+huge(1d0)
  call d_scalar_infcheck(NUM,INFO)
  ! check info
  if (info.NE.1) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 5)

  ! set variables
  num = -huge(1d0)
  num = num-huge(1d0)
  call d_scalar_infcheck(NUM,INFO)
  ! check info
  if (info.NE.1) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 6)

  ! set variables
  nul = 0d0
  num = nul/nul
  call d_scalar_infcheck(NUM,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_d_scalar_infcheck
