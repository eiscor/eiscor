#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_scalar_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_scalar_check. 
! The following tests are run:
!
! 1) test  0d0 + i0d0     
!
! 2) test 1.5d1 + i0d0            
!
! 3) test -2.42d24 + i1d0             
!
! 4) test  INF + i0d0
!
! 5) test  0d0 + iINF
!
! 6) test -INF + i0d0
!
! 7) test  0d0 - iINF
!
! 8) test  NAN + i0d0
!
! 9) test  0d0 + iNAN
!
! 10) test  huge(1d0) + 0d0
!
! 11) test  -huge(1d0) + 0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_scalar_check

  implicit none

  ! compute variables
  complex(8) :: num
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
  num = cmplx(0d0,0d0,kind=8)
  call z_scalar_check(num,flag)
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  num = cmplx(1.5d1,0d0,kind=8)
  call z_scalar_check(num,flag)
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  num = cmplx(-2.42d24,1d0,kind=8)
  call z_scalar_check(num,flag)
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  num = cmplx(inf,0d0,kind=8)
  call z_scalar_check(num,flag)
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 5)

  ! set variables
  num = cmplx(0d0,inf,kind=8)
  call z_scalar_check(num,flag)
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 6)

  ! set variables
  num = cmplx(-inf,0d0,kind=8)
  call z_scalar_check(num,flag)
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 7)

  ! set variables
  num = cmplx(0d0,-inf,kind=8)
  call z_scalar_check(num,flag)
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 8)

  ! set variables
  num = cmplx(nan,0d0,kind=8)
  call z_scalar_check(num,flag)
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 9)

  ! set variables
  num = cmplx(0d0,nan,kind=8)
  call z_scalar_check(num,flag)
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 10)

  ! set variables
  num = cmplx(huge(1d0),0d0,kind=8)
  call z_scalar_check(num,flag)
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 11)

  ! set variables
  num = cmplx(-huge(1d0),0d0,kind=8)
  call z_scalar_check(num,flag)
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_scalar_check
