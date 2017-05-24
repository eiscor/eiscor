#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_1Darray_random_normal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_d_1Darray_random_normal.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_1Darray_random_normal

  implicit none

  ! compute variables
  real(8) :: A(3)
  integer :: info
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)

  call u_fixedseed_initialize(info)

  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  call d_1Darray_random_normal(3,A)

  if (abs(A(1)-0.47761357716530684d0)>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if

  if (abs(A(2)-7.3805546426916641d-2)>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if

  if (abs(A(3)-0.75081696949437793d0)>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_d_1Darray_random_normal
