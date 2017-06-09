#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_2Darray_random_normal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_z_2Darray_random_normal.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_2Darray_random_normal

  implicit none

  ! compute variables
  complex(8) :: A(3,2)
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

  call z_2Darray_random_normal(3,2,A)
  
  if (abs(A(1,1)-cmplx(-0.66808819570369726d0, 0.50753608835290620d0, kind=8))>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if
  
  if (abs(A(1,2)-cmplx(-0.39181453064203980d0, 0.90969551610559185d0, kind=8))>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if
  
  if (abs(A(2,1)-cmplx(-2.5355156920362356d0, -0.94358596649153170d0, kind=8))>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if

  if (abs(A(2,2)-cmplx(0.60792575305475460d0, -1.6421048016782014d0, kind=8))>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if
  
  if (abs(A(3,1)-cmplx(0.27405052576540079d0, 0.92278425003955999d0, kind=8))>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if
  
  if (abs(A(3,2)-cmplx(0.77196141468304724d0, 0.10360858696343152d0, kind=8))>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
end program test_z_2Darray_random_normal
