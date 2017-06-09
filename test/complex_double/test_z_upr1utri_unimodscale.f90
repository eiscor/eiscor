#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1utri_unimodscale
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1utri_unimodscale. 
! The following tests are run (both for row and column scaling):
!
! 1)  D = 0-1i, C = B = [1+0i,1]/sqrt(2) and SCL = 0+1i
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1utri_unimodscale

  implicit none

  ! compute variables
  real(8) :: rt2, D(2), C(3), B(3)
  complex(8) :: SCL

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
    D(1) = 0d0
    D(2) = -1d0

    rt2 = 1d0/sqrt(2d0)
    B(1) = rt2
    B(2) = 0d0
    B(3) = rt2
    
    C = B
    
    SCL = cmplx(0d0,1d0,kind=8)
    
    ! call 
    call z_upr1utri_unimodscale(.FALSE.,D,C,B,SCL)
    
    ! check D
    if ( (D(1).NE.1d0).OR.(D(2).NE.0d0) ) then
      call u_test_failed(__LINE__)
    end if

    ! check C
    if ( (C(1).NE.0d0).OR.(C(2).NE.-rt2).OR.(C(3).NE.rt2) ) then
      call u_test_failed(__LINE__)
    end if

    ! check B
    if ( (B(1).NE.0d0).OR.(B(2).NE.rt2).OR.(B(3).NE.rt2) ) then
      call u_test_failed(__LINE__)
    end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_upr1utri_unimodscale
