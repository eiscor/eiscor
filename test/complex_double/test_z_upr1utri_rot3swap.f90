#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1utri_rot3swap
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1utri_rot3swap. 
! The following tests are run (once with 'L2R' and once with 'R2L'):
!
! 1)  D = I, C = B = all [0, -1; 1, 0] and G = [1/sqrt(2) 0 1/sqrt(2)]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1utri_rot3swap

  implicit none

  ! parameter
  real(8), parameter :: tol = 2d0*EISCOR_DBL_EPS ! accuracy (tolerance)

  ! compute variables
  real(8) :: Dold(4), Cold(6), Bold(6), Gold(3)
  real(8) :: D(4), C(6), B(6), G(3)

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
    D(1) = 1d0
    D(2) = 0d0
    D(3) = 1d0
    D(4) = 0d0

    C(1) = 0d0
    C(2) = 0d0
    C(3) = 1d0
    C(4) = 0d0
    C(5) = 0d0
    C(6) = 1d0
    
    B = C
    
    G(1) = 1d0/sqrt(2d0)
    G(2) = 0d0
    G(3) = G(1)
    
    Dold = D
    Cold = C
    Bold = B
    Gold = G
    
    ! call 
    call z_upr1utri_rot3swap(.TRUE.,D,C,B,G)
    
    ! call 
    call z_upr1utri_rot3swap(.FALSE.,D,C,B,G)

    ! check results
    if (maxval(abs(Dold-D)) > tol) then
      call u_test_failed(__LINE__)
    end if
    if (maxval(abs(Cold-C)) > tol) then
      call u_test_failed(__LINE__)
    end if  
    if (maxval(abs(Bold-B)) > tol) then
      call u_test_failed(__LINE__)
    end if  
    if (maxval(abs(Gold-G)) > tol) then
      call u_test_failed(__LINE__)
    end if    

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_upr1utri_rot3swap
