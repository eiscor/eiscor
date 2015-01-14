#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_mergebulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_mergebulge. 
! The following tests are run:
!
! 1) Q = [0, -1; 1, 0], P = FALSE, D = I, K = 1, JOB = L and 
!    G = [0, 1; -1, 0]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_mergebulge

  implicit none
  
  ! compute variables
  integer, parameter :: N = 3
  integer :: ii
  logical :: P(N-2)
  real(8) :: Q(3*(N-1)),D(2*(N+1)), G(3)
  integer :: INFO
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! check 1)
    ! set INFO
    INFO = 0
  
    ! set Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii) = 1d0
    end do
    
    ! set D
    D = 0d0
    do ii=1,(N+1)
      D(2*ii-1) = 1d0
    end do
    
    ! set G
    G = 0d0
    G(3) = -1d0
    
    ! set P
    P = .FALSE.
    
    ! call merge bulge
    call z_upr1fact_mergebulge('L',N,1,N-1,1,P,Q,D,G,INFO)
  
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_upr1fact_mergebulge
