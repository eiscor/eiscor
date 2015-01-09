#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_twistedqz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_twistedqz. 
! The following tests are run:
!
! 1) check roots of unity with upperhess QR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_twistedqz

  implicit none
  
  ! compute variables
  integer, parameter :: N = 3
  integer :: ii, INFO, ITS(N-1)
  logical :: P(N-2), HESS
  real(8) :: Q(3*(N-1)), D(2,2*(N+1)), R(4,3*N)
  complex(8) :: V(N,N), W(N,N)
  
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
    
    ! set P
    P = .FALSE.
    
    ! set valid Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii) = 1d0
    end do     
  
    ! set valid D
    D = 0d0
    do ii=1,(N+1)
      D(:,2*ii-1) = 1d0
    end do

    ! set valid R
    R = 0d0
    do ii=1,N
      R(1,3*ii) = -1d0
      R(2,3*ii) = 1d0
      R(3,3*ii) = -1d0
      R(4,3*ii) = 1d0
    end do
    
    ! call twisted QZ
    call z_upr1fact_twistedqz('QR','I',N,P,HESS,Q,D,R,V,W,ITS,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if
  ! end check 1)
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_upr1fact_twistedqz
