#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_deflationcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_deflationcheck. 
! The following tests are run:
!
! 1) upper hessenberg
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_deflationcheck

  implicit none
  
  ! compute variables
  integer, parameter :: N = 6
  integer :: ii, ind, INFO, ITCNT, STR, STP, ZERO, ITS(N-1)
  logical :: P(N-1)
  real(8) :: Q(3*(N-1)), D(2,2*(N+1))
  complex(8) :: vec1(N+1), vec2(N+1)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! initialize INFO
  INFO = 0
  
  ! initialize ITS
  ITS = 0
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! check 1)
    ! initialize P
    P = .FALSE.
    
    ! initialize Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii-2) = 1d0/sqrt(2d0)
      Q(3*ii) = 1d0/sqrt(2d0)
    end do
    ind = N/2
    Q(3*ind-2) = 0d0
    Q(3*ind-1) = 1d0
    Q(3*ind) = 0d0
    
    ! initialize D
    D = 0d0
    do ii=1,(N+1)
      D(1,2*ii-1) = 1d0
    end do
    
! print
print*,""
print*,"Q"
do ii=1,(N-1)
print*,Q(3*ii-2),Q(3*ii-1),Q(3*ii)
end do
print*,""

! print
print*,"D"
do ii=1,(N+1)
print*,D(1,2*ii-1),D(1,2*ii)
end do
print*,""
    
    ! set ITCNT
    ITCNT = 10
    
    ! set STR, STP, ZERO
    STR = 1
    STP = N-1
    ZERO = 0
    
    ! call z_upr1fact_deflationcheck
    call z_upr1fact_deflationcheck(N,STR,STP,ZERO,P,Q,D,ITCNT,ITS,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

! print
print*,""
print*,"Q"
do ii=1,(N-1)
print*,Q(3*ii-2),Q(3*ii-1),Q(3*ii)
end do
print*,""

! print
print*,"D"
do ii=1,(N+1)
print*,D(1,2*ii-1),D(1,2*ii)
end do
print*,""

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
end program test_z_upr1fact_deflationcheck
