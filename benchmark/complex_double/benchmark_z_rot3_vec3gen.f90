#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! benchmark_z_rot3_vec3gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program benchmarks the subroutine z_rot3_vec3gen. 
! The following benchmarks are run:
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program benchmark_z_rot3_vec3gen

  implicit none

  ! parameter
  integer, parameter :: time_trials = 10**1
  integer, parameter :: error_trials = 10**6

  ! compute variables
  integer :: ii
  real(8) :: AR, AI, B, CR, CI, S, NRM
  real(8) :: XR, XI, Y, ERROR
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  !!!!!!!!!!!!!!!!!!!
  ! Configuration 0
    ! set AR, AI, B
    AR = 0d0; AI = 0d0; B = 0d0
    
    ! start timer
    call system_clock(count_rate=c_rate)
    call system_clock(count=c_start) 
    
    ! loop
    do ii=1,time_trials
      call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
    end do  

    ! stop timer
    call system_clock(count=c_stop)
    
    ! set XR, XI, Y
    XR = AR; XI = AI; Y = B
    
    ! set ERROR
    ERROR = 0d0
    
    ! loop
    do ii=1,error_trials
      call z_rot3_vec3gen(XR,XI,Y,CR,CI,S,NRM)
      XR = NRM*CR
      XI = NRM*CI
      Y = NRM*S
      ERROR = max(ERROR,abs(XR-AR))
      ERROR = max(ERROR,abs(XI-AI))
      ERROR = max(ERROR,abs(Y-B))
    end do  
    
    ! print time
    print*, ""
    print*, "Configuration 0"
    print*, "  Average time:", (dble(c_stop-c_start)/dble(c_rate))/dble(time_trials)
    print*, "  Worst error:", ERROR
    print*, ""
  
  !!!!!!!!!!!!!!!!!!!
  ! Configuration 1
    ! set AR, AI, B
    AR = 1d0; AI = 100d0*sqrt(EISCOR_DBL_EPS); B = EISCOR_DBL_EPS
    
    ! start timer
    call system_clock(count_rate=c_rate)
    call system_clock(count=c_start) 
    
    ! loop
    do ii=1,time_trials
      call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
    end do  

    ! stop timer
    call system_clock(count=c_stop)
    
    ! set XR, XI, Y
    XR = AR; XI = AI; Y = B
    
    ! set ERROR
    ERROR = 0d0
    
    ! loop
    do ii=1,error_trials
      call z_rot3_vec3gen(XR,XI,Y,CR,CI,S,NRM)
      XR = NRM*CR
      XI = NRM*CI
      Y = NRM*S
      ERROR = max(ERROR,abs(XR-AR))
      ERROR = max(ERROR,abs(XI-AI))
      ERROR = max(ERROR,abs(Y-B))
    end do  
    
    ! print time
    print*, ""
    print*, "Configuration 1"
    print*, "  Average time:", (dble(c_stop-c_start)/dble(c_rate))/dble(time_trials)
    print*, "  Worst error:", ERROR
    print*, ""
    
  !!!!!!!!!!!!!!!!!!!
  ! Configuration 2
    ! set AR, AI, B
    AR = sqrt(EISCOR_DBL_EPS); AI = 1d0; B = EISCOR_DBL_EPS
    
    ! start timer
    call system_clock(count_rate=c_rate)
    call system_clock(count=c_start) 
    
    ! loop
    do ii=1,time_trials
      call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
    end do  

    ! stop timer
    call system_clock(count=c_stop)
    
    ! set XR, XI, Y
    XR = AR; XI = AI; Y = B
    
    ! set ERROR
    ERROR = 0d0
    
    ! loop
    do ii=1,error_trials
      call z_rot3_vec3gen(XR,XI,Y,CR,CI,S,NRM)
      XR = NRM*CR
      XI = NRM*CI
      Y = NRM*S
      ERROR = max(ERROR,abs(XR-AR))
      ERROR = max(ERROR,abs(XI-AI))
      ERROR = max(ERROR,abs(Y-B))
    end do  
    
    ! print time
    print*, ""
    print*, "Configuration 2"
    print*, "  Average time:", (dble(c_stop-c_start)/dble(c_rate))/dble(time_trials)
    print*, "  Worst error:", ERROR
    print*, ""
    
  !!!!!!!!!!!!!!!!!!!
  ! Configuration 3
    ! set AR, AI, B
    AR = sqrt(EISCOR_DBL_EPS); AI = EISCOR_DBL_EPS; B = 1d0
    
    ! start timer
    call system_clock(count_rate=c_rate)
    call system_clock(count=c_start) 
    
    ! loop
    do ii=1,time_trials
      call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
    end do  

    ! stop timer
    call system_clock(count=c_stop)
    
    ! set XR, XI, Y
    XR = AR; XI = AI; Y = B
    
    ! set ERROR
    ERROR = 0d0
    
    ! loop
    do ii=1,error_trials
      call z_rot3_vec3gen(XR,XI,Y,CR,CI,S,NRM)
      XR = NRM*CR
      XI = NRM*CI
      Y = NRM*S
      ERROR = max(ERROR,abs(XR-AR))
      ERROR = max(ERROR,abs(XI-AI))
      ERROR = max(ERROR,abs(Y-B))
    end do  
    
    ! print time
    print*, ""
    print*, "Configuration 3"
    print*, "  Average time:", (dble(c_stop-c_start)/dble(c_rate))/dble(time_trials)
    print*, "  Worst error:", ERROR
    print*, ""

end program benchmark_z_rot3_vec3gen
