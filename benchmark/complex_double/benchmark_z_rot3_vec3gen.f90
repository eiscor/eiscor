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
  integer, parameter :: num_trials = 10**9

  ! compute variables
  integer :: ii, n
  integer, allocatable :: seed(:)
  real(8) :: AR, AI, B, CR, CI, S, NRM
  real(8) :: base_time, total_time, ERROR
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! print banner
  call u_test_banner(__FILE__) 
  
  ! get size of seed        
  call random_seed(size = n)
  
  ! allocate memory for seed
  allocate(seed(n))
  
  ! set seed        
  seed = 1
  
  ! set the generator
  call random_seed(put = seed)
  
  ! free memory        
  deallocate(seed)
  
  ! base line for random number generate and error calculation
    ! set ERROR
    ERROR = 0d0
    
    ! set generators
    CR = 1d0; CI = 0d0; S = 0d0; NRM = 1d0
    
    ! start timer
    call system_clock(count_rate=c_rate)
    call system_clock(count=c_start) 
    
    ! loop
    do ii=1,num_trials
      
      ! set AR, AI, B
      call random_number(AR)
      call random_number(AI)
      call random_number(B)
      
      ! compute worse ERROR
      ERROR = max(ERROR,abs(AR-NRM*CR))
      ERROR = max(ERROR,abs(AI-NRM*CI))
      ERROR = max(ERROR,abs(B-NRM*S))
      
    end do  

    ! stop timer
    call system_clock(count=c_stop)
    
    ! set base time
    base_time = dble(c_stop-c_start)/dble(c_rate)
    
  ! total time for random number generate and error calculation
    ! reset ERROR
    ERROR = 0d0
    
    ! start timer
    call system_clock(count_rate=c_rate)
    call system_clock(count=c_start) 
    
    ! loop
    do ii=1,num_trials
      
      ! set AR, AI, B
      call random_number(AR)
      call random_number(AI)
      call random_number(B)
      
      ! compute core transformation
      call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
      
      ! compute worse ERROR
      ERROR = max(ERROR,abs(AR-NRM*CR))
      ERROR = max(ERROR,abs(AI-NRM*CI))
      ERROR = max(ERROR,abs(B-NRM*S))
      
    end do  

    ! stop timer
    call system_clock(count=c_stop)
    
    ! set total time
    total_time = dble(c_stop-c_start)/dble(c_rate)
    
    ! subract off base time and average
    total_time = (total_time - base_time)/dble(num_trials)
    
    ! print results
    call u_benchmark_print(total_time,ERROR)

end program benchmark_z_rot3_vec3gen
