#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! benchmark_d_rot2_turnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program benchmarks the subroutine d_rot2_turnover. 
! The following benchmarks are run:
!
! 1) Compute ten million turnovers from uniformly generated 
!    user input. The average compute time and the worst error are printed.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program benchmark_d_rot2_turnover

  implicit none

  ! parameter
  integer, parameter :: num_trials = 10**7 ! ten million trials

  ! compute variables
  integer :: ii, n
  integer, allocatable :: seed(:)
  real(8) :: G1(2), G2(2), G3(2), nrm
  real(8) :: base_time, total_time, ERROR
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! print banner
  call u_test_banner(__FILE__) 
  
  ! get size of seed        
  call random_seed(size = n)
  
  ! allocate memory for seed
  allocate(seed(n))
  
  ! check allocation
  if (.NOT.allocated(seed)) then
    call u_test_failed(__LINE__)
  end if
  
  ! set seed        
  seed = 1
  
  ! set the generator
  call random_seed(put = seed)
  
  ! free memory        
  deallocate(seed)
  
  ! base line for random number generate and error calculation
    ! set ERROR
    ERROR = 0d0
    
    ! start timer
    call system_clock(count_rate=c_rate)
    call system_clock(count=c_start) 
    
    ! loop
    do ii=1,num_trials
      
      ! set G1
      call random_number(G1(1))
      call random_number(G1(2))
      call d_rot2_vec2gen(G1(1),G1(2),G1(1),G1(2),nrm)
      
      ! set G2 and G3
      G2 = G1; G3 = G1
      
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
      
      ! set G1
      call random_number(G1(1))
      call random_number(G1(2))
      call d_rot2_vec2gen(G1(1),G1(2),G1(1),G1(2),nrm)
      
      ! set G2 and G3
      G2 = G1; G3 = G1
      
      ! compute 10 turnovers
      call d_rot2_turnover(G1,G2,G3)
      call d_rot2_turnover(G1,G2,G3)
      call d_rot2_turnover(G1,G2,G3)
      call d_rot2_turnover(G1,G2,G3)
      call d_rot2_turnover(G1,G2,G3)
      call d_rot2_turnover(G1,G2,G3)
      call d_rot2_turnover(G1,G2,G3)
      call d_rot2_turnover(G1,G2,G3)
      call d_rot2_turnover(G1,G2,G3)
      call d_rot2_turnover(G1,G2,G3)
      
    end do  

    ! stop timer
    call system_clock(count=c_stop)
    
    ! set total time
    total_time = dble(c_stop-c_start)/dble(c_rate)
    
    ! subract off base time and average
    total_time = (total_time - base_time)/dble(num_trials)/10d0
    
    ! print results
    call u_benchmark_print(total_time,ERROR)

end program benchmark_d_rot2_turnover
