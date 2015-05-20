#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! benchmark_z_rot3_vs_rot4_turnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program compares the subroutine z_rot3_turnover + z_rot3_swapdiag with
! z_rot4_turnover. 
!
! The following benchmarks are run:
!
! 1) Compute 100,000 rot3 swapdiags and turnovers with a 100,000 rot4 turnovers.
!    Choose two random rotations and a 3x3 random diagonal matrix. move a bulge
!    around all the time to simulate bulge chasing.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program benchmark_z_rot3_vs_rot4_turn

  implicit none

  ! parameter
  integer, parameter :: num_trials = 100000 ! 100,000 trials

  ! compute variables
  integer :: ii, n
  integer, allocatable :: seed(:)
  real(8) :: G1(4), G2(4), B(3), nrm
  real(8) :: D(6)
  real(8) :: rot3_time, rot4_time
  
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
  
  call random_number(G1(1))
  call random_number(G1(2))
  call random_number(G1(3))
  call z_rot3_vec3gen(G1(1),G1(2),G1(3),G1(1),G1(2),G1(3),nrm)

  call random_number(G2(1))
  call random_number(G2(2))
  call random_number(G2(3))
  call z_rot3_vec3gen(G2(1),G2(2),G2(3),G2(1),G2(2),G2(3),nrm)

  call random_number(B(1))
  call random_number(B(2))
  call random_number(B(3))
  call z_rot3_vec3gen(B(1),B(2),B(3),B(1),B(2),B(3),nrm)

  call random_number(D(1))
  call random_number(D(2))
  call d_rot2_vec2gen(D(1),D(2),D(1),D(2),nrm)

  call random_number(D(3))
  call random_number(D(4))
  call d_rot2_vec2gen(D(3),D(4),D(3),D(4),nrm)

  call random_number(D(5))
  call random_number(D(6))
  call d_rot2_vec2gen(D(5),D(6),D(5),D(6),nrm)

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start) 

  ! loop
  do ii=1,num_trials
     ! through diag
     call z_rot3_swapdiag(.FALSE.,D(1:4),B)
     
     ! through Q
     call z_rot3_turnover(G1(1:3),G2(1:3),B)     
      
  end do

  ! stop timer
  call system_clock(count=c_stop)
  
  ! set base time
  rot3_time = dble(c_stop-c_start)/dble(c_rate)
  
  ! subract off base time and average
  rot3_time = rot3_time/dble(num_trials)/10d0
    

  ! print banner
  call u_test_banner(__FILE__) 

  call random_number(G1(1))
  call random_number(G1(2))
  call random_number(G1(3))
  call random_number(G1(4))
  call z_rot4_vec4gen(G1(1),G1(2),G1(3),G1(4),G1(1),G1(2),G1(3),G1(4),nrm)

  call random_number(G2(1))
  call random_number(G2(2))
  call random_number(G2(3))
  call random_number(G2(4))
  call z_rot4_vec4gen(G2(1),G2(2),G2(3),G2(4),G2(1),G2(2),G2(3),G2(4),nrm)

  call random_number(B(1))
  call random_number(B(2))
  call random_number(B(3))
  call z_rot3_vec3gen(B(1),B(2),B(3),B(1),B(2),B(3),nrm)

  ! start timer
  call system_clock(count=c_start) 

  ! loop
  do ii=1,num_trials
     ! through Q
     call z_rot4_turnover(G1(1:4),G2(1:4),B)     
      
  end do

  ! stop timer
  call system_clock(count=c_stop)
  
  ! set base time
  rot4_time = dble(c_stop-c_start)/dble(c_rate)
  
  ! subract off base time and average
  rot4_time = rot4_time/dble(num_trials)/10d0

  ! print results
  call u_benchmark_print(rot4_time,0d0)
  
end program benchmark_z_rot3_vs_rot4_turn
