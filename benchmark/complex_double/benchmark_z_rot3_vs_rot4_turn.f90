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
! 1) Compute 1,000,000 rot3 swapdiags and turnovers with a 1,000,000 rot4
!    turnovers.
!    Choose 100 random rotations and a 3x3 random diagonal matrix. Move a bulge
!    around all the time to simulate bulge chasing.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program benchmark_z_rot3_vs_rot4_turn

  implicit none

  ! parameter
  integer, parameter :: num_trials = 1000000 ! 1 mio trials

  ! compute variables
  integer :: ii, n, ind, jj
  integer, allocatable :: seed(:)
  real(8) :: Q(4*99), B(4), nrm
  real(8) :: D(2*100)
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
  
  do ii=1,99
     ind = 3*ii-2
     call random_number(Q(ind))
     call random_number(Q(ind+1))
     call random_number(Q(ind+2))
     call z_rot3_vec3gen(Q(ind),Q(ind+1),Q(ind+2),Q(ind),Q(ind+1),Q(ind+2),nrm)
  end do

  call random_number(B(1))
  call random_number(B(2))
  call random_number(B(3))
  call z_rot3_vec3gen(B(1),B(2),B(3),B(1),B(2),B(3),nrm)

  do ii=1,100
     ind = 2*ii-1
     call random_number(D(ind))
     call random_number(D(ind+1))
     call d_rot2_vec2gen(D(ind),D(ind+1),D(ind),D(ind+1),nrm)
  end do

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start) 

  ! loop
  do ii=1,num_trials
     jj = mod(ii,98)+1
     ind = 2*jj-1;
     ! through diag
     call z_rot3_swapdiag(.FALSE.,D(ind:ind+3),B)
     
     ind = 3*jj-2;
     ! through Q
     call z_rot3_turnover(Q(ind:ind+2),Q(ind+3:ind+5),B)     
      
  end do

  ! stop timer
  call system_clock(count=c_stop)
  
  ! set base time
  rot3_time = dble(c_stop-c_start)/dble(c_rate)
  
  ! subract off base time and average
  rot3_time = rot3_time/dble(num_trials)/10d0

  ! print results
  call u_benchmark_print(rot3_time,0d0)

  print*, ""
  print*, ""


  ! print banner
  call u_test_banner(__FILE__) 

  do ii=1,99
     ind = 4*ii-3
     call random_number(Q(ind))
     call random_number(Q(ind+1))
     call random_number(Q(ind+2))
     call random_number(Q(ind+3))
     call z_rot4_vec4gen(Q(ind),Q(ind+1),Q(ind+2),Q(ind+3),Q(ind),Q(ind+1),Q(ind+2),Q(ind+3),nrm)
  end do

  call random_number(B(1))
  call random_number(B(2))
  call random_number(B(3))
  call z_rot3_vec3gen(B(1),B(2),B(3),B(1),B(2),B(3),nrm)

  ! start timer
  call system_clock(count=c_start) 

  ! loop
  do ii=1,num_trials
     ! through Q
     jj = mod(ii,98)+1
     ind = 4*jj-3;
     call z_rot4_turnover(Q(ind:ind+3),Q(ind+4:ind+7),B)     
      
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
