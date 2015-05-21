#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! benchmark_z_rot3_vs_rot4_vec3/4gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program benchmark_z_rot3_vs_rot4_vec

  implicit none

  ! parameter
  integer, parameter :: num_trials = 1000000 ! 1 mio trials

  ! compute variables
  integer, parameter :: M = 1000
  integer :: ii, n, ind, jj
  integer, allocatable :: seed(:)
  real(8) :: Q(4*M), W(3*M), nrm
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
  
  do ii=1,M
     ind = 3*ii-2
     call random_number(W(ind))
     call random_number(W(ind+1))
     call random_number(W(ind+2))
  end do

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start) 

  ! loop
  do ii=1,num_trials
     jj = mod(ii,M)+1

     ind = 3*jj-2;
     ! vec3gen W
     call z_rot3_vec3gen(W(ind),W(ind+1),W(ind+2),W(ind),W(ind+1),W(ind+2),nrm)
      
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

  do ii=1,M
     ind = 4*ii-3
     call random_number(Q(ind))
     call random_number(Q(ind+1))
     call random_number(Q(ind+2))
     call random_number(Q(ind+3))
  end do

  ! start timer
  call system_clock(count=c_start) 

  ! loop
  do ii=1,num_trials
     ! vecgen Q
     jj = mod(ii,M)+1
     ind = 4*jj-3;
     call z_rot4_vec4gen(Q(ind),Q(ind+1),Q(ind+2),Q(ind+3),Q(ind),Q(ind+1),Q(ind+2),Q(ind+3),nrm)
      
  end do

  ! stop timer
  call system_clock(count=c_stop)
  
  ! set base time
  rot4_time = dble(c_stop-c_start)/dble(c_rate)
  
  ! subract off base time and average
  rot4_time = rot4_time/dble(num_trials)/10d0

  ! print results
  call u_benchmark_print(rot4_time,0d0)
  
end program benchmark_z_rot3_vs_rot4_vec
