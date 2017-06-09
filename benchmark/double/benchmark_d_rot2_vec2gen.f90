#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! benchmark_d_rot2_vec2gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program computes two benchmarks for subroutine d_rot2_vec2gen (generating
! rotations).  
!
! 1) Time required for 100 million runs.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program benchmark_d_rot2_vec2gen

  implicit none

  ! parameter
  integer, parameter :: notests = 10**8 ! 100 million
  
  ! compute variables
  integer :: ii, n
  integer, allocatable :: seed(:)
  real(8) :: nrm
  real(8) :: a, b
  real(8) :: c, s
  real(8) :: t1, t2, ERROR

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

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! just random number generation
  do ii=1,notests
    call random_number(a)
    call random_number(b)
     
    ! compute worst ERROR
    ERROR = max(ERROR,abs(a-nrm*c))
    ERROR = max(ERROR,abs(b-nrm*s))
  end do

  ! stop timer
  call system_clock(count=c_stop)
  t1 = dble(c_stop-c_start)/dble(c_rate)

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! set ERROR
  ERROR = 0d0

  ! random number generation + rotation generation
  do ii=1,notests
    call random_number(a)
    call random_number(b)
    call d_rot2_vec2gen(a,b,c,s,nrm)
     
    ! compute worst ERROR
    ERROR = max(ERROR,abs(a-nrm*c))
    ERROR = max(ERROR,abs(b-nrm*s))
    
    ! 9 more times to build up the average
    call d_rot2_vec2gen(a,b,c,s,nrm)
    call d_rot2_vec2gen(a,b,c,s,nrm)
    call d_rot2_vec2gen(a,b,c,s,nrm)
    call d_rot2_vec2gen(a,b,c,s,nrm)
    call d_rot2_vec2gen(a,b,c,s,nrm)
    call d_rot2_vec2gen(a,b,c,s,nrm)
    call d_rot2_vec2gen(a,b,c,s,nrm)
    call d_rot2_vec2gen(a,b,c,s,nrm)
    call d_rot2_vec2gen(a,b,c,s,nrm)
  end do
  
  ! stop timer
  call system_clock(count=c_stop)
  t2 = dble(c_stop-c_start)/dble(c_rate)
  t2 = (t2 - t1)/dble(notests)/10d0

  ! print success
  call u_benchmark_print(t2,ERROR)

end program benchmark_d_rot2_vec2gen
