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
! 1) Time required for 1 billion runs.
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program benchmark_d_rot2_vec2gen

  implicit none

  ! parameter
  !integer, parameter :: notests = 1000000000 ! 1 billion
  integer, parameter :: notests = 10000000 ! 10 million
  !integer, parameter :: notests = 1000000 ! 1 million
  
  ! compute variables
  real(8) :: nrm
  real(8) :: a, b
  real(8) :: c, s
  real(8) :: t1, t2

    
  ! timing variables
  integer:: c_start, c_stop, c_rate, ii


  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! just random number generation
  do ii=1,notests
     call random_number(a)
     call random_number(b)
  end do

  ! stop timer
  call system_clock(count=c_stop)
  t1 = dble(c_stop-c_start)/dble(c_rate)

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! random number generation + rotation generation
  do ii=1,notests
     call random_number(a)
     call random_number(b)
     call d_rot2_vec2gen(a,b,c,s,nrm)
  end do
  
  ! stop timer
  call system_clock(count=c_stop)
  t2 = dble(c_stop-c_start)/dble(c_rate) - t1

  ! print success
  call u_benchmark(t2/notests,0d0)

end program benchmark_d_rot2_vec2gen
