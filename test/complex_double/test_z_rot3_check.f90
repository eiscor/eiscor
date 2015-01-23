#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_rot3_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_rot3_check. 
! The following tests are run:
!
! 0) contains an INF or NAN
!
! 1) CR = 1; CI = 0; S = 0                
!    and replace 0 by EISCOR_DBL_EPS 
!
! 2) CR = cos(t1)*cos(t2); CI = sin(t1)*cos(t2); S = sin(t2)                
!    for t1, t2 uniform random in [0,2pi) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rot3_check

  implicit none

  ! parameter
  real(8), parameter :: tol = EISCOR_DBL_EPS ! accuracy (tolerance)
  real(8), parameter :: two_pi = 2d0*EISCOR_DBL_PI

  ! compute variables
  logical :: FLAG
  integer :: ii, n, num_tests = 10**6
  integer, allocatable :: seed(:)
  real(8) :: CR, CI, S, NRM
  real(8) :: t1, t2, inf, nan
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)
  
  ! set inf
  inf = EISCOR_DBL_INF
  inf = 10d0*inf
  
  ! set nan
  nan = 0d0
  nan = 0d0/nan
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 0)

    ! one INF
    CR = 1d0; CI = 1d0; S = inf
    call z_rot3_check(CR,CI,S,FLAG)
    if (FLAG) then
      call u_test_failed(__LINE__)
    end if
    
    ! one NAN
    CR = 1d0; CI = nan; S = 0d0
    call z_rot3_check(CR,CI,S,FLAG)
    if (FLAG) then
      call u_test_failed(__LINE__)
    end if
    
  !!!!!!!!!!!!!!!!!!!!
  ! check 1)

    ! identity
    CR = 1d0; CI = 0d0; S = 0d0
    call z_rot3_check(CR,CI,S,FLAG)
    if (.NOT.FLAG) then
      call u_test_failed(__LINE__)
    end if
    
    ! nearly identity
    CR = 1d0; CI = tol; S = tol
    call z_rot3_check(CR,CI,S,FLAG)
    if (.NOT.FLAG) then
      call u_test_failed(__LINE__)
    end if
    
  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  
    ! get size of see        
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
    
    ! random tests
    do ii=1,num_tests
      call random_number(t1)
      call random_number(t2)
      CR = cos(t1*two_pi)*cos(t2*two_pi)
      CI = sin(t1*two_pi)*cos(t2*two_pi)
      S = sin(t2*two_pi)
      call z_rot3_vec3gen(CR,CI,S,CR,CI,S,NRM)
      call z_rot3_check(CR,CI,S,FLAG)
      if (.NOT.FLAG) then
        call u_test_failed(__LINE__)
      end if
    end do
  
  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rot3_check
