#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! benchmark_z_rfr3_turnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program benchmarks the subroutine z_rfr3_turnover. 
! The following benchmarks are run:
!
! 1) Compute ten million turnovers from uniformly generated 
!    user input. The average compute time and the worst error are printed.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program benchmark_z_rfr3_turnover

  implicit none

  ! parameter
  integer, parameter :: num_trials = 10**7 ! ten million trials

  ! compute variables
  integer :: ii, n
  integer, allocatable :: seed(:)
  complex(8) :: w, u, rho
  real(8) :: cc, ss, vv, nn
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
      
      ! set W
      call random_number(cc)
      call random_number(ss)
      w = cmplx(cc,ss,kind=8)
      w = w/abs(w)

      ! set U and VV
      call random_number(cc)
      call random_number(ss)
      call random_number(vv)
      vv = abs(vv)
      u = cmplx(cc,ss,kind=8)
      nn = abs(u)**2 + vv
      u = u/sqrt(nn)
      vv = vv/nn

      ! set RHO
      call random_number(cc)
      call random_number(ss)
      rho = cmplx(cc,ss,kind=8)
      rho = rho/abs(rho)

      ! set CC and SS
      call random_number(cc)
      call random_number(ss)
      cc = abs(cc)
      ss = abs(ss)
      nn = cc + ss
      cc = cc/nn
      ss = ss/nn

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
      
      ! set W
      call random_number(cc)
      call random_number(ss)
      w = cmplx(cc,ss,kind=8)
      w = w/abs(w)

      ! set U and VV
      call random_number(cc)
      call random_number(ss)
      call random_number(vv)
      vv = abs(vv)
      u = cmplx(cc,ss,kind=8)
      nn = abs(u)**2 + vv
      u = u/sqrt(nn)
      vv = vv/nn

      ! set RHO
      call random_number(cc)
      call random_number(ss)
      rho = cmplx(cc,ss,kind=8)
      rho = rho/abs(rho)

      ! set CC and SS
      call random_number(cc)
      call random_number(ss)
      cc = abs(cc)
      ss = abs(ss)
      nn = cc + ss
      cc = cc/nn
      ss = ss/nn

      ! compute 10 turnovers
      call z_rfr3_turnover(w,cc,ss,u,vv,rho)
      call z_rfr3_turnover(w,cc,ss,u,vv,rho)
      call z_rfr3_turnover(w,cc,ss,u,vv,rho)
      call z_rfr3_turnover(w,cc,ss,u,vv,rho)
      call z_rfr3_turnover(w,cc,ss,u,vv,rho)
      call z_rfr3_turnover(w,cc,ss,u,vv,rho)
      call z_rfr3_turnover(w,cc,ss,u,vv,rho)
      call z_rfr3_turnover(w,cc,ss,u,vv,rho)
      call z_rfr3_turnover(w,cc,ss,u,vv,rho)
      call z_rfr3_turnover(w,cc,ss,u,vv,rho)
      call z_rfr3_turnover(w,cc,ss,u,vv,rho)
      
    end do  

    ! stop timer
    call system_clock(count=c_stop)
    
    ! set total time
    total_time = dble(c_stop-c_start)/dble(c_rate)
    
    ! subract off base time and average
    total_time = (total_time - base_time)/dble(num_trials)/10d0
    
    ! print results
    call u_benchmark_print(total_time,ERROR)

end program benchmark_z_rfr3_turnover
