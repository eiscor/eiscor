#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_factorcompmat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_factorcompmat. 
! The following tests are run:
!
! 1) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_factorcompmat

  implicit none
  
  ! compute variables
  integer, parameter :: N = 3
  integer :: ii, INFO, ITS(N-1)
  logical :: P(N-2)
  real(8) :: Q(3*N), D1(2*(N+1)), D2(2*(N+1))
  real(8) :: C1(3*N), B1(3*N), C2(3*N), B2(3*N)
  complex(8) :: A1(N), A2(N)
  complex(8) :: V(N,N), W(N,N)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! set INFO
  INFO = 0
  
  ! check 1)
    
    ! set A1   
    A1 = cmplx(1d0,0d0,kind=8)
    
    ! set P
    P = .FALSE.
    
    ! call factor compmat
    call z_upr1fact_factorcompmat('QR','N',N,A1,A2,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if
  ! end check 1)
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_upr1fact_factorcompmat
