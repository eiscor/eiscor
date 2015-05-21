#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_comppenc_factor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_comppenc_factor. 
! The following tests are run:
!
! 1) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_comppenc_factor

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2**3
  real(8) :: tol
  integer :: ii, INFO
  logical :: P(N-2)
  real(8) :: Q(3*(N-1)), D1(2*(N+1)), C1(3*N), B1(3*N)
  real(8) :: D2(2*(N+1)), C2(3*N), B2(3*N)
  complex(8) :: V(N), W(N)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! set tolerance
  tol = 1d2*dble(N)*EISCOR_DBL_EPS

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! check 1)
    ! set INFO
    INFO = 0
    
    ! set P
    P = .FALSE.
    
    ! set valid V
    V = cmplx(2d0,1d0,kind=8)
    do ii=1,N
    end do     
  
    ! call twisted QZ
    call z_comppenc_factor(.FALSE.,N,P,V,W,Q,D1,C1,B1,D2,C2,B2,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

  ! end check 1)

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_comppenc_factor
