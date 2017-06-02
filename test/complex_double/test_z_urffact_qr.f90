#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_urffact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_urffact_qr. The following tests are run:
!
! 1) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_urffact_qr

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2**13
  integer :: ii, jj, INFO, ITCNT(N-1)
  complex(8) :: U(N)
  real(8) :: VV(N)
  real(8) :: tol
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  
 
  ! Check 1)
  ! initialize U and VV
  U = cmplx(0d0,0d0,kind=8)
  U(N) = cmplx(sign(1d0,(-1d0)**(N-1)),0d0,kind=8)
  VV = 1d0
  VV(N) = 0d0
!print*,""
!print*,"U,VV"
!do ii = 1,N
!  print*,U(ii),VV(ii)
!end do
    
  ! call singlestep
  call z_urffact_rfqr(N,U,VV,ITCNT,INFO)
!print*,""
!print*,"U,VV"
!print*,U(1),VV(1)
!do ii = 2,N
!  print*,U(ii),VV(ii),ITCNT(ii-1)
!end do


  ! set tolerance
  tol = 10d0*EISCOR_DBL_EPS
    
  ! check maximum entry
!  if (maxval(abs(H(1:2,1:2))) >= tol) then
!    call u_test_failed(__LINE__)
!  end if  

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_urffact_qr
