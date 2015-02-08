#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_poly_roots
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_poly_roots. 
! The following tests are run:
!
! 1) check roots of unity with upperhess QR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_poly_roots

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2**6
  real(8) :: tol, a, b,  RESIDUALS(N)
  integer :: ii
  complex(8) :: COEFFS(N+1), ROOTS(N)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! set tolerance
  tol = 1d1*dble(N)*EISCOR_DBL_EPS

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! check 1)
  
    ! set valid COEFFS
!    COEFFS = cmplx(0d0,0d0,kind=8)
!    COEFFS(1) = cmplx(1d0,0d0,kind=8)
!    COEFFS(N+1) = cmplx(-1d0,0d0,kind=8)

!    COEFFS = cmplx(1d0,0d0,kind=8)

    call u_randomseed_initialize()
    do ii=1,(N+1)
      call random_number(a)
      call random_number(b)
      COEFFS(ii) = cmplx(a,b,kind=8)
    end do

    ! extract diagonal
    call z_poly_roots(N,COEFFS,ROOTS,RESIDUALS)
   
    print*,""
    print*,"ROOTS"
    do ii=1,(N)
      print*,ROOTS(ii),abs(ROOTS(ii)),RESIDUALS(ii)
    end do
    print*,""

!    print*,""
!    print*,"worst residual"
!    print*,maxval(RESIDUALS)
!    print*,""

  ! end check 1)

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_poly_roots
