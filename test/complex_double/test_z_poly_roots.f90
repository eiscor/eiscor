#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_poly_roots
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_comppenc_factor. 
! The following tests are run:
!
! 1) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_poly_roots

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2**5
  real(8) :: tol
  integer :: ii, INFO
  real(8) :: RESIDUALS(N)
  complex(8) :: COEFFS(N+1),ROOTS(N)
  
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
    ! set valid COEFFS
    COEFFS = cmplx(0d0,0d0,kind=8)
    COEFFS(1) = cmplx(1d0,0d0,kind=8)
    COEFFS(N+1) = cmplx(-1d0,0d0,kind=8)
  
    ! call roots
    call z_poly_roots(N,COEFFS,ROOTS,RESIDUALS)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

  ! end check 1)

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_poly_roots
