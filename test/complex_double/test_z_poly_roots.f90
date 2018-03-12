#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_poly_roots
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_comppenc_factor. 
! The following tests are run:
!
! 1) roots of unity
! 2) roots of unity modulo 1
! 3) degree 8 wilkinson polynomial
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_poly_roots

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2**3
  real(8), parameter :: pi = 3.14159265358979323846264338327950d0
  real(8) :: tol, scl
  integer :: ii, INFO, ARGS(N)
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
  call z_poly_roots(N,COEFFS,ROOTS,RESIDUALS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    call u_test_failed(__LINE__)
  end if
 
  ! check residuals
  do ii=1,N
    if (RESIDUALS(ii) >= tol) then
      call u_test_failed(__LINE__)
    end if
  end do
 
  ! end check 1)




  ! check 2)
  ! set valid COEFFS
  COEFFS = cmplx(1d0,0d0,kind=8)
 
  ! call roots
  call z_poly_roots(N,COEFFS,ROOTS,RESIDUALS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    call u_test_failed(__LINE__)
  end if

  ! check residuals
  do ii=1,N
    if (RESIDUALS(ii) >= tol) then
      call u_test_failed(__LINE__)
    end if
  end do
 
  ! end check 2)
 
 
  ! check 3)
  ! set valid COEFFS
  COEFFS(1) = cmplx(1d0,0d0,kind=8)         
  COEFFS(2) = cmplx(-36d0,0d0,kind=8)         
  COEFFS(3) = cmplx(546d0,0d0,kind=8)       
  COEFFS(4) = cmplx(-4536d0,0d0,kind=8)       
  COEFFS(5) = cmplx(22449d0,0d0,kind=8)      
  COEFFS(6) = cmplx(-67284d0,0d0,kind=8)      
  COEFFS(7) = cmplx(118124d0,0d0,kind=8)     
  COEFFS(8) = cmplx(-109584d0,0d0,kind=8)       
  COEFFS(9) = cmplx(40320d0,0d0,kind=8)

  ! call roots
  call z_poly_roots(N,COEFFS,ROOTS,RESIDUALS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    call u_test_failed(__LINE__)
  end if

  ! check residuals
  do ii=1,N
    if (RESIDUALS(ii) >= 1d3*tol) then
      call u_test_failed(__LINE__)
    end if
  end do

  ! end check 3)


  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_poly_roots
