#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_poly_badpoly
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_poly_badpoly. 
! The following tests are run:
!
! 1) check roots of unity with upperhess QR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_poly_badpoly

  implicit none
  
  ! compute variables
  integer :: ii, N
  real(8) :: tol, a, b
  real(8), allocatable :: RESIDUALS(:)
  complex(8), allocatable :: COEFFS(:), ROOTS(:)
  
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

    ! open file
    open(unit=10, file='/Users/jared/badpoly.txt', status='unknown')

    ! read in degree
    read(10,*) N

    ! allocate memory
    allocate(RESIDUALS(N),ROOTS(N),COEFFS(N+1))

    ! read in coefficients
    do ii=1,(N+1)
      read(10,*) COEFFS(ii)
    end do

    ! close file 
    close(10)
  
    ! compute roots
    call z_poly_roots(N,COEFFS,ROOTS,RESIDUALS)
   
!    print*,""
!    print*,"ROOTS"
!    do ii=1,(N)
!      print*,ROOTS(ii),abs(ROOTS(ii)),RESIDUALS(ii)
!    end do
!    print*,""

  ! free memory
  deallocate(RESIDUALS,ROOTS,COEFFS)

  ! end check 1)

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_poly_badpoly
