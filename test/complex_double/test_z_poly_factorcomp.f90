#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_poly_factorcomp
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_poly_factorcomp. 
! The following tests are run:
!
! 1) check roots of unity with upperhess QR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_poly_factorcomp

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2**2
  real(8) :: tol
  integer :: ii
  logical :: P(N-2)
  real(8) :: Q(3*(N-1)), D1(2*(N+1)), C1(3*N), B1(3*N)
  real(8) :: D2(2*(N+1)), C2(3*N), B2(3*N)
  complex(8) :: COEFFS(N+1),V(N,N), W(N,N)
  
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
    COEFFS = cmplx(0d0,0d0,kind=8)
    COEFFS(1) = cmplx(1d0,0d0,kind=8)
    COEFFS(N+1) = cmplx(-1d0,0d0,kind=8)
    
    ! create factored companion matrix
    call z_poly_factorcomp(.FALSE.,.FALSE.,.FALSE.,N,COEFFS,P,Q,D1,C1,B1 &
    ,D2,C2,B2,V,W)
    
!    print*,""
!    print*,"Q"
!    do ii=1,(N-1)
!      print*,Q(3*ii-2),Q(3*ii-1),Q(3*ii)
!    end do
!    print*,""

!    print*,"D1"
!    do ii=1,(N+1)
!      print*,D1(2*ii-1),D1(2*ii)
!    end do
!    print*,""

!    print*,"C1"
!    do ii=1,(N)
!      print*,C1(3*ii-2),C1(3*ii-1),C1(3*ii)
!    end do
!    print*,""

!    print*,"B1"
!    do ii=1,(N)
!      print*,B1(3*ii-2),B1(3*ii-1),B1(3*ii)
!    end do
!    print*,""

  ! end check 1)

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
end program test_z_poly_factorcomp
