#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_extracttri
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_extracttri. 
! The following tests are run:
!
! 1) check roots of unity with upperhess QR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_extracttri

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2**2-1
  real(8) :: tol
  integer :: ii
  real(8) :: D(2*(N+1)), C(3*N), B(3*N)
  complex(8) :: T(N,N)
  
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
  
    ! set valid D
    D = 0d0
    do ii=1,(N+1)
      D(2*ii) = 1d0
    end do

    ! set valid C and B
    C = 0d0
    do ii=1,N
      C(3*ii) = -1d0
    end do
    B = -C
    
    ! extract diagonal
    call z_upr1fact_extracttri(.TRUE.,N,D,C,B,T)
   
!    print*,""
!    print*,"T" 
!    do ii=1,(N)
!      print*,T(ii,:)
!    end do
!    print*,""

    ! extract entire upper triangular part
    call z_upr1fact_extracttri(.FALSE.,N,D,C,B,T)
   
!    print*,""
!    print*,"T" 
!    do ii=1,(N)
!      print*,T(ii,:)
!    end do
!    print*,""

  ! end check 1)

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_upr1fact_extracttri
