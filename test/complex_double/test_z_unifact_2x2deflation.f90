#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_unifact_2x2deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_unifact_2x2deflation. 
! The following tests are run:
!
! 1) 1st factor in Q is [0, -1; 1, 0] and D is the identity
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_unifact_2x2deflation

  implicit none
  
  ! compute variables
  real(8), parameter :: tol = 1d1*EISCOR_DBL_EPS
  real(8) :: Q(3), D(4)
  complex(8) :: Z(2,2)
  complex(8) :: temp(2,2), A(2,2)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! check 1)
    ! set valid Q
    Q = 0d0
    Q(3) = 1d0
    temp(1,1) = cmplx(Q(1),Q(2),kind=8)
    temp(2,1) = cmplx(Q(3),0d0,kind=8)
    temp(1,2) = -temp(2,1)
    temp(2,2) = conjg(temp(1,1))
  
    ! set valid D
    D = 0d0
    D(1) = 1d0; D(3) = 1d0

    ! set valid Z
    Z = cmplx(0d0,0d0,kind=8)
    Z(1,1) = cmplx(1d0,0d0,kind=8)
    Z(2,2) = cmplx(1d0,0d0,kind=8)
    
    ! 2x2 deflation
    call z_unifact_2x2deflation(.TRUE.,Q,D,2,Z)

    ! check results
    A(1,1) = cmplx(Q(1),Q(2),kind=8)
    A(2,1) = cmplx(Q(3),0d0,kind=8)
    A(1,2) = -A(2,1)
    A(2,2) = conjg(A(1,1))

    A(:,1) = A(:,1)*cmplx(D(1),D(2),kind=8)
    A(:,2) = A(:,2)*cmplx(D(3),D(4),kind=8)
    
    A = matmul(Z,A)
    A = matmul(A,transpose(conjg(Z)))
    
    if (maxval(abs(A-temp)) > tol) then
      call u_test_failed(__LINE__)
    end if
  ! end check 1)
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
end program test_z_unifact_2x2deflation
