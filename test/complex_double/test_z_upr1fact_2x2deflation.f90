#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_2x2deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_2x2deflation. 
! The following tests are run:
!
! 1) 1st factor in Q is [0, -1; 1, 0] and both triangular parts are I
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_2x2deflation

  implicit none
  
  ! compute variables
  real(8), parameter :: tol = 1d1*EISCOR_DBL_EPS
  real(8) :: Q(3), D1(4), C1(6), B1(6)
  real(8) :: D2(4), C2(6), B2(6)
  complex(8) :: V(2,2), W(2,2)
  complex(8) :: temp(2,2), A(2,2), B(2,2)
  
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
  
    ! set valid D1 and D2
    D1 = 0d0
    D1(1) = 1d0; D1(3) = 1d0
    D2 = D1

    ! set valid C1, B1, C2 and B2
    C1 = 0d0
    C1(3) = 1d0; C1(6) = 1d0
    C2 = C1; B1 = -C1; B2 = B1
 
    ! set valid V
    V = cmplx(0d0,0d0,kind=8)
    V(1,1) = cmplx(1d0,0d0,kind=8)
    V(2,2) = cmplx(1d0,0d0,kind=8)
    
    ! set valid W
    W = V
    
    ! call twisted QZ
    call z_upr1fact_2x2deflation(.TRUE.,.TRUE.,Q,D1,C1,B1,D2,C2,B2,2,V,W)
    
    ! check results
    A(1,1) = cmplx(Q(1),Q(2),kind=8)
    A(2,1) = cmplx(Q(3),0d0,kind=8)
    A(1,2) = -A(2,1)
    A(2,2) = conjg(A(1,1))

    A(:,1) = A(:,1)*cmplx(D1(1),D1(2),kind=8)
    A(:,2) = A(:,2)*cmplx(D1(3),D1(4),kind=8)
    
    A = matmul(V,A)
    A = matmul(A,transpose(conjg(W)))
    
    B = cmplx(0d0,0d0,kind=8)
    B(1,1) = cmplx(D2(1),D2(2),kind=8)
    B(2,2) = cmplx(D2(3),D2(4),kind=8)
    
    B = matmul(V,B)
    B = matmul(B,transpose(conjg(W)))
    A = matmul(transpose(conjg(B)),A)
   
    if (maxval(abs(A-temp)) > tol) then
      call u_test_failed(__LINE__)
    end if
  ! end check 1)
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
end program test_z_upr1fact_2x2deflation
