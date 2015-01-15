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
  integer, parameter :: N = 3
  real(8), parameter :: tol = 1d1*epsilon(1d0)
  integer :: ii, INFO
  logical :: P(N-2)
  real(8) :: Q(3*(N-1)), D(2,2*(N+1)), R(4,3*N)
  complex(8) :: V(N,N), W(N,N)
  complex(8) :: temp(2,2), A(2,2), B(2,2)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
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
    
    ! set valid Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii-2) = 1d0
    end do
    Q(1) = 0d0
    Q(2) = 0d0
    Q(3) = 1d0   
    temp(1,1) = cmplx(Q(1),Q(2),kind=8)
    temp(2,1) = cmplx(Q(3),0d0,kind=8)
    temp(1,2) = -temp(2,1)
    temp(2,2) = conjg(temp(1,1))
  
    ! set valid D
    D = 0d0
    do ii=1,(N+1)
      D(:,2*ii-1) = 1d0
    end do

    ! set valid R
    R = 0d0
    do ii=1,N
      R(1,3*ii) = -1d0
      R(2,3*ii) = 1d0
      R(3,3*ii) = -1d0
      R(4,3*ii) = 1d0
    end do
    
    ! set valid V
    V = cmplx(0d0,0d0,kind=8)
    do ii=1,N
      V(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
    
    ! set valid W
    W = V
    
    ! call twisted QZ
    call z_upr1fact_2x2deflation('QZ','I',N,1,P,Q,D,R,V,W,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if
    
    ! check results
    A(1,1) = cmplx(Q(1),Q(2),kind=8)
    A(2,1) = cmplx(Q(3),0d0,kind=8)
    A(1,2) = -A(2,1)
    A(2,2) = conjg(A(1,1))
    A(:,1) = A(:,1)*cmplx(D(1,1),D(1,2),kind=8)
    A(:,2) = A(:,2)*cmplx(D(1,3),D(1,4),kind=8)
    
    A = matmul(V(1:2,1:2),A)
    A = matmul(A,transpose(conjg(W(1:2,1:2))))
    
    B = cmplx(0d0,0d0,kind=8)
    B(1,1) = cmplx(D(2,1),D(2,2),kind=8)
    B(2,2) = cmplx(D(2,3),D(2,4),kind=8)
    
    B = matmul(V(1:2,1:2),B)
    B = matmul(B,transpose(conjg(W(1:2,1:2))))
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
