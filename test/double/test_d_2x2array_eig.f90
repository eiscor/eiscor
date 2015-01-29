#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_2x2array_eig
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_d_2x2array_eig. 
! The following tests are run:
!
! 1) A = [1,2;2,1], B = I
!
! 2) A = [0.5, 0.3; 0.2, 0.3], B = I
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_2x2array_eig

  implicit none
  
  ! parameter
  real(8), parameter :: tol = 1d1*EISCOR_DBL_EPS ! accuracy (tolerance)
  
  ! compute variables
  real(8) :: A(2,2), B(2,2), Q(2,2), Z(2,2)
  real(8) :: Aold(2,2), Bold(2,2)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 1)
 
    ! set input
    A(1,1) = 1d0
    A(2,1) = 2d0
    A(1,2) = 2d0
    A(2,2) = 1d0
    Aold = A
    
    call d_2x2array_eig(.FALSE.,A,B,Q,Z)
    
print*,""
print*,"Aold"
print*,Aold(1,:)
print*,Aold(2,:)
print*,""

print*,"A"
print*,A(1,:)
print*,A(2,:)
print*,""

print*,"Q"
print*,Q(1,:)
print*,Q(2,:)
print*,""
    
    ! check results
    Aold = A - matmul(transpose(Q),matmul(Aold,Q))
    
print*,"Aold"
print*,Aold(1,:)
print*,Aold(2,:)
print*,""
    
    if (maxval(abs(Aold))>tol) then
      call u_test_failed(__LINE__)
    end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  
    ! set input
    A(1,1) = 5d0
    A(2,1) = -2d0
    A(1,2) = 3d0
    A(2,2) = 3d0
    Aold = A
    
    call d_2x2array_eig(.FALSE.,A,B,Q,Z)
    
print*,""
print*,"Aold"
print*,Aold(1,:)
print*,Aold(2,:)
print*,""

print*,"A"
print*,A(1,:)
print*,A(2,:)
print*,""

print*,"Q"
print*,Q(1,:)
print*,Q(2,:)
print*,""
    
    ! check results
    Aold = matmul(Aold,Q) 
    A = matmul(Q,A)
    
print*,"Aold"
print*,Aold(1,:)
print*,Aold(2,:)
print*,""

print*,"A"
print*,A(1,:)
print*,A(2,:)
print*,""
    
    if (maxval(abs(Aold))>tol) then
      call u_test_failed(__LINE__)
    end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_d_2x2array_eig
