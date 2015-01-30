#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_2x2diagblocks
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_z_upr1fact_2x2diagblocks.
! The following tests are run:
!
! 1) P = 0, Q, D, R = I and ALG = QZ top and bottom block
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_2x2diagblocks

  implicit none
  
  ! parameter
  real(8), parameter :: tol = 1d1*EISCOR_DBL_EPS
  
  ! compute variables
  logical :: P
  real(8) :: Q(6), D1(4), C1(6), B1(6)
  real(8) :: D2(4), C2(6), B2(6)
  complex(8) :: A(2,2), B(2,2), eye(2,2)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! set eye
  eye(1,1) = cmplx(1d0,0d0,kind=8)
  eye(1,2) = cmplx(0d0,0d0,kind=8)
  eye(2,1) = cmplx(0d0,0d0,kind=8)
  eye(2,2) = cmplx(1d0,0d0,kind=8)
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 1)
  P = .FALSE.
  
  Q = 0d0
  Q(1) = 1d0; Q(4) = 1d0

  D1 = 0d0
  D1(1) = 1d0; D1(3) = 1d0
  D2 = D1
  
  C1 = 0d0
  C1(3) = 1d0; C1(6) = 1d0
  B1 = -C1; C2 = C1; B2 = B1

  ! top block
  call z_upr1fact_2x2diagblocks(.TRUE.,.TRUE.,.TRUE.,P,Q,D1,C1,B1,D2,C2,B2,A,B)

  ! check output
  if (maxval(abs(A-eye)) > tol) then
     call u_test_failed(__LINE__)
  end if
  if (maxval(abs(B-eye)) > tol) then
     call u_test_failed(__LINE__)
  end if
  
  ! bottom block
  call z_upr1fact_2x2diagblocks(.FALSE.,.TRUE.,.TRUE.,P,Q,D1,C1,B1,D2,C2,B2,A,B)
  
  ! check output
  if (maxval(abs(A-eye)) > tol) then
     call u_test_failed(__LINE__)
  end if
  if (maxval(abs(B-eye)) > tol) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_upr1fact_2x2diagblocks
