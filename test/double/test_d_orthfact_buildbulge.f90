#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_orthfact_buildbulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_d_orthfact_buildbulge
! The following tests are run:
!
! 1) check with orthogonal Q and unit digaonal D
!    Q1=[0.8, 0.6], Q2=[0.5,sqrt(3)/2], Q3=[-1,eps/2], 
!    Q4=[-1/sqrt(2),1/sqrt(2)],
!    D = diag([1, -1, -1, 1, 1])
!    single shift 1.0
!
! 2) check with orthogonal Q and unit digaonal D
!    Q1=[0.8, 0.6], Q2=[0.5,sqrt(3)/2], Q3=[-1,eps/2], 
!    Q4=[-1/sqrt(2),1/sqrt(2)],
!    D = diag([1, -1, -1, 1, 1])
!    double shift 0.5 +- i*sqrt(3)/2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Notes:
!
! A = [ 0.8   0.3     ... ]
!     [ 0.6  -0.4     ... ]
!     [  0  sqrt(3)/2 ... ]
!     [ ...    0      ... ]
!
! mu = alpha + i beta
! A^2e_1 - 2alpha Ae_1 + (alpha^2 + beta^2)e_1 = 
! = [  1.02        ]
!   [ -0.36        ]
!   [  0.3*sqrt(3) ]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_orthfact_buildbulge

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*epsilon(1d0) ! accuracy (tolerance)
  
  ! compute variables
  integer :: info, str, stp, zero, itcnt, its(4)
  real(8) :: Q(8), Qs(8) ! N = 5
  real(8) :: D(10), Ds(10)
  real(8) :: E(2), Es(2)
  real(8) :: B1(2), B2(2)

  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 1)
  Q(1) = 8d-1
  Q(2) = 6d-1
  Q(3) = 5d-1
  Q(4) = -sqrt(3d0)/2d0
  Q(5) = -1d0
  Q(6) = tol/2d1
  Q(7) = -1d0/sqrt(2d0)
  Q(8) = 1d0/sqrt(2d0)
  
  D(1) = 1d0
  D(2) = 0d0
  D(3) = -1d0
  D(4) = 0d0
  D(5) = -1d0
  D(6) = 0d0
  D(7) = 1d0
  D(8) = 0d0
  D(9) = 1d0
  D(10) = 0d0

  E(1) = 1d0
  E(2) = 0d0

  Qs = Q
  Ds = D
  Es = E
  call d_orthfact_buildbulge('S',5,1,Q,D,E,B1,B2,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  
  ! check B1
  if (abs(B1(1)+1d0/sqrt(10d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(B1(2)-3d0/sqrt(10d0))>tol) then
     call u_test_failed(__LINE__)
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  Q = Qs
  D = Ds
  E(1) = 5d-1
  E(2) = sqrt(3d0)/2d0
  call d_orthfact_buildbulge('D',5,1,Q,D,E,B1,B2,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  
  ! check B1
  if (abs(B1(1)+0.36d0/sqrt(0.3996d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(B1(2)-0.3d0/sqrt(0.3996d0)*sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  ! check B2
  if (abs(B2(1)-1.02d0/1.2d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(B2(2)-sqrt(0.3996d0)/1.2d0)>tol) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_d_orthfact_buildbulge
