#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_unifact_factorcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_z_unifact_factorcheck.
! The following tests are run:
!
! Matrix:
! Q1=[0.8, 0.6, 0.0], Q2=[0.3,0.4,-sqrt(3)/2], 
! Q3=[-1/sqrt(3),1/sqrt(3), 1/sqrt(3)]
! D = diag([1, 0.8+i0,6, 1/sqrt(2)+i/sqrt(2), 1])
!
! 1) Matrix
! 2) N too small
! 3) Q with nan or inf
! 4) Q not orthogonal
! 5) D with nan or inf
! 6) D not unitary diagonal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_unifact_factorcheck

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*EISCOR_DBL_EPS ! accuracy (tolerance)
  
  ! compute variables
  integer :: N = 4, info
  real(8) :: Q(9), Qs(9), D(8), Ds(8), nul
  
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
  Q(3) = 0d0
  Q(4) = 3d-1
  Q(5) = 4d-1
  Q(6) = -sqrt(3d0)/2d0
  Q(7) = -1d0/sqrt(3d0)
  Q(8) = 1d0/sqrt(3d0)
  Q(9) = 1d0/sqrt(3d0)
  
  D(1) = 1d0
  D(2) = 0d0
  D(3) = 0.8d0
  D(4) = 0.6d0
  D(5) = 1d0/sqrt(2d0)
  D(6) = 1d0/sqrt(2d0)
  D(7) = 1d0
  D(8) = 0d0

  nul = 0d0

  Qs = Q
  Ds = D

  call z_unifact_factorcheck(N,Q,D,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  
  call z_unifact_factorcheck(1,Q,D,INFO)

  ! check info
  if (INFO.NE.-1) then
     call u_test_failed(__LINE__)
  end if

  call z_unifact_factorcheck(-10,Q,D,INFO)

  ! check info
  if (INFO.NE.-1) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)
  Q(1) = 1d0/nul
  call z_unifact_factorcheck(N,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  Q(2) = -1d0/nul
  call z_unifact_factorcheck(N,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  Q(3) = nul/nul
  call z_unifact_factorcheck(N,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  Q(4) = -nul/nul
  call z_unifact_factorcheck(N,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)
  Q = Qs
  Q(1) = 2d0
  call z_unifact_factorcheck(N,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  Q(2) = -0d0
  call z_unifact_factorcheck(N,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  Q(3) = 0.1d0
  call z_unifact_factorcheck(N,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  Q(4) = 20d0
  call z_unifact_factorcheck(N,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 5)
  Q = Qs
  D(1) = -1d0/nul
  call z_unifact_factorcheck(N,Q,D,INFO)

  ! check info
  if (INFO.NE.-3) then
     call u_test_failed(__LINE__)
  end if

  D = Ds
  D(2) = +1d0/nul
  call z_unifact_factorcheck(N,Q,D,INFO)

  ! check info
  if (INFO.NE.-3) then
     call u_test_failed(__LINE__)
  end if

  D = Ds
  D(3) = nul/nul
  call z_unifact_factorcheck(N,Q,D,INFO)

  ! check info
  if (INFO.NE.-3) then
     call u_test_failed(__LINE__)
  end if
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 6)
  D = Ds
  D(1) = 0.9
  call z_unifact_factorcheck(N,Q,D,INFO)

  ! check info
  if (INFO.NE.-3) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_unifact_factorcheck
