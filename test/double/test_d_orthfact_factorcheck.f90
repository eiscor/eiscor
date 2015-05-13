#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_orthfact_factorcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_d_orthfact_factorcheck.
! The following tests are run:
!
! 1) check with orthogonal Q and unit digaonal D
!    Q1=[0.8, 0.6], Q2=[0.5,sqrt(3)/2], Q3=[-1,eps/2], 
!    Q4=[-1/sqrt(2),1/sqrt(2)],
!    D = diag([1, -1, -1, 1, 1])
! 
! 2) perturb N
!
! 3) perturb Q1, then Q2, then Q3, then Q4
!
! 4) perturb D(1), D(2), D(3), D(4), D(5), D(6), D(7), D(8), D(9), D(10)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_orthfact_factorcheck

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*EISCOR_DBL_EPS ! accuracy (tolerance)
  
  ! compute variables
  integer :: INFO
  real(8) :: Q(8), Qs(8) ! N = 5
  real(8) :: D(5), Ds(5)
  
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
  D(2) = -1d0
  D(3) = -1d0
  D(4) = 1d0
  D(5) = 1d0

  Qs = Q
  Ds = D
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  Q = Qs
  D = Ds
  call d_orthfact_factorcheck(1,Q,D,INFO)

  ! check info
  if (INFO.NE.-1) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)
  Q = Qs
  D = Ds
  Q(1) = Q(1)+10*tol
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  D = Ds
  Q(2) = Q(2)-10*tol
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  D = Ds
  Q(3) = Q(3)+10*tol
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  D = Ds
  Q(4) = Q(4)-10*tol
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  D = Ds
  Q(5) = Q(5)+10*tol
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  D = Ds
  Q(6) = 0.1
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  D = Ds
  Q(7) = Q(7)+10*tol
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  D = Ds
  Q(8) = Q(8)-10*tol
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.-2) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)
  Q = Qs
  D = Ds
  D(1) = 0.2d0
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.-3) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  D = Ds
  D(2) = -2d0
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.-3) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  D = Ds
  D(3) = -.9d0
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.-3) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  D = Ds
  D(4) = 1d-300
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.-3) then
     call u_test_failed(__LINE__)
  end if

  Q = Qs
  D = Ds
  D(5) = 5d0
  call d_orthfact_factorcheck(5,Q,D,INFO)

  ! check info
  if (INFO.NE.-3) then
     call u_test_failed(__LINE__)
  end if


  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_d_orthfact_factorcheck
