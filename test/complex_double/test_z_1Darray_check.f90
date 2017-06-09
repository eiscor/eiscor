#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_1Darray_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_1Darray_check (inf+nancheck). 
! The following tests are run:
!
! 1) test  0d0 + i0d0, 0d0 + i0d0      
!
! 2) test 1.5d1 + i0d0, 0d0 + i1.5d1            
!
! 3) test -2.42d24 + i1d0, 1d0+ i-2.42d24             
!
! 4) test  INF + i0d0, 0d0 + i0d0
!
! 5) test  0d0 + iINF, 0d0 + i0d0
!
! 6) test  0d0 + i0d0, -INF + i0d0
!
! 7) test  0d0 + i0d0, 0d0 - iINF
!
! 8) test  0d0 + i0d0, NAN + i0d0
!
! 9) test  0d0 + iNAN, 0d0 + i0d0
!
! 10) test  0d0 + i0d0, huge(1d0) + 0d0
!
! 11) test  -huge(1d0) + 0d0, 0d0 + i0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_1Darray_check

  implicit none

  ! compute variables
  real(8) :: a,b,nul
  complex(8) :: C(2)
  complex(8) :: num
  logical :: flag
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)

  !!!!!!!!!!!!!!!!!!!!
  ! check 1)

  ! set variables
  a = 0d0
  b = 0d0
  C(1) = cmplx(a,b,kind=8)
  C(2) = cmplx(a,b,kind=8)
  call z_1Darray_check(2,C,flag)
  ! check flag
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  a = 1.5d1
  b = 0d0
  C(1) = cmplx(a,b,kind=8)
  a = 0d0
  b = 1.5d1
  C(2) = cmplx(a,b,kind=8)
  call z_1Darray_check(2,C,flag)
  ! check flag
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  a = -2.42d24
  b = 1d0
  C(1) = cmplx(a,b,kind=8)
  a = 1d0
  b = -2.42d24
  C(2) = cmplx(a,b,kind=8)
  call z_1Darray_check(2,C,flag)
  ! check flag
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  a = +huge(1d0)
  a = a+huge(1d0)
  b = 0d0
  C(1) = cmplx(a,b,kind=8)
  a = 0d0
  b = 0d0
  C(2) = cmplx(a,b,kind=8)
  call z_1Darray_check(2,C,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 5)

  ! set variables
  a = 0d0
  b = +huge(1d0)
  b = b+huge(1d0)
  C(1) = cmplx(a,b,kind=8)
  a = 0d0
  b = 0d0
  C(2) = cmplx(a,b,kind=8)
  call z_1Darray_check(2,C,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 6)

  ! set variables
  a = 0d0
  b = 0d0
  C(1) = cmplx(a,b,kind=8)
  a = -huge(1d0)
  a = a-huge(1d0)
  b = 0d0
  C(2) = cmplx(a,b,kind=8)
  call z_1Darray_check(2,C,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 7)

  ! set variables
  a = 0d0
  b = 0d0
  C(1) = cmplx(a,b,kind=8)
  a = 0d0
  b = -huge(1d0)
  b = b-huge(1d0)
  C(2) = cmplx(a,b,kind=8)
  call z_1Darray_check(2,C,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 8)

  ! set variables 
  a = 0d0
  b = 0d0
  C(1) = cmplx(a,b,kind=8)
  nul = 0d0
  a = nul/nul
  b = 0d0
  C(2) = cmplx(a,b,kind=8)
  call z_1Darray_check(2,C,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 9)

  ! set variables
  nul = 0d0
  a = 0d0
  b = nul/nul
  C(1) = cmplx(a,b,kind=8)
  a = 0d0
  b = 0d0
  C(2) = cmplx(a,b,kind=8)
  call z_1Darray_check(2,C,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 10)

  ! set variables
  a = 0d0
  b = 0d0
  C(1) = cmplx(a,b,kind=8)
  a = huge(1d0)
  b = 0d0
  C(2) = cmplx(a,b,kind=8)
  call z_1Darray_check(2,C,flag)
  ! check flag
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 11)

  ! set variables
  a = -huge(1d0)
  b = 0d0
  C(1) = cmplx(a,b,kind=8)
  a = 0d0
  b = 0d0
  C(2) = cmplx(a,b,kind=8)
  call z_1Darray_check(2,C,flag)
  ! check flag
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_1Darray_check
