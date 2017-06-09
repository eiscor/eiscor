#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_2Darray_check
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_2Darray_check (inf+nancheck). 
! The following tests are run:
!
! 1) test  0d0 + i0d0, 0d0 + i0d0      
!        1.5d1 + i0d0, 0d0 + i1.5d1            
!     -2.42d24 + i1d0, 1d0+ i-2.42d24             
!
! 2) test with one INF
!
! 3) test with one -INF
!
! 4) test with one NAN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_2Darray_check

  implicit none

  ! compute variables
  real(8) :: a,b,nul, num
  complex(8) :: C(3,2), C2(3,2)
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
  C(1,1) = cmplx(a,b,kind=8)
  C(1,2) = cmplx(a,b,kind=8)
  a = 1.5d1
  b = 0d0
  C(2,1) = cmplx(a,b,kind=8)
  a = 0d0
  b = 1.5d1
  C(2,2) = cmplx(a,b,kind=8)
  a = -2.42d24
  b = 1d0
  C(3,1) = cmplx(a,b,kind=8)
  a = 1d0
  b = -2.42d24
  C(3,2) = cmplx(a,b,kind=8)
  C2 = C
  call z_2Darray_check(3,2,C,flag)
  ! check flag
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  C = C2
  num = +huge(1d0)
  C(1,1) = cmplx(num+huge(1d0),0d0,kind=8)
  call z_2Darray_check(3,2,C,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  C = C2
  C(2,1) = cmplx(0d0,num+huge(1d0),kind=8)
  call z_2Darray_check(3,2,C,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  C = C2
  num = -huge(1d0)
  C(2,2) = cmplx(num-huge(1d0),0d0,kind=8)
  call z_2Darray_check(3,2,C,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  C = C2
  C(1,2) = cmplx(0d0,num-huge(1d0),kind=8)
  call z_2Darray_check(3,2,C,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  C = C2
  nul = 0d0
  C(3,2) = cmplx(nul/nul,0d0,kind=8)
  call z_2Darray_check(3,2,C,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  C = C2
  C(2,1) = cmplx(0d0,nul/nul,kind=8)
  call z_2Darray_check(3,2,C,flag)
  ! check flag
  if (flag) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_2Darray_check
