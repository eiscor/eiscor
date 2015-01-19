#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_scalar_nancheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_scalar_nancheck (nancheck). 
! The following tests are run:
!
! 1) test  0d0 + i0d0     
!
! 2) test 1.5d1 + i0d0            
!
! 3) test -2.42d24 + i1d0             
!
! 4) test  INF + i0d0
!
! 5) test  0d0 + iINF
!
! 6) test -INF + i0d0
!
! 7) test  0d0 - iINF
!
! 8) test  NAN + i0d0
!
! 9) test  0d0 + iNAN
!
! 10) test  NAN + i2d0
!
! 11) test  -2d0 + iNAN
!
! 12) test  huge(1d0) + 0d0
!
! 13) test  -huge(1d0) + 0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_scalar_nancheck

  implicit none

  ! compute variables
  real(8) :: a,b,nul
  complex(8) :: num
  integer :: info
  
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
  num = cmplx(a,b,kind=8)
  call z_scalar_nancheck(NUM,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  a = 1.5d1
  b = 0d0
  num = cmplx(a,b,kind=8)
  call z_scalar_nancheck(NUM,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  a = -2.42d24
  b = 1d0
  num = cmplx(a,b,kind=8)
  call z_scalar_nancheck(NUM,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  a = +huge(1d0)
  a = a+huge(1d0)
  b = 0d0
  num = cmplx(a,b,kind=8)
  call z_scalar_nancheck(NUM,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 5)

  ! set variables
  a = 0d0
  b = +huge(1d0)
  b = b+huge(1d0)
  num = cmplx(a,b,kind=8)
  call z_scalar_nancheck(NUM,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 6)

  ! set variables
  a = -huge(1d0)
  a = a-huge(1d0)
  b = 0d0
  num = cmplx(a,b,kind=8)
  call z_scalar_nancheck(NUM,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 7)

  ! set variables
  a = 0d0
  b = -huge(1d0)
  b = b-huge(1d0)
  num = cmplx(a,b,kind=8)
  call z_scalar_nancheck(NUM,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 8)

  ! set variables
  nul = 0d0
  a = nul/nul
  b = 0d0
  num = cmplx(a,b,kind=8)
  call z_scalar_nancheck(NUM,INFO)
  ! check info
  if (info.NE.-1) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 9)

  ! set variables
  nul = 0d0
  a = 0d0
  b = nul/nul
  num = cmplx(a,b,kind=8)
  call z_scalar_nancheck(NUM,INFO)
  ! check info
  if (info.NE.-1) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 10)

  ! set variables
  nul = 0d0
  a = nul/nul
  b = 2d0
  num = cmplx(a,b,kind=8)
  call z_scalar_nancheck(NUM,INFO)
  ! check info
  if (info.NE.-1) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 11)

  ! set variables
  nul = 0d0
  a = -2d0
  b = nul/nul
  num = cmplx(a,b,kind=8)
  call z_scalar_nancheck(NUM,INFO)
  ! check info
  if (info.NE.-1) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 12)

  ! set variables
  a = huge(1d0)
  b = 0d0
  num = cmplx(a,b,kind=8)
  call z_scalar_nancheck(NUM,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 13)

  ! set variables
  a = -huge(1d0)
  b = 0d0
  num = cmplx(a,b,kind=8)
  call z_scalar_nancheck(NUM,INFO)
  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_scalar_nancheck
