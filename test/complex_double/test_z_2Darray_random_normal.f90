#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_2Darray_random_normal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_z_2Darray_random_normal.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_2Darray_random_normal

  implicit none

  ! compute variables
  complex(8) :: A(3,2)
  integer :: info
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)

  call u_fixedseed_initialize(info)

  ! check info
  if (info.NE.0) then
     call u_test_failed(__LINE__)
  end if

  call z_2Darray_random_normal(3,2,A)

  if (abs(A(1,1)-cmplx(0.47761357716530684d0, -1.0710694263009257d0, kind=8))>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if
  
  if (abs(A(2,1)-cmplx(0.75081696949437793d0, 0.88892563777904299d0, kind=8))>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if
  
  if (abs(A(3,1)-cmplx(1.0977344437214227d0, -1.0886786789374205d0, kind=8))>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if

  if (abs(A(1,2)-cmplx(7.38055464269166406d-2, 1.2244927626472537d0, kind=8))>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if
  
  if (abs(A(2,2)-cmplx(-0.79659862878327770d0, 1.0930703380490570d0, kind=8))>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if
  
  if (abs(A(3,2)-cmplx(-1.6229037583856123d0, -0.44858989349582140d0, kind=8))>100*EISCOR_DBL_EPS) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
end program test_z_2Darray_random_normal
