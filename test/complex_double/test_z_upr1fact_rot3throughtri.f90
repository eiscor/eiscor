#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_rot3throughtri
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_rot3throughtri. 
! The following tests are run (once with 'L2R' and once with 'R2L'):
!
! 1)  D = I, C = B = all [0, -1; 1, 0] and G = [1/sqrt(2) 0 1/sqrt(2)]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_rot3throughtri

  implicit none

  ! parameter
  integer, parameter :: N = 2
  real(8), parameter :: tol = 2d0*epsilon(1d0) ! accuracy (tolerance)

  ! compute variables
  integer :: ii, INFO
  real(8) :: Dold(2*(N+1)), Cold(3*N), Bold(3*N), Gold(3)
  real(8) :: D(2*(N+1)), C(3*N), B(3*N), G(3)

  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)
  
  ! set INFO
  INFO = 0

  !!!!!!!!!!!!!!!!!!!!
  ! check 1)
 
  ! set variables
  D = 0d0
  do ii=1,(N+1)
    D(2*ii-1) = 1d0
  end do

  C = 0d0
  do ii=1,N
    C(3*ii) = 1d0
  end do
  
  B = C
  
  G(1) = 1d0/sqrt(2d0)
  G(2) = 0d0
  G(3) = G(1)
  
  Dold = D
  Cold = C
  Bold = B
  Gold = G
  
  ! call 
  call z_upr1fact_rot3throughtri('L2R',N,1,D,C,B,G,INFO)
  
  ! check info
  if (INFO.NE.0) then
    call u_test_failed(__LINE__)
  end if
  
  ! call 
  call z_upr1fact_rot3throughtri('R2L',N,1,D,C,B,G,INFO)
  
  ! check info
  if (INFO.NE.0) then
    call u_test_failed(__LINE__)
  end if

  ! check results
  if (maxval(abs(Dold-D)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(Cold-C)) > tol) then
    call u_test_failed(__LINE__)
  end if  
  if (maxval(abs(Bold-B)) > tol) then
    call u_test_failed(__LINE__)
  end if  
  if (maxval(abs(Gold-G)) > tol) then
    call u_test_failed(__LINE__)
  end if    

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_upr1fact_rot3throughtri
