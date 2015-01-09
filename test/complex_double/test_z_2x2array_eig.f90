#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_2x2array_eig
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_z_2x2array_eig. 
! The following tests are run:
!
! 1) H = [1, 2; 2, 1] +i0d0 (example from test_d_2x2array_eig)
!
! 2) H = [0.5, 0.3; 0.2, 0.3] +i0d0 (example from test_d_2x2array_eig)
!
! 3) H = [1+i, 1; i, 3]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_2x2array_eig

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*epsilon(1d0) ! accuracy (tolerance)
  
  ! compute variables
  integer :: info
  complex(8) :: H(2,2)
  complex(8) :: E(2), Z(2,2)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 1)

  H(1,1) = cmplx(1d0,0d0,kind=8)
  H(2,1) = cmplx(2d0,0d0,kind=8)
  H(1,2) = cmplx(2d0,0d0,kind=8)
  H(2,2) = cmplx(1d0,0d0,kind=8)
  
  call z_2x2array_eig(H,E,Z,INFO)
  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  
  ! check results
  ! check E
  if ((abs(E(1)-3d0)>tol).AND.(abs(E(2)-3d0)>tol)) then
    call u_test_failed(__LINE__)
  end if
  if ((abs(E(1)+1d0)>tol).AND.(abs(E(2)+1d0)>tol)) then
    call u_test_failed(__LINE__)
  end if

  ! check Z
  if (abs(E(1)-3d0)>tol) then ! E(2) = 3
     if (abs(Z(1,1)-cmplx(sqrt(2d0)/2d0,0d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if
     if (abs(Z(2,1)+cmplx(sqrt(2d0)/2d0,0d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if
     if (abs(Z(1,2)-cmplx(sqrt(2d0)/2d0,0d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if     
     if (abs(Z(2,2)-cmplx(sqrt(2d0)/2d0,0d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if     
  else
     if (abs(Z(1,1)-cmplx(sqrt(2d0)/2d0,0d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if
     if (abs(Z(2,1)-cmplx(sqrt(2d0)/2d0,0d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if
     if (abs(Z(1,2)+cmplx(sqrt(2d0)/2d0,0d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if     
     if (abs(Z(2,2)-cmplx(sqrt(2d0)/2d0,0d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if     
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  H(1,1) = cmplx(5d-1,0d0,kind=8)
  H(2,1) = cmplx(-2d-1,0d0,kind=8)
  H(1,2) = cmplx(3d-1,0d0,kind=8)
  H(2,2) = cmplx(3d-1,0d0,kind=8)

  call z_2x2array_eig(H,E,Z,INFO)
  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
 
  ! check results
  ! check E
  if ((abs(E(1)-cmplx(4d-1,sqrt(5d-2),kind=8))>tol).AND.(abs(E(2)-cmplx(4d-1,sqrt(5d-2),kind=8))>tol)) then
    call u_test_failed(__LINE__)
  end if
  if ((abs(E(1)-cmplx(4d-1,-sqrt(5d-2),kind=8))>tol).AND.(abs(E(2)-cmplx(4d-1,-sqrt(5d-2),kind=8))>tol)) then
    call u_test_failed(__LINE__)
  end if

  ! check Z
  if (abs(E(1)-(4d-1+sqrt(5d-2)))>tol) then ! E(2) = (4d-1+sqrt(5d-2))
     if (abs(Z(1,1)-cmplx(1d0/sqrt(1d1),-sqrt(2d0)/2d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if
     if (abs(Z(1,2)-cmplx(sqrt(2d0)/sqrt(5d0),0d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if     
     if (abs(Z(2,1)-cmplx(-sqrt(2d0)/sqrt(5d0),0d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if
     if (abs(Z(2,2)-cmplx(1d0/sqrt(1d1),sqrt(2d0)/2d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if     
  else
     if (abs(Z(1,1)-cmplx(1d0/sqrt(1d1),-sqrt(2d0)/2d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if
     if (abs(Z(1,2)-cmplx(sqrt(2d0)/sqrt(5d0),0d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if     
     if (abs(Z(2,1)-cmplx(-sqrt(2d0)/sqrt(5d0),0d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if
     if (abs(Z(2,2)-cmplx(1d0/sqrt(1d1),sqrt(2d0)/2d0,kind=8))>tol) then
        call u_test_failed(__LINE__)
     end if     
  end if


  !check 3) 
  H(1,1) = cmplx(1d0,1d0,kind=8)
  H(1,2) = cmplx(2d0,0d0,kind=8)
  H(2,1) = cmplx(0d0,1d0,kind=8)
  H(2,2) = cmplx(3d0,0d0,kind=8)

   call z_2x2array_eig(H,E,Z,INFO)
  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
 
  ! check results
  ! check E
  if ((abs(E(1)-cmplx(1d0,0.0d0,kind=8))>tol)&
       &.AND.(abs(E(2)-cmplx(1d0,0.0d0,kind=8))>tol)) then
    call u_test_failed(__LINE__)
  end if
  if ((abs(E(1)-cmplx(3d0,1d0,kind=8))>tol)&
       &.AND.(abs(E(2)-cmplx(3d0,1d0,kind=8))>tol)) then
    call u_test_failed(__LINE__)
  end if

  !check Z
  if (abs(abs(Z(1,1))-cmplx(sqrt(2d0)/2d0,0d0,kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(abs(Z(1,2))-cmplx(sqrt(2d0)/2d0,0d0,kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(abs(Z(2,1))-cmplx(sqrt(2d0)/2d0,0d0,kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(abs(Z(2,2))-cmplx(sqrt(2d0)/2d0,0d0,kind=8))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_2x2array_eig
