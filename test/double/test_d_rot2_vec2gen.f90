#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_rot2_vec2gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine d_rot2_vec2gen (generating rotations). 
! The following tests are run:
!
! 1)          
!    [ 1 ] = [ 1 ] [ 1 ]
!    [ 0 ]   [ 0 ]
!                
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
! 2)              
!    [ 1 ] = [ sqrt(2)/2 ] [ sqrt(2) ]
!    [ 1 ]   [ sqrt(2)/2 ]
!
! 3)              
!    [ 0 ] = [ 0 ] [ 1 ]
!    [ 1 ]   [ 1 ]
!
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_rot2_vec2gen

  implicit none

  ! parameter
  real(8) :: tol = 1d0*epsilon(1d0) ! accuracy (tolerance)
  real(8) :: alpha = 1d-18 ! small perturbations

  ! compute variables
  real(8) :: nrm
  real(8) :: a, b
  real(8) :: Q1(2)
  integer :: info, ii
  
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
  a = 1d0
  b = 0d0

  call d_rot2_vec2gen(a,b,Q1(1),Q1(2),nrm,info)
   
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (Q1(1).NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q1(2).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.NE.1d0) then
    call u_test_failed(__LINE__)
  end if

  if (DEBUG) then
     print*, ""
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, ""
  end if
  

  ! set variables
  a = 1d0
  b = -alpha

  call d_rot2_vec2gen(a,b,Q1(1),Q1(2),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (DEBUG) then
     print*, Q1, nrm  
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, ""
  end if


  ! set variables
  a = 1d0
  b = tol

  call d_rot2_vec2gen(a,b,Q1(1),Q1(2),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (DEBUG) then
     print*, Q1, nrm  
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, ""
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  a = 1d0
  b = 1d0

  call d_rot2_vec2gen(a,b,Q1(1),Q1(2),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  if (DEBUG) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, ""
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  a = 0d0
  b = 1d0

  call d_rot2_vec2gen(a,b,Q1(1),Q1(2),nrm,info)
   
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (Q1(1).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q1(2).NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.NE.1d0) then
    call u_test_failed(__LINE__)
  end if

  if (DEBUG) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, ""
  end if
  
  ! set variables
  a = alpha
  b = 1d0

  call d_rot2_vec2gen(a,b,Q1(1),Q1(2),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (DEBUG) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, ""
  end if

  ! set variables
  a = -alpha
  b = 1d0

  call d_rot2_vec2gen(a,b,Q1(1),Q1(2),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (DEBUG) then
     print*, Q1, nrm  
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, ""
  end if


  ! set variables
  a = tol
  b = 1d0

  call d_rot2_vec2gen(a,b,Q1(1),Q1(2),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (DEBUG) then
     print*, Q1, nrm  
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, ""
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_d_rot2_vec2gen
