#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_rot3_vec3gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_rot3_vec3gen (generating rotations). 
! The following tests are run:
!
! 1)          
!    [ 1 + 0i ] = [ 1 ] [ 1 ]
!    [   0    ]   [ 0 ]
!                
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
! 2)              
!    [ 1 + i ] = [ sqrt(2)/2 + isqrt(2)/2  ] [ sqrt(2) ]
!    [   0   ]   [           0             ]
!
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
! 3)              
!    [ 1 + 0i ] = [ sqrt(2)/2 ] [ sqrt(2) ]
!    [   1    ]   [ sqrt(2)/2 ]
!
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
! 4)              
!    [ 1 + i ] = [ 1/sqrt(3) + i 1/sqrt(3)] [ sqrt(3) ]
!    [   1   ]   [ 1/sqrt(3)              ]
!
! 5)              
!    [ 0 + i ] = [ isqrt(2)/2 ] [ sqrt(2) ]
!    [   1   ]   [  sqrt(2)/2 ]
!
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
! 6)              
!    [ 0 + i ] = [ i ] [ 1 ]
!    [   0   ]   [ 0 ]
!
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
! 7)              
!    [ 0 + 0i ] = [ 0 ] [ 1 ]
!    [   1    ]   [ 1 ]
!
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
! In DEBUG mode additionally checks with incorrect input are run. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rot3_vec3gen

  implicit none

  ! parameter
  real(8) :: tol = 1d0*epsilon(1d0) ! accuracy (tolerance)
  real(8) :: alpha = 1d-18 ! small perturbations

  ! compute variables
  real(8) :: rp,rm,nrm,pi = 3.141592653589793239d0
  real(8) :: a, b, c
  real(8) :: Q1(3), nul = 0d0
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
  c = 0d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
   
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
  if (Q1(3).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.NE.1d0) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, ""
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if
  
  ! set variables
  a = 1d0
  b = 0d0
  c = alpha

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
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
  if (abs(Q1(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  ! set variables
  a = 1d0
  b = -alpha
  c = alpha

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
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
  if (abs(Q1(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm  
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3) 
     print*, ""
  end if


  ! set variables
  a = 1d0
  b = tol
  c = -tol

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
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
  if (abs(Q1(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm  
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  a = 1d0
  b = 1d0
  c = 0d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
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
  if (Q1(3).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  ! set variables
  a = 1d0
  b = 1d0
  c = alpha

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
    
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
  if (abs(Q1(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  ! set variables
  a = 1d0
  b = 1d0
  c = -tol

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
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
  if (abs(Q1(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  a = 1d0
  b = 0d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (Q1(2).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(3)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  ! set variables
  a = 1d0
  b = alpha
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
    
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(3)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  ! set variables
  a = 1d0
  b = -tol
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(3)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  a = 1d0
  b = 1d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1)-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2)-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(3)-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 5)

  ! set variables
  a = 0d0
  b = 1d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (Q1(1).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(3)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  ! set variables
  a = alpha
  b = 1d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
    
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(3)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  ! set variables
  a = -tol
  b = 1d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(3)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 6)

  ! set variables
  a = 0d0
  b = 1d0
  c = 0d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
   
  if (Q1(1).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q1(2).NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q1(3).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.NE.1d0) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if
  
  ! set variables
  a = alpha
  b = 1d0
  c = 0d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
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
  if (abs(Q1(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  ! set variables
  a = -alpha
  b = 1d0
  c = alpha

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
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
  if (abs(Q1(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm  
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3) 
     print*, ""
  end if


  ! set variables
  a = tol
  b = 1d0
  c = -tol

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
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
  if (abs(Q1(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm  
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 7)

  ! set variables
  a = 0d0
  b = 0d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
   
  if (Q1(1).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q1(2).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q1(3).NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.NE.1d0) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if
  
  ! set variables
  a = alpha
  b = 0d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(3)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  ! set variables
  a = -alpha
  b = alpha
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(3)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm  
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3) 
     print*, ""
  end if


  ! set variables
  a = tol
  b = -tol
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q1(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q1(3)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q1, nrm  
     print*, a, nrm*Q1(1), a - nrm*Q1(1)
     print*, b, nrm*Q1(2), b - nrm*Q1(2)
     print*, c, nrm*Q1(3), c - nrm*Q1(3)
     print*, ""
  end if

  if (DEBUG) then

     a = a/0d0
     call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)     
     ! check INFO
     if (INFO.NE.-1) then
        call u_test_failed(__LINE__)
     end if

     a = 2d0
     b = a/0d0
     call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)  
     ! check INFO
     if (INFO.NE.-2) then
        call u_test_failed(__LINE__)
     end if

     b = 3d0
     c = a/0d0
     call z_rot3_vec3gen(a,b,c,Q1(1),Q1(2),Q1(3),nrm,info)     
     ! check INFO
     if (INFO.NE.-3) then
        call u_test_failed(__LINE__)
     end if     

  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rot3_vec3gen
