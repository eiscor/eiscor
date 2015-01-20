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
! 8)
!    a = NAN, b = 1, c = 1 should give Q(1) = Q(2) = Q(3) = nrm = NAN
!
! 9)
!    a = INF, b = 1, c = 1 should give Q(1) = 1, Q(2) = Q(3) = 0 and nrm = INF
!    b = INF, a = 1, c = 1 should give Q(2) = 1, Q(1) = Q(3) = 0 and nrm = INF
!    c = INF, a = 1, b = 1 should give Q(3) = 1, Q(1) = Q(2) = 0 and nrm = INF
!    a = b = c = INF should give Q(1) = Q(2) = Q(3) = 1/sqrt(3) and nrm = INF
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rot3_vec3gen

  implicit none

  ! parameter
  real(8), parameter :: tol = 1d0*EISCOR_DBL_EPS ! accuracy (tolerance)
  real(8), parameter :: alpha = 1d-18 ! small perturbations

  ! compute variables
  real(8) :: rp,rm,nrm,pi = 3.141592653589793239d0
  real(8) :: a, b, c
  real(8) :: Q(3)
  integer :: ii
  
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

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
   

  ! check results
  if (Q(1).NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q(2).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q(3).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.NE.1d0) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, ""
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if
  
  ! set variables
  a = 1d0
  b = 0d0
  c = alpha

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if

  ! set variables
  a = 1d0
  b = -alpha
  c = alpha

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm  
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3) 
     print*, ""
  end if


  ! set variables
  a = 1d0
  b = tol
  c = -tol

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm  
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  a = 1d0
  b = 1d0
  c = 0d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (Q(3).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if

  ! set variables
  a = 1d0
  b = 1d0
  c = alpha

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
    

  ! check results
  if (abs(Q(1)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if

  ! set variables
  a = 1d0
  b = 1d0
  c = -tol

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  a = 1d0
  b = 0d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (Q(2).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if

  ! set variables
  a = 1d0
  b = alpha
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
    

  ! check results
  if (abs(Q(1)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if

  ! set variables
  a = 1d0
  b = -tol
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  a = 1d0
  b = 1d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1)-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 5)

  ! set variables
  a = 0d0
  b = 1d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (Q(1).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if

  ! set variables
  a = alpha
  b = 1d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
    

  ! check results
  if (abs(Q(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if

  ! set variables
  a = -tol
  b = 1d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 6)

  ! set variables
  a = 0d0
  b = 1d0
  c = 0d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
   
  if (Q(1).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q(2).NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q(3).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.NE.1d0) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if
  
  ! set variables
  a = alpha
  b = 1d0
  c = 0d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if

  ! set variables
  a = -alpha
  b = 1d0
  c = alpha

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm  
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3) 
     print*, ""
  end if


  ! set variables
  a = tol
  b = 1d0
  c = -tol

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm  
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 7)

  ! set variables
  a = 0d0
  b = 0d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
   
  if (Q(1).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q(2).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q(3).NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.NE.1d0) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if
  
  ! set variables
  a = alpha
  b = 0d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if

  ! set variables
  a = -alpha
  b = alpha
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm  
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3) 
     print*, ""
  end if


  ! set variables
  a = tol
  b = -tol
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  

  ! check results
  if (abs(Q(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm  
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 8)

  ! set variables
  a = 0d0
  a = a/a
  b = 1d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
   
  if (Q(1).EQ.Q(1)) then
    call u_test_failed(__LINE__)
  end if
  if (Q(2).EQ.Q(2)) then
    call u_test_failed(__LINE__)
  end if
  if (Q(3).EQ.Q(3)) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.EQ.nrm) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 9)

  ! set variables
  a = EISCOR_DBL_INF
  a = a*10d0
  b = 1d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  
  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if
   
  if (Q(1).NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q(2).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q(3).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm <= EISCOR_DBL_INF) then
    call u_test_failed(__LINE__)
  end if
  
  ! set variables
  b = EISCOR_DBL_INF
  b = b*10d0
  a = 1d0
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  
  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if
   
  if (Q(1).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q(2).NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q(3).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm <= EISCOR_DBL_INF) then
    call u_test_failed(__LINE__)
  end if
  
  ! set variables
  c = EISCOR_DBL_INF
  c = c*10d0
  a = 1d0
  b = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  
  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if
   
  if (Q(1).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q(2).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (Q(3).NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm <= EISCOR_DBL_INF) then
    call u_test_failed(__LINE__)
  end if
  
  ! set variables
  c = EISCOR_DBL_INF
  c = c*10d0
  a = c
  b = c
  c = 1d0

  call z_rot3_vec3gen(a,b,c,Q(1),Q(2),Q(3),nrm)
  
  if (VERBOSE) then
     print*, Q, nrm
     print*, a, nrm*Q(1), a - nrm*Q(1)
     print*, b, nrm*Q(2), b - nrm*Q(2)
     print*, c, nrm*Q(3), c - nrm*Q(3)
     print*, ""
  end if
   
  if (Q(1).NE.1d0/sqrt(3d0)) then
    call u_test_failed(__LINE__)
  end if
  if (Q(2).NE.1d0/sqrt(2d0)) then
    call u_test_failed(__LINE__)
  end if
  if (Q(3).NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm <= EISCOR_DBL_INF) then
    call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rot3_vec3gen
