#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_rot3_vec4gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_rot3_vec4gen (generating rotations). 
! The following tests are run:
!
! 1)          
!    [ 1 + 0i ] = [ 1 ] [ 1 ]
!    [ 0 + 0i ]   [ 0 ]
!                
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
! 2)              
!    [ 1 + i  ] = [ sqrt(2)/2 + isqrt(2)/2  ] [ sqrt(2) ]
!    [ 0 + 0i ]   [           0             ]
!
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
! 3)              
!    [ 1 + 0i ] = [ sqrt(2)/2 ] [ sqrt(2) ]
!    [ 1 + 0i ]   [ sqrt(2)/2 ]
!
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
! 4)              
!    [ 1 + i  ] = [ 1/sqrt(3) + i 1/sqrt(3)] [ sqrt(3) ]
!    [ 1 + 0i ]   [ 1/sqrt(3)              ]
!
! 5)              
!    [ 0 + i  ] = [ isqrt(2)/2 ] [ sqrt(2) ]
!    [ 1 + 0i ]   [  sqrt(2)/2 ]
!
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
! 6)              
!    [ 0 + i  ] = [ i ] [ 1 ]
!    [ 0 + 0i ]   [ 0 ]
!
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
! 7)              
!    [ 0 + 0i ] = [ 0 ] [ 1 ]
!    [ 1 + 0i ]   [ 1 ]
!
!    and replace 0 by (+-) 1e-18 and (+-) eps 
!
! 8)              
!    [ 3 + 2i ] = [ 0 ] [ 1 ]
!    [ 1 - 4i ]   [ 1 ]
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rot3_vec4gen

  implicit none

  ! parameter
  real(8) :: tol = 2d0*epsilon(1d0) ! accuracy (tolerance)
  real(8) :: alpha = 1d-18 ! small perturbations
  real(8) :: eps = 1d0*epsilon(1d0) ! accuracy (tolerance)

  ! compute variables
  real(8) :: rp,rm,nrm,pi = 3.141592653589793239d0
  real(8) :: a, b, c, d
  real(8) :: Q(3)
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
  d = 0d0
  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
   
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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

  ! set variables
  a = 1d0
  b = 0d0
  c = alpha
  d = alpha

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q(1)-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)+sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = 1d0
  b = -alpha
  c = alpha
  d = -alpha
  
  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = 1d0
  b = eps
  c = -eps
  d = eps
  
  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q(1)+sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)+sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
  print*, tol, Q(3)
  
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  a = 1d0
  b = 1d0
  c = 0d0
  d = 0d0

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)

  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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

  ! set variables
  a = 1d0
  b = 1d0
  c = alpha
  d = -alpha

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
      
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = 1d0
  b = 1d0
  c = -eps
  d = eps

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q(1))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)+1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  a = 1d0
  b = 0d0
  c = 1d0
  d = 0d0

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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

  ! set variables
  a = 1d0
  b = alpha
  c = 1d0
  d = alpha

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)

  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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

  ! set variables
  a = 1d0
  b = -eps
  c = 1d0
  d = eps

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  a = 1d0
  b = 1d0
  c = 1d0
  d = 0d0

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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

  !!!!!!!!!!!!!!!!!!!!
  ! check 5)

  ! set variables
  a = 0d0
  b = 1d0
  c = 1d0
  d = 0d0

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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
  
  ! set variables
  a = alpha
  b = 1d0
  c = 1d0
  d = -alpha

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
    
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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

  ! set variables
  a = -eps
  b = 1d0
  c = 1d0
  d = -eps
  
  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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

  !!!!!!!!!!!!!!!!!!!!
  ! check 6)

  ! set variables
  a = 0d0
  b = 1d0
  c = 0d0
  d = 0d0

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
   
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

  ! set variables
  a = alpha
  b = 1d0
  c = 0d0
  d = alpha

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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

  ! set variables
  a = -alpha
  b = 1d0
  c = alpha
  d = alpha

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = eps
  b = 1d0
  c = -eps
  d = 0d0
  
  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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

  !!!!!!!!!!!!!!!!!!!!
  ! check 7)

  ! set variables
  a = 0d0
  b = 0d0
  c = 1d0
  d = 0d0

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
   
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

  ! set variables
  a = alpha
  b = 0d0
  c = 1d0
  d = alpha
  
  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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

  ! set variables
  a = -alpha
  b = alpha
  c = 1d0
  d = alpha

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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

  ! set variables
  a = eps
  b = -eps
  c = 1d0
  d = eps

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
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

  !!!!!!!!!!!!!!!!!!!!
  ! check 8)
  ! set variables
  a = 3d0
  b = 2d0
  c = 1d0
  d = -4d0

  call z_rot3_vec4gen(a,b,c,d,Q(1),Q(2),Q(3),nrm,info)
  
  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  ! check results
  if (abs(Q(1)+sqrt(5d0/17d0/6d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-sqrt(17d0/30d0)*2+sqrt(80d0/102d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-sqrt(17d0/30d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(30d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rot3_vec4gen
