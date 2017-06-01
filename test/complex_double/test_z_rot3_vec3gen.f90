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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rot3_vec3gen

  implicit none

  ! parameter
  real(8), parameter :: tol = 1d0*EISCOR_DBL_EPS ! accuracy (tolerance)
  real(8), parameter :: alpha = 1d-18 ! small perturbations
  real(8), parameter :: pi = EISCOR_DBL_PI

  ! compute variables
  real(8) :: rp, rm
  real(8) :: AR, AI, B, CR, CI, S, NRM
  real(8) :: inf, nan
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
  AR = 1d0
  AI = 0d0
  B = 0d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
   

  ! check results
  if (CR.NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (CI.NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (S.NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (NRM.NE.1d0) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, ""
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if
  
  ! set variables
  AR = 1d0
  AI = 0d0
  B = alpha

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  ! set variables
  AR = 1d0
  AI = -alpha
  B = alpha

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM  
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S 
     print*, ""
  end if


  ! set variables
  AR = 1d0
  AI = tol
  B = -tol

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM  
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  AR = 1d0
  AI = 1d0
  B = 0d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (S.NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  ! set variables
  AR = 1d0
  AI = 1d0
  B = alpha

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
    

  ! check results
  if (abs(CR-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  ! set variables
  AR = 1d0
  AI = 1d0
  B = -tol

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  AR = 1d0
  AI = 0d0
  B = 1d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (CI.NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  ! set variables
  AR = 1d0
  AI = alpha
  B = 1d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
    

  ! check results
  if (abs(CR-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  ! set variables
  AR = 1d0
  AI = -tol
  B = 1d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  AR = 1d0
  AI = 1d0
  B = 1d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S-1d0/sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-sqrt(3d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 5)

  ! set variables
  AR = 0d0
  AI = 1d0
  B = 1d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (CR.NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  ! set variables
  AR = alpha
  AI = 1d0
  B = 1d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
    

  ! check results
  if (abs(CR)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  ! set variables
  AR = -tol
  AI = 1d0
  B = 1d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 6)

  ! set variables
  AR = 0d0
  AI = 1d0
  B = 0d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
   
  if (CR.NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (CI.NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (S.NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (NRM.NE.1d0) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if
  
  ! set variables
  AR = alpha
  AI = 1d0
  B = 0d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  ! set variables
  AR = -alpha
  AI = 1d0
  B = alpha

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM  
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S 
     print*, ""
  end if


  ! set variables
  AR = tol
  AI = 1d0
  B = -tol

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM  
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 7)

  ! set variables
  AR = 0d0
  AI = 0d0
  B = 1d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
   
  if (CR.NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (CI.NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (S.NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (NRM.NE.1d0) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if
  
  ! set variables
  AR = alpha
  AI = 0d0
  B = 1d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  ! set variables
  AR = -alpha
  AI = alpha
  B = 1d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM  
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S 
     print*, ""
  end if


  ! set variables
  AR = tol
  AI = -tol
  B = 1d0

  call z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  if (VERBOSE) then
     print*, CR, CI, S, NRM  
     print*, AR, NRM*CR, AR - NRM*CR
     print*, AI, NRM*CI, AI - NRM*CI
     print*, B, NRM*S, B - NRM*S
     print*, ""
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rot3_vec3gen
