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
! 0) checks every exceptional case in comment block
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
  real(8), parameter :: tol = 10d0*EISCOR_DBL_EPS ! accuracy (tolerance)
  real(8), parameter :: alpha = 1d-18 ! small perturbations
  real(8), parameter :: eps = 1d0*EISCOR_DBL_EPS ! accuracy (tolerance)
  real(8), parameter :: pi = EISCOR_DBL_PI

  ! compute variables
  real(8) :: rp,rm,NRM
  real(8) :: AR, AI, BR, BI
  real(8) :: CR, CI, S
  real(8) :: inf, nan
  integer :: ii
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)
  
  ! set inf
  inf = EISCOR_DBL_INF
  inf = 10d0*inf
  
  ! set nan
  nan = 0d0
  nan = 0d0/nan
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 0)
    ! all zeros
    AR = 0; AI = 0; BR = 0; BI = 0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.NE.1d0).OR.(CI.NE.0d0).OR.(S.NE.0d0).OR.(NRM.NE.0d0)) then
      call u_test_failed(__LINE__)
    end if
    
    ! one INF
    AR = inf; AI = 1d0; BR = 1d0; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)

    if ((CR.NE.1d0).OR.(CI.NE.0d0).OR.(S.NE.0d0).OR.(NRM.NE.inf)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = inf; BR = 1d0; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.NE.0d0).OR.(CI.NE.1d0).OR.(S.NE.0d0).OR.(NRM.NE.inf)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = 1d0; BR = inf; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.NE.0d0).OR.(CI.NE.0d0).OR.(S.NE.1d0).OR.(NRM.NE.inf)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = 1d0; BR = 1d0; BI = inf
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.NE.0d0).OR.(CI.NE.0d0).OR.(S.NE.1d0).OR.(NRM.NE.inf)) then
      call u_test_failed(__LINE__)
    end if
    
    ! one -INF
    AR = -inf; AI = 1d0; BR = 1d0; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.NE.-1d0).OR.(CI.NE.0d0).OR.(S.NE.0d0).OR.(NRM.NE.inf)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = -inf; BR = 1d0; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.NE.0d0).OR.(CI.NE.-1d0).OR.(S.NE.0d0).OR.(NRM.NE.inf)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = 1d0; BR = -inf; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.NE.0d0).OR.(CI.NE.0d0).OR.(S.NE.1d0).OR.(NRM.NE.inf)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = 1d0; BR = 1d0; BI = -inf
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.NE.0d0).OR.(CI.NE.0d0).OR.(S.NE.1d0).OR.(NRM.NE.inf)) then
      call u_test_failed(__LINE__)
    end if
    
    ! two INFs
    AR = inf; AI = inf; BR = 1d0; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = inf; BR = inf; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = 1d0; BR = inf; BI = inf
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = inf; AI = 1d0; BR = inf; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = inf; BR = 1d0; BI = inf
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = inf; AI = 1d0; BR = 1d0; BI = inf
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    ! three INFs
    AR = inf; AI = inf; BR = inf; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = inf; BR = inf; BI = inf
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = inf; AI = 1d0; BR = inf; BI = inf
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = inf; AI = inf; BR = 1d0; BI = inf
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if 
    
    ! four INFs
    AR = inf; AI = inf; BR = inf; BI = inf
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
     ! one NAN
    AR = nan; AI = 1d0; BR = 1d0; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = nan; BR = 1d0; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = 1d0; BR = nan; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = 1d0; BR = 1d0; BI = nan
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    ! two NANs
    AR = nan; AI = nan; BR = 1d0; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = nan; BR = nan; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = 1d0; BR = nan; BI = nan
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = nan; AI = 1d0; BR = nan; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = nan; BR = 1d0; BI = nan
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = nan; AI = 1d0; BR = 1d0; BI = nan
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    ! three NANs
    AR = nan; AI = nan; BR = nan; BI = 1d0
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = 1d0; AI = nan; BR = nan; BI = nan
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = nan; AI = 1d0; BR = nan; BI = nan
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
    
    AR = nan; AI = nan; BR = 1d0; BI = nan
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if 
    
    ! four NANs
    AR = nan; AI = nan; BR = nan; BI = nan
    call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
    if ((CR.EQ.CR).OR.(CI.EQ.CI).OR.(S.EQ.S).OR.(NRM.EQ.NRM)) then
      call u_test_failed(__LINE__)
    end if
  

  !!!!!!!!!!!!!!!!!!!!
  ! check 1)

  ! set variables
  AR = 1d0
  AI = 0d0
  BR = 0d0
  BI = 0d0
  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)

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

  ! set variables
  AR = 1d0
  AI = 0d0
  BR = alpha
  BI = alpha

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  
  ! check results
  if (abs(CR-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI+sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  AR = 1d0
  AI = -alpha
  BR = alpha
  BI = -alpha
  
  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  
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
  if (abs(NRM-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  AR = 1d0
  AI = eps
  BR = -eps
  BI = eps
  
  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR+sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI+sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S)>tol) then
  print*, tol, S
  
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  AR = 1d0
  AI = 1d0
  BR = 0d0
  BI = 0d0

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)


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

  ! set variables
  AR = 1d0
  AI = 1d0
  BR = alpha
  BI = -alpha

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
   
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
  if (abs(NRM-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  AR = 1d0
  AI = 1d0
  BR = -eps
  BI = eps

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  

  ! check results
  if (abs(CR)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI+1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  AR = 1d0
  AI = 0d0
  BR = 1d0
  BI = 0d0

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  

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

  ! set variables
  AR = 1d0
  AI = alpha
  BR = 1d0
  BI = alpha

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)


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

  ! set variables
  AR = 1d0
  AI = -eps
  BR = 1d0
  BI = eps

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  

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

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  ! set variables
  AR = 1d0
  AI = 1d0
  BR = 1d0
  BI = 0d0

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  

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

  !!!!!!!!!!!!!!!!!!!!
  ! check 5)

  ! set variables
  AR = 0d0
  AI = 1d0
  BR = 1d0
  BI = 0d0

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  

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
  
  ! set variables
  AR = alpha
  AI = 1d0
  BR = 1d0
  BI = -alpha

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
 
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

  ! set variables
  AR = -eps
  AI = 1d0
  BR = 1d0
  BI = -eps
  
  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  

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

  !!!!!!!!!!!!!!!!!!!!
  ! check 6)

  ! set variables
  AR = 0d0
  AI = 1d0
  BR = 0d0
  BI = 0d0

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
   
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

  ! set variables
  AR = alpha
  AI = 1d0
  BR = 0d0
  BI = alpha

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  

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

  ! set variables
  AR = -alpha
  AI = 1d0
  BR = alpha
  BI = alpha

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  

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
  if (abs(NRM-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  AR = -eps
  AI = 1d0
  BR = eps
  BI = 0d0
  
  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)

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

  !!!!!!!!!!!!!!!!!!!!
  ! check 7)

  ! set variables
  AR = 0d0
  AI = 0d0
  BR = 1d0
  BI = 0d0

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
   
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

  ! set variables
  AR = alpha
  AI = 0d0
  BR = 1d0
  BI = alpha
  
  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  

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

  ! set variables
  AR = -alpha
  AI = alpha
  BR = 1d0
  BI = alpha

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  

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

  ! set variables
  AR = eps
  AI = -eps
  BR = 1d0
  BI = eps

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)
  

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

  !!!!!!!!!!!!!!!!!!!!
  ! check 8)
  ! set variables
  AR = 3d0
  AI = 2d0
  BR = 1d0
  BI = -4d0

  call z_rot3_vec4gen(AR,AI,BR,BI,CR,CI,S,NRM)

  ! check results
  if (abs(CR+sqrt(5d0/17d0/6d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(CI-sqrt(17d0/30d0)*2+sqrt(80d0/102d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(S-sqrt(17d0/30d0))>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(NRM-sqrt(30d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rot3_vec4gen
