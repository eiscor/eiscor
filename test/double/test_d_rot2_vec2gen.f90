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
! 4) [ 0 ] = [ 1 ] [ 0 ] 
!    [ 0 ]   [ 0 ] 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_rot2_vec2gen

  implicit none

  ! parameter
  real(8) :: alpha = 1d-18 ! small perturbations
  
  ! compute variables
  real(8) :: nrm, nul=0d0
  real(8) :: a, b
  real(8) :: c, s
    
  ! timing variables
  integer:: c_start, c_stop, c_rate, ii


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

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (VERBOSE) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm
  end if   
  ! check results
  if (c.NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (s.NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.NE.1d0) then
    call u_test_failed(__LINE__)
  end if  

  ! set variables
  a = 1d0
  b = -alpha

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (VERBOSE) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm
  end if     
  ! check results
  if (abs(c-1d0).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = 1d0
  b = EISCOR_DBL_EPS

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (VERBOSE) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm
  end if     
  ! check results
  if (abs(c-1d0).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! set variables
  a = 1d0
  b = 1d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (VERBOSE) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm
  end if   
  ! check results
  if (abs(c-sqrt(2d0)/2d0).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s-sqrt(2d0)/2d0).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0)).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if

  
  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  a = 0d0
  b = 1d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (VERBOSE) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm
  end if   
  ! check results
  if (c.NE.0d0) then
    call u_test_failed(__LINE__)
  end if
  if (s.NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.NE.1d0) then
    call u_test_failed(__LINE__)
  end if
  
  ! set variables
  a = alpha
  b = 1d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (VERBOSE) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm
  end if     
  ! check results
  if (abs(c).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s-1d0).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = -alpha
  b = 1d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (VERBOSE) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm
  end if     
  ! check results
  if (abs(c).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s-1d0).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = EISCOR_DBL_EPS
  b = 1d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (VERBOSE) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm
  end if    
  ! check results
  if (abs(c).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s-1d0).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 4)
  ! set variables
  a = 0d0
  b = 0d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (VERBOSE) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm
  end if     
  ! check results
  if (abs(c-1d0).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm).GT.EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if

 
  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_d_rot2_vec2gen
