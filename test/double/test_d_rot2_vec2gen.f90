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
  real(8), parameter :: infdef = huge(1d0)
  !integer, parameter :: notests = 1000000000 ! 1 billion
  integer, parameter :: notests = 1000000 ! 1 million
  
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
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
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
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if     
  ! check results
  if (abs(c-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = 1d0
  b = tol

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if     
  ! check results
  if (abs(c-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s)>tol) then
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

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if   
  ! check results
  if (abs(c-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s-sqrt(2d0)/2d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-sqrt(2d0))>tol) then
    call u_test_failed(__LINE__)
  end if

  
  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! set variables
  a = 0d0
  b = 1d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
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
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if     
  ! check results
  if (abs(c)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = -alpha
  b = 1d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if     
  ! check results
  if (abs(c)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = tol
  b = 1d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if    
  ! check results
  if (abs(c)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(nrm-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! Test with INF

  ! set variables
  a = 1d0/nul ! INF
  b = 0d0     ! 0d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (abs(c-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (.NOT.(nrm>infdef)) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = 1d0/nul ! INF
  b = 2d0     ! 2d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (abs(c-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (.NOT.(nrm>infdef)) then
    call u_test_failed(__LINE__)
  end if
 
  ! set variables
  a = 1d0/nul  ! INF
  b = -2d0     ! -2d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (abs(c-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (.NOT.(nrm>infdef)) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = -1d0/nul ! -INF
  b = 0d0      ! 0d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (abs(c+1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (.NOT.(nrm>infdef)) then
    call u_test_failed(__LINE__)
  end if
  
  ! set variables
  a = -1d0/nul ! -INF
  b = 2d0      ! 2d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (abs(c+1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (.NOT.(nrm>infdef)) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = -1d0/nul ! -INF
  b = -2d0     ! -2d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (abs(c+1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(s)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (.NOT.(nrm>infdef)) then
    call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!

  ! set variables
  b = 1d0/nul ! INF
  a = 0d0     ! 0d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (abs(s-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(c)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (.NOT.(nrm>infdef)) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  b = 1d0/nul ! INF
  a = 2d0     ! 2d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (abs(s-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(c)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (.NOT.(nrm>infdef)) then
    call u_test_failed(__LINE__)
  end if
 
  ! set variables
  b = 1d0/nul  ! INF
  a = -2d0     ! -2d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (abs(s-1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(c)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (.NOT.(nrm>infdef)) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  b = -1d0/nul ! -INF
  a = 0d0      ! 0d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (abs(s+1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(c)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (.NOT.(nrm>infdef)) then
    call u_test_failed(__LINE__)
  end if
  
  ! set variables
  b = -1d0/nul ! -INF
  a = 2d0      ! 2d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (abs(s+1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(c)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (.NOT.(nrm>infdef)) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  b = -1d0/nul ! -INF
  a = -2d0     ! -2d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (abs(s+1d0)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (abs(c)>tol) then
    call u_test_failed(__LINE__)
  end if
  if (.NOT.(nrm>infdef)) then
    call u_test_failed(__LINE__)
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! set variables
  a = 1d0/nul ! INF
  b = 1d0/nul ! INF

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (c.EQ.c) then
    call u_test_failed(__LINE__)
  end if
  if (s.EQ.s) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.EQ.nrm) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = -1d0/nul ! -INF
  b = 1d0/nul  ! INF

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (c.EQ.c) then
    call u_test_failed(__LINE__)
  end if
  if (s.EQ.s) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.EQ.nrm) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = 1d0/nul  ! INF
  b = -1d0/nul ! -INF

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (c.EQ.c) then
    call u_test_failed(__LINE__)
  end if
  if (s.EQ.s) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.EQ.nrm) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = -1d0/nul ! -INF
  b = -1d0/nul ! -INF

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (c.EQ.c) then
    call u_test_failed(__LINE__)
  end if
  if (s.EQ.s) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.EQ.nrm) then
    call u_test_failed(__LINE__)
  end if
 

  !!!!!!!!!!!!!!!!!!!!
  ! Test with NAN

  ! set variables
  a = 0d0/nul ! NAN
  b = 1d0     ! 1d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (c.EQ.c) then
    call u_test_failed(__LINE__)
  end if
  if (s.EQ.s) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.EQ.nrm) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = 1d0     ! 1d0
  b = 0d0/nul ! NAN

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (c.EQ.c) then
    call u_test_failed(__LINE__)
  end if
  if (s.EQ.s) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.EQ.nrm) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = 0d0/nul ! NAN
  b = 0d0/nul ! NAN

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (c.EQ.c) then
    call u_test_failed(__LINE__)
  end if
  if (s.EQ.s) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.EQ.nrm) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = 0d0     ! 0d0
  b = 0d0/nul ! NAN

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (c.EQ.c) then
    call u_test_failed(__LINE__)
  end if
  if (s.EQ.s) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.EQ.nrm) then
    call u_test_failed(__LINE__)
  end if

  ! set variables
  a = 0d0/nul ! NAN
  b = 0d0     ! 0d0

  call d_rot2_vec2gen(a,b,c,s,nrm)
  if (DEBUG) then
     print*, "a", a, "b", b, "c", c, "s", s, "nrm", nrm, "infdef", infdef, "nrm>infdef", (nrm>infdef)
  end if
  ! check results
  if (c.EQ.c) then
    call u_test_failed(__LINE__)
  end if
  if (s.EQ.s) then
    call u_test_failed(__LINE__)
  end if
  if (nrm.EQ.nrm) then
    call u_test_failed(__LINE__)
  end if

  do ii=1,notests
     call random_number(a)
     call random_number(b)
     call d_rot2_vec2gen(a,b,c,s,nrm)
  end do
  
  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_d_rot2_vec2gen
