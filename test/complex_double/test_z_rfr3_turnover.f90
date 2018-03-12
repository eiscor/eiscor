#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_rfr3turnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks some limiting cases in the root free turnover.
!
! 1) The case when w = u.
! 
! 2) The case when w + u is strictly real and close to 0.
!
! 3) The case when w + u is complex and close to 0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rfr3_turnover

  implicit none

  ! compute variables
  real(8), parameter :: eps = (EISCOR_DBL_EPS)
  real(8) :: vv, cc, ss, nrm
  complex(8) :: u, w, rho
  complex(8) :: uold

  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)




  !!!!!!!!!!!!!!!!!!!!
  ! check 1)

  ! create u,vv close to 1
  u = cmplx(1d0,0d0,kind=8)
  vv = 0d0
 
  ! set w to -1
  w = cmplx(-1d0,0d0,kind=8)

  ! set cc, ss and rho
  cc = 5d-1
  ss = 5d-1
  rho = cmplx(1d0,0d0,kind=8)

  ! save u
  uold = u

  ! perform turnover
  call z_rfr3_turnover(w,cc,ss,u,vv,rho)

  ! check error
  nrm = abs(-conjg(rho)*cc*w + rho*ss*u - uold)
  if (nrm > 10d0*eps) then
     call u_test_failed(__LINE__)
  end if



  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! create u,vv close to 1
  u = cmplx(1d0-eps,0d0,kind=8)
  vv = eps*(2d0-eps)
 
  ! set w to -1
  w = cmplx(-1d0,0d0,kind=8)

  ! set cc, ss and rho
  cc = 5d-1
  ss = 5d-1
  rho = cmplx(1d0,0d0,kind=8)

  ! save u
  uold = u

  ! perform turnover
  call z_rfr3_turnover(w,cc,ss,u,vv,rho)

  ! check error
  nrm = abs(-conjg(rho)*cc*w + rho*ss*u - uold)
  if (nrm > 10d0*eps) then
     call u_test_failed(__LINE__)
  end if



  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! create u,vv close to 1
  u = cmplx(1d0-eps,0d0,kind=8)
  vv = eps*(2d0-eps)
 
  ! set w close to -1
  call d_rot2_vec2gen(-1d0,2d0*eps,cc,ss,nrm)
  w = cmplx(cc,ss,kind=8)

  ! set cc, ss and rho
  cc = 5d-1
  ss = 5d-1
  rho = cmplx(1d0,0d0,kind=8)

  ! save u
  uold = u

  ! perform turnover
  call z_rfr3_turnover(w,cc,ss,u,vv,rho)

  ! check error
  nrm = abs(-conjg(rho)*cc*w + rho*ss*u - uold)
  if (nrm > 10d0*eps) then
     call u_test_failed(__LINE__)
  end if



  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rfr3_turnover
