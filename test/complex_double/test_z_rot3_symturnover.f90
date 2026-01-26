#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_rot3_symturnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks some limiting cases in the symmetric turnover.
!
! 1) Non-trivial, no deflation
! 
! 2) Diagonal example
! 
! 3) Non-trivial with deflation
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rot3_symturnover

  implicit none

  ! compute variables
  real(8), parameter :: eps = (EISCOR_DBL_EPS)
  real(8) :: nrm, v, s
  complex(8) :: g, c, u, r
  real(8) :: v_in, s_in
  complex(8) :: g_in, c_in, u_in
  real(8) :: v_out, s_out
  complex(8) :: g_out, c_out, u_out

  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)




  !!!!!!!!!!!!!!!!!!!!
  ! check 1)

  ! known inputs
  g_in = cmplx(7.071067811865476d-01,7.071067811865476d-01,kind=8)
  c_in = cmplx(5.773502691896257d-01,5.773502691896257d-01,kind=8)
  s_in =       5.773502691896257d-01
  u_in = cmplx(5.773502691896257d-01,5.773502691896257d-01,kind=8)
  v_in =       5.773502691896257d-01
  r    = cmplx(8.944271909999159d-01,4.472135954999579d-01,kind=8)

  ! known outputs
  g_out = cmplx(-9.486832980505138d-01,-3.162277660168379d-01,kind=8)
  c_out = cmplx(-7.410692737960474d-01,-4.690917827740838d-01,kind=8)
  s_out =        4.803844614152614d-01
  u_out = cmplx( 4.690174004250531d-01,-5.463892354512888d-01,kind=8)
  v_out =        6.938886664887109d-01

  ! store inputs for comparison
  g = g_in
  c = c_in
  s = s_in
  u = u_in
  v = v_in

  ! perform symmetric turnover
  call z_rot3_symturnover(g,c,s,u,v,r)

  ! check error
  nrm = abs(g-g_out) + abs(c-c_out) + abs(s-s_out) + abs(u-u_out) + abs(v-v_out)
  if (nrm > 10d0*eps) then
     call u_test_failed(__LINE__)
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  ! known inputs
  g_in = cmplx(6.000000000000000d-01,-8.000000000000000d-01,kind=8)
  c_in = cmplx(0.000000000000000d+00,+0.000000000000000d+00,kind=8)
  s_in =       1.000000000000000d+00
  u_in = cmplx(4.472135954999579d-01,-8.944271909999159d-01,kind=8)
  v_in =       0.000000000000000d+00
  r    = cmplx(8.944271909999159d-01,+4.472135954999579d-01,kind=8)

  ! known outputs
  g_out = cmplx(-1.788854381999832d-01,+9.838699100999074d-01,kind=8)
  c_out = cmplx( 1.000000000000000d+00,+0.000000000000000d+00,kind=8)
  s_out =        0.000000000000000d+00
  u_out = cmplx( 0.000000000000000d+00,-1.000000000000000d+00,kind=8)
  v_out =        0.000000000000000d+00

  ! store inputs for comparison
  g = g_in
  c = c_in
  s = s_in
  u = u_in
  v = v_in

  ! perform symmetric turnover
  call z_rot3_symturnover(g,c,s,u,v,r)

  ! check error
  nrm = abs(g-g_out) + abs(c-c_out) + abs(s-s_out) + abs(u-u_out) + abs(v-v_out)
  if (nrm > 10d0*eps) then
     call u_test_failed(__LINE__)
  end if


  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  ! known inputs
  g_in = cmplx(1.000000000000000d+00,+0.000000000000000d+00,kind=8)
  c_in = cmplx(0.000000000000000d+00,-7.071067811865476d-01,kind=8)
  s_in =       7.071067811865476d-01
  u_in = cmplx(1.000000000000000d+00,+0.000000000000000d+00,kind=8)
  v_in =       0.000000000000000d+00
  r    = cmplx(0.000000000000000d+00,+1.000000000000000d+00,kind=8)

  ! known outputs
  g_out = cmplx(0.000000000000000d+00,+1.000000000000000d+00,kind=8)
  c_out = cmplx(1.000000000000000d+00,+0.000000000000000d+00,kind=8)
  s_out =       0.000000000000000d+00
  u_out = cmplx(0.000000000000000d+00,-1.000000000000000d+00,kind=8)
  v_out =       0.000000000000000d+00

  ! store inputs for comparison
  g = g_in
  c = c_in
  s = s_in
  u = u_in
  v = v_in

  ! perform symmetric turnover
  call z_rot3_symturnover(g,c,s,u,v,r)

  ! check error
  nrm = abs(g-g_out) + abs(c-c_out) + abs(s-s_out) + abs(u-u_out) + abs(v-v_out)
  if (nrm > 10d0*eps) then
     call u_test_failed(__LINE__)
  end if



  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rot3_symturnover
