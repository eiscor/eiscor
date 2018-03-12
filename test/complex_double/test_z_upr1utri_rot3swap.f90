#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1utri_rot3swap
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1utri_rot3swap. 
!
! 1)  semi-random D, C and B
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1utri_rot3swap

  implicit none

  ! parameter
  real(8), parameter :: tol = 4d0*EISCOR_DBL_EPS ! accuracy (tolerance)

  ! compute variables
  real(8) :: nrm, D(4), C(6), B(6), G(3)
  complex(8) :: T(2,2), Told(2,2)
  complex(8) :: A(2,2), Aold(2,2)

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
  call d_rot2_vec2gen(1d0,-3d0,D(1),D(2),nrm)
  call d_rot2_vec2gen(-17d0,3d0,D(3),D(4),nrm)

  call z_rot3_vec3gen(1d0,2d0,3d0,C(1),C(2),C(3),nrm)
  call z_rot3_vec3gen(9d0,2d0,5d0,C(4),C(5),C(6),nrm)
  
  call z_rot3_vec3gen(11d0,-2d0,3d0,B(1),B(2),B(3),nrm)
  call z_rot3_vec3gen(-3d0,47d0,5d0,B(4),B(5),B(6),nrm)
  
  call z_rot3_vec3gen(11d0,19d0,-7d0,G(1),G(2),G(3),nrm)
  
  ! fill Aold
  Aold(1,1) = cmplx(G(1),G(2),kind=8)
  Aold(2,1) = cmplx(G(3),0d0,kind=8)
  Aold(1,2) = -Aold(2,1)
  Aold(2,2) = conjg(Aold(1,1))
  
  ! fill Told
  call z_upr1utri_decompress(.FALSE.,2,D,C,B,Told)

  ! call 
  call z_upr1utri_rot3swap(.FALSE.,D,C,B,G)
  
  ! fill A
  A(1,1) = cmplx(G(1),G(2),kind=8)
  A(2,1) = cmplx(G(3),0d0,kind=8)
  A(1,2) = -A(2,1)
  A(2,2) = conjg(A(1,1))
  
  ! fill T
  call z_upr1utri_decompress(.FALSE.,2,D,C,B,T)

  ! check results
  A = abs(matmul(Told,Aold)-matmul(A,T))
  if (maxval(dble(A)) > tol) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_upr1utri_rot3swap
