#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rot3vec_vec2block

  implicit none

  ! compute variables
  integer, parameter :: N = 2
  integer :: ii
  real(8) :: ar, ai, b, nrm
  real(8) :: Q(3*N*N)
  complex(8) :: Z(2*N,2*N)

  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start) 

  ! print banner
  call u_test_banner(__FILE__)

  ! fill Q
  do ii = 1,N*N
    call random_number(ar)
    call random_number(ai)
    call random_number(b)
    call z_rot3_vec3gen(ar,ai,b,Q(3*ii-2),Q(3*ii-1),Q(3*ii),nrm)
  end do

  ! fill Z
  call z_rot3vec_vec2block(.TRUE.,N,Q,2*N,Z)

print*,""
print*," Z"
do ii = 1,2*N
print*,Z(ii,:)
end do
print*,""

Z = matmul(conjg(transpose(Z)),Z)
print*,""
print*," Z^H*Z"
do ii = 1,2*N
print*,Z(ii,:)
end do
print*,""

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rot3vec_vec2block
