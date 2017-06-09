#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_orthhess_factor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_d_orthhess_factor.
! The following tests are run:
!
! 1)  Test cyclic for N = 2**i, i = 1, ..., 12
!     [ 0 0 ... 0 1]
!     [ 1 0 ... 0 0]
!     [ 0 1 ... 0 0]
!     [ ... ... ...]
!     [ 0 0 ... 1 0]
!
! 2) Q, D given => H d_orthhess_factor has to 
!    reproduce Q and D, N = 17
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_orthhess_factor

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*EISCOR_DBL_EPS ! accuracy (tolerance)
  integer, parameter :: MMAX = 12
  integer, parameter :: NMAX = 2**12
  
  ! compute variables
  integer :: N, M, ii, jj, INFO
  real(8) :: Q(2*(NMAX-1))
  real(8) :: D(NMAX), nrm
  real(8), allocatable :: H(:,:), A(:,:), B(:,:), H2(:,:), Ho(:,:)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 1)
  do ii=1,MMAX
     N = 2**ii
     allocate(H(N,N))
     H = 0d0
     do jj=1,N-1
        H(jj+1,jj) = 1d0
     end do
     H (1,N) = 1d0
     
     call d_orthhess_factor(N,H,Q,D,INFO)
         
     ! check info
     if (INFO.NE.0) then
        call u_test_failed(__LINE__)
     end if
  
     do jj=1,N-1
        ! check Q
        if (abs(Q(2*jj-1))>tol) then
           call u_test_failed(__LINE__)
        end if
        if (abs(Q(2*jj)-1d0)>tol) then
           call u_test_failed(__LINE__)
        end if
        ! check D
        if (abs(D(jj)-1d0)>tol) then
           call u_test_failed(__LINE__)
        end if
     end do
     if (abs(D(N)+1d0)>tol) then
        call u_test_failed(__LINE__)
     end if
     deallocate(H)
  end do


  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  N = 17
  allocate(H(N,N),H2(N,N),Ho(N,N),A(N,N),B(N,N))
  H = 0d0
  do jj=1,N
     H(jj,jj) = 1d0
  end do
  A = H
  
  D = 0d0
  do jj=1,8
     D(2*jj-1) =  1d0
     D(2*jj) = -1d0
  end do
  D(N) = 1d0
  
  call d_rot2_vec2gen(1d0,2d0,Q(1),Q(2),nrm)
  call d_rot2_vec2gen(2d0,1d0,Q(3),Q(4),nrm)
  call d_rot2_vec2gen(3d0,-5d0,Q(5),Q(6),nrm)
  call d_rot2_vec2gen(-1d0,6d0,Q(7),Q(8),nrm)
  call d_rot2_vec2gen(2d0,2d0,Q(9),Q(10),nrm)
  call d_rot2_vec2gen(3d0,-1d0,Q(11),Q(12),nrm)
  call d_rot2_vec2gen(1d0,5d0,Q(13),Q(14),nrm)
  call d_rot2_vec2gen(-2d0,6d0,Q(15),Q(16),nrm)
  call d_rot2_vec2gen(3d0,-2d0,Q(17),Q(18),nrm)
  call d_rot2_vec2gen(9d0,1d0,Q(19),Q(20),nrm)
  call d_rot2_vec2gen(8d0,5d0,Q(21),Q(22),nrm)
  call d_rot2_vec2gen(-7d0,-6d0,Q(23),Q(24),nrm)
  call d_rot2_vec2gen(6d0,2d0,Q(25),Q(26),nrm)
  call d_rot2_vec2gen(5d0,1d0,Q(27),Q(28),nrm)
  call d_rot2_vec2gen(4d0,-5d0,Q(29),Q(30),nrm)
  call d_rot2_vec2gen(-3d0,6d0,Q(31),Q(32),nrm)
  do jj=1,N-1
     B = A
     B(jj,jj) = Q(2*jj-1)
     B(jj+1,jj) = Q(2*jj)
     B(jj,jj+1) = -Q(2*jj)
     B(jj+1,jj+1) = Q(2*jj-1)
     H = matmul(H,B)
  end do
  do jj=1,N
     H(:,jj) = H(:,jj)*D(jj)
  end do
  Ho = H
  
  call d_orthhess_factor(N,H,Q,D,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  
  ! check D
  do jj=1,N
     if (abs(abs(D(jj))-1d0)>tol) then
        call u_test_failed(__LINE__)
     end if
  end do

  H2 = A
  do jj=1,N-1
     B = A
     B(jj,jj) = Q(2*jj-1)
     B(jj+1,jj) = Q(2*jj)
     B(jj,jj+1) = -Q(2*jj)
     B(jj+1,jj+1) = Q(2*jj-1)
     H2 = matmul(H2,B)
  end do
  do jj=1,N
     H2(:,jj) = H2(:,jj)*D(jj)
  end do

  ! compare H2 with Ho
  do ii=1,N-1
     do jj=1,ii+1
        if (abs(Ho(jj,ii)-H2(jj,ii))>tol) then
           if (DEBUG) then
              print*, "H2 does not agree with H in position", jj,",", ii
           end if
           print*, "H2 does not agree with H in position", jj,",", ii
           call u_test_failed(__LINE__)
        end if
     end do
  end do

  deallocate(H,H2,Ho,A,B)  

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_d_orthhess_factor
