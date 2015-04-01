#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_unihess_factor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_z_unihess_factor.
! The following tests are run:
!
! 1)  Test cyclic for N = 2**i, i = 1, ..., 12
!     [ 0 0 ... 0 1]
!     [ 1 0 ... 0 0]
!     [ 0 1 ... 0 0]
!     [ ... ... ...]
!     [ 0 0 ... 1 0]
!
! 2) Q, D given => H z_unihess_factor has to 
!    reproduce Q and D, N = 17
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_unihess_factor

  implicit none
  
  ! parameter
  real(8), parameter :: tol = 1d1*EISCOR_DBL_EPS ! accuracy (tolerance)
  integer, parameter :: MMAX = 12
  integer, parameter :: NMAX = 2**12
  
  ! compute variables
  integer :: N, M, ii, jj, INFO
  real(8) :: Q(3*(NMAX-1))
  real(8) :: D(2*(NMAX)), nrm
  complex(8), allocatable :: H(:,:), A(:,:), B(:,:), H2(:,:), Ho(:,:)
  
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
     H = cmplx(0d0,0d0,kind=8)
     do jj=1,N-1
        H(jj+1,jj) = cmplx(1d0,0d0,kind=8)
     end do
     H(1,N) = cmplx(1d0,0d0,kind=8)
     
     call z_unihess_factor(N,H,Q,D,INFO)
         
     ! check info
     if (INFO.NE.0) then
        call u_test_failed(__LINE__)
     end if
  
     do jj=1,N-1
        ! check Q
        if (abs(Q(3*jj-2))>tol) then
           call u_test_failed(__LINE__)
        end if
        if (abs(Q(3*jj-1))>tol) then
           call u_test_failed(__LINE__)
        end if
        if (abs(Q(3*jj)-1d0)>tol) then
           call u_test_failed(__LINE__)
        end if
        ! check D
        if (abs(D(2*jj-1)-1d0)>tol) then
           call u_test_failed(__LINE__)
        end if
        if (abs(D(2*jj))>tol) then
           call u_test_failed(__LINE__)
        end if
     end do
     if (abs(D(2*N-1)+1d0)>tol) then
        call u_test_failed(__LINE__)
     end if
     if (abs(D(2*N))>tol) then
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
     H(jj,jj) = cmplx(1d0,0d0,kind=8)
  end do
  A = H
  
  D = 0d0
  do jj=1,8
     D(4*jj-3) =  1d0
     D(4*jj-1) = -1d0
  end do
  D(33) = 1d0
  
  call z_rot3_vec3gen(1d0,2d0,4d0,Q(1),Q(2),Q(3),nrm)
  call z_rot3_vec3gen(2d0,1d0,3d0,Q(4),Q(5),Q(6),nrm)
  call z_rot3_vec3gen(3d0,-5d0,2d0,Q(7),Q(8),Q(9),nrm)
  call z_rot3_vec3gen(-1d0,6d0,1d0,Q(10),Q(11),Q(12),nrm)
  call z_rot3_vec3gen(2d0,2d0,4d0,Q(13),Q(14),Q(15),nrm)
  call z_rot3_vec3gen(3d0,-1d0,3d0,Q(16),Q(17),Q(18),nrm)
  call z_rot3_vec3gen(1d0,5d0,2d0,Q(19),Q(20),Q(21),nrm)
  call z_rot3_vec3gen(-2d0,6d0,1d0,Q(22),Q(23),Q(24),nrm)
  call z_rot3_vec3gen(3d0,-2d0,4d0,Q(25),Q(26),Q(27),nrm)
  call z_rot3_vec3gen(9d0,1d0,3d0,Q(28),Q(29),Q(30),nrm)
  call z_rot3_vec3gen(8d0,5d0,2d0,Q(31),Q(32),Q(33),nrm)
  call z_rot3_vec3gen(-7d0,-6d0,1d0,Q(34),Q(35),Q(36),nrm)
  call z_rot3_vec3gen(6d0,2d0,4d0,Q(37),Q(38),Q(39),nrm)
  call z_rot3_vec3gen(5d0,1d0,3d0,Q(40),Q(41),Q(42),nrm)
  call z_rot3_vec3gen(4d0,-5d0,2d0,Q(43),Q(44),Q(45),nrm)
  call z_rot3_vec3gen(-3d0,6d0,1d0,Q(46),Q(47),Q(48),nrm)
  do jj=1,N-1
     B = A
     B(jj,jj) = cmplx(Q(3*jj-2),Q(3*jj-1),kind=8)
     B(jj+1,jj) = cmplx(Q(3*jj),0d0,kind=8)
     B(jj,jj+1) = -cmplx(Q(3*jj),0d0,kind=8)
     B(jj+1,jj+1) = cmplx(Q(3*jj-2),-Q(3*jj-1),kind=8)
     H = matmul(H,B)
  end do
  do jj=1,N
     H(:,jj) = H(:,jj)*D(2*jj-1)
  end do
  Ho = H
  
  call z_unihess_factor(N,H,Q,D,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  
  ! check D
  do jj=1,N
     if (abs(abs(D(2*jj-1))-1d0)>tol) then
        call u_test_failed(__LINE__)
     end if
     if (abs(D(2*jj))>tol) then
        call u_test_failed(__LINE__)
     end if
  end do

  H2 = A
  do jj=1,N-1
     B = A
     B(jj,jj) = cmplx(Q(3*jj-2),Q(3*jj-1),kind=8)
     B(jj+1,jj) = cmplx(Q(3*jj),0d0,kind=8)
     B(jj,jj+1) = -cmplx(Q(3*jj),0d0,kind=8)
     B(jj+1,jj+1) = cmplx(Q(3*jj-2),-Q(3*jj-1),kind=8)
     H2 = matmul(H2,B)
  end do
  do jj=1,N
     H2(:,jj) = H2(:,jj)*D(2*jj-1)
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
     
end program test_z_unihess_factor
