#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_orthhess_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine d_orthhess_qr. The following tests are run:
!
! 1) Compute roots of unity and checks the residuals for various powers of 2
!
! In DEBUG additional test with wrong COMPZ, N, H and Z are run.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_orthhess_qr

  implicit none
  
  ! compute variables
  integer, parameter :: MPOW = 4
  integer, parameter :: N = 2**MPOW
  integer :: ii, INFO, jj, M
  real(8) :: WORK(4*N), Hold(N,N), H(N,N), Z(N,N)
  integer :: ITS(N-1)
  real(8) :: tol, swap, nul=0d0
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! loop through powers of 2
  do jj=1,MPOW
  
    ! set current degree
    M = 2**jj
  
    ! initialize H to be an upper hessenberg permutation matrix
    H = 0d0
    do ii=1,M-1
      H(ii+1,ii) = 1d0
    end do
    H(1,M) = 1d0
    Hold = H
    
    ! call dohfqr
    print*, M
    call d_orthhess_qr('I',M,H(1:M,1:M),Z(1:M,1:M),ITS,WORK,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if
    
    ! compute residual matrix
    Hold(1:M,1:M) = matmul(Hold(1:M,1:M),Z(1:M,1:M))-matmul(Z(1:M,1:M),H(1:M,1:M))
    
    ! set tolerance
    tol = max(10d0,dble(M))*epsilon(1d0)
    
    ! check maximum entry
    if (maxval(abs(Hold(1:M,1:M))) >= tol) then
      call u_test_failed(__LINE__)
    end if  
 
  end do
  
  if (DEBUG) then

    call d_orthhess_qr('Q',M,H(1:M,1:M),Z(1:M,1:M),ITS,WORK,INFO)   
    ! check INFO
    if (INFO.NE.-1) then
      call u_test_failed(__LINE__)
    end if

    call d_orthhess_qr('I',1,H(1:M,1:M),Z(1:M,1:M),ITS,WORK,INFO)   
    ! check INFO
    if (INFO.NE.-2) then
      call u_test_failed(__LINE__)
    end if

    swap = H(1,1)
    H(1,1) = 1d0/nul
    call d_orthhess_qr('I',M,H(1:M,1:M),Z(1:M,1:M),ITS,WORK,INFO)   
    ! check INFO
    if (INFO.NE.-3) then
      call u_test_failed(__LINE__)
    end if

    H(1,1) = swap
    Z(1,1) = 1d0/nul
    call d_orthhess_qr('I',M,H(1:M,1:M),Z(1:M,1:M),ITS,WORK,INFO)   
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

    Z(1,1) = 1d0/nul
    call d_orthhess_qr('V',M,H(1:M,1:M),Z(1:M,1:M),ITS,WORK,INFO)   
    ! check INFO
    if (INFO.NE.-4) then
      call u_test_failed(__LINE__)
    end if
     
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_d_orthhess_qr
