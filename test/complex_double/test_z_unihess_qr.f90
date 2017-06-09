#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_unihess_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_unihess_qr. The following tests are run:
!
! 1) Compute roots of unity and checks the residuals for various powers of 2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_unihess_qr

  implicit none
  
  ! compute variables
  integer, parameter :: MPOW = 4
  integer, parameter :: N = 2**MPOW
  integer :: ii, INFO, jj, M
  real(8) :: WORK(5*N)
  complex(8) :: Hold(N,N), H(N,N), Z(N,N)
  integer :: ITS(N-1)
  real(8) :: tol
  
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
    H = cmplx(0d0,0d0,kind=8)
    do ii=1,M-1
      H(ii+1,ii) = cmplx(1d0,0d0,kind=8)
    end do
    H(1,M) = cmplx(1d0,0d0,kind=8)
    Hold = H
    
    ! call dohfqr
    call z_unihess_qr(.TRUE.,.TRUE.,M,H(1:M,1:M),WORK,M,Z(1:M,1:M),ITS,INFO)
    
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
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_unihess_qr
