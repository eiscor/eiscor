#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_orthfact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine d_orthfact_qr. The following tests are run:
!
! 1) Compute roots of unity and checks the residuals for various powers of 2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_orthfact_qr

  implicit none
  
  ! compute variables
  integer, parameter :: MPOW = 4
  integer, parameter :: N = 2**MPOW
  integer :: ii, INFO, jj, M
  real(8) :: Q(2*(N-1)), D(N)
  real(8) :: H(N,N), Z(N,N), T(N,N), block(2,2)
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
    H = 0d0
    do ii=1,M-1
      H(ii+1,ii) = 1d0
    end do
    H(1,M) = 1d0
   
    ! initialize Q
    Q = 0d0
    do ii=1,M-1
      Q(2*ii) = 1d0
    end do
 
    ! initialize D
    D = 1d0
    D(M) = (-1d0)**(M-1)

    ! call dohfqr
    call d_orthfact_qr(.TRUE.,.TRUE.,M,Q,D,M,Z(1:M,1:M),ITS,INFO)

    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

    ! initialize T
    T = 0d0
    do ii=1,M
      T(ii,ii) = 1d0
    end do
    do ii=1,(M-1)
      block(1,1) = Q(2*ii-1)
      block(2,1) = Q(2*ii)
      block(1,2) = -Q(2*ii)
      block(2,2) = Q(2*ii-1)
      T(1:M,ii:(ii+1)) = matmul(T(1:M,ii:(ii+1)),block)
    end do
    do ii=1,M
      T(1:M,ii) = T(1:M,ii)*D(ii)
    end do
    
    ! compute residual matrix
    H(1:M,1:M) = matmul(H(1:M,1:M),Z(1:M,1:M)) - matmul(Z(1:M,1:M),T(1:M,1:M))
    
    ! set tolerance
    tol = max(10d0,dble(M))*EISCOR_DBL_EPS
    
    ! check maximum entry
    if (maxval(abs(H(1:M,1:M))) >= tol) then
      call u_test_failed(__LINE__)
    end if  
 
  end do
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_d_orthfact_qr
