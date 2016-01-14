#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_unifact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_unifact_qr. The following tests are run:
!
! 1) Compute roots of unity and checks the residuals for various powers of 2
!
! 2) Compute of the 2x2 matrix with Q = (0,0,1) and D = (0,1,0,-1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_unifact_qr

  implicit none
  
  ! compute variables
  integer, parameter :: MPOW = 4
  integer, parameter :: N = 2**MPOW
  integer :: ii, INFO, jj, M
  real(8) :: Q(3*(N-1)), D(2*N)
  complex(8) :: H(N,N), Z(N,N)
  integer :: ITS(N-1)
  real(8) :: tol
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! Check 1)
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
   
    ! initialize Q
    Q = 0d0
    do ii=1,M-1
      Q(3*ii) = 1d0
    end do
 
    ! initialize D
    D = 0d0
    do ii=1,M-1
      D(2*ii-1) = 1d0
    end do
    D(2*M-1) = (-1d0)**(M-1)

    ! call dohfqr
    call z_unifact_qr(.TRUE.,.TRUE.,M,Q,D,M,Z(1:M,1:M),ITS,INFO)

    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if
    
    ! compute residual matrix
    H(1:M,1:M) = matmul(H(1:M,1:M),Z(1:M,1:M))
    do ii=1,M
      H(:,ii) = H(:,ii) - Z(:,ii)*cmplx(D(2*ii-1),D(2*ii),kind=8)
    end do
    
    ! set tolerance
    tol = max(10d0,dble(M))*EISCOR_DBL_EPS
    
    ! check maximum entry
    if (maxval(abs(H(1:M,1:M))) >= tol) then
      call u_test_failed(__LINE__)
    end if  
 
  end do
  
 
  ! Check 2)
  ! initialize H to be an upper hessenberg permutation matrix
  H = cmplx(0d0,0d0,kind=8)
  H(2,1) = cmplx(0d0,1d0,kind=8)
  H(1,2) = cmplx(0d0,1d0,kind=8)
   
  ! initialize Q
  Q = 0d0
  Q(3) = 1d0
 
  ! initialize D
  D = 0d0
  D(2) = 1d0
  D(4) = -1d0

  ! call dohfqr
  call z_unifact_qr(.TRUE.,.TRUE.,2,Q(1:3),D(1:4),2,Z(1:2,1:2),ITS,INFO)

  ! check INFO
  if (INFO.NE.0) then
    call u_test_failed(__LINE__)
  end if
    
  ! compute residual matrix
  H(1:2,1:2) = matmul(H(1:2,1:2),Z(1:2,1:2))
  do ii=1,2
    H(:,ii) = H(:,ii) - Z(:,ii)*cmplx(D(2*ii-1),D(2*ii),kind=8)
  end do
    
  ! set tolerance
  tol = 10d0*EISCOR_DBL_EPS
    
  ! check maximum entry
  if (maxval(abs(H(1:2,1:2))) >= tol) then
    call u_test_failed(__LINE__)
  end if  

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_unifact_qr
