#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_urffact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_urffact_qr. The following tests are run:
!
! 1) Compute roots of unity and check forward error for various powers of 2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_urffact_qr

  implicit none
  
  ! compute variables
  integer, parameter :: MPOW = 10
  integer, parameter :: N = 2**MPOW
  real(8), parameter :: twopi = 2d0*EISCOR_DBL_PI
  integer :: ii, INFO, jj, kk, M, id
  complex(8) :: U(N), E(N), swap
  real(8) :: VV(N)
  integer :: ITS(N-1), A(N)
  real(8) :: tol, small
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! Check 1)
  ! loop through powers of 2
  do kk=1,MPOW
  
    ! set current degree
    M = 2**kk
  
    ! initialize U and VV
    U = cmplx(0d0,0d0,kind=8)
    U(M) = cmplx(sign(1d0,(-1d0)**(M-1)),0d0,kind=8)
    VV = 1d0
    VV(M) = 0d0

    ! call dohfqr
    call z_urffact_qr(M,U,VV,ITS,INFO)

    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if
    
    ! compute argument
    do ii = 1,M
      A(ii) = nint(dble(M)*(aimag(log(U(ii)))/twopi))
    end do
  
    ! sort by argument
    do ii = 1,M
      small = 2*M+1
      id = ii
      do jj = ii,M
        if ( A(jj) < small ) then
          id = jj
          small = A(id)
        end if
      end do
      A(id) = A(ii)
      A(ii) = small
      swap = U(id)
      U(id) = U(ii)
      U(ii) = swap
    end do
      
    ! set tolerance
    tol = dble(M)*EISCOR_DBL_EPS
   
    ! true eigenvalues
    do ii = 1,M
      small = twopi*dble(A(ii))/dble(M)
      E(ii) = cmplx(cos(small),sin(small),kind=8)
    end do
  
    ! compute maximum forward error
    small = 0d0
    do ii = 1,M
      if (abs(U(ii)-E(ii)) > small) then
        small = abs(U(ii)-E(ii))
      end if
    end do

    ! check maximum error
    if (small >= tol) then
      call u_test_failed(__LINE__)
    end if  
 
  end do

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_urffact_qr
