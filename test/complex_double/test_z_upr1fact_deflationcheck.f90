#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_deflationcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_deflationcheck. 
! The following tests are run:
!
! 1) deflation at top with P(1) = .FALSE.
! 1) deflation at bottom with P(N-2) = .TRUE.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_deflationcheck

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2**3
  real(8) :: tol
  integer :: ii, ZERO
  logical :: P(N-2)
  real(8) :: Q(3*(N-1)), D(2*N), C(3*N), B(3*N)
  complex(8) :: temp,V(N,N)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! set tolerance
  tol = 1d2*dble(N)*EISCOR_DBL_EPS

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! check 1)
    ! set P
    P = .FALSE.
    
    ! set valid Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii) = 1d0
    end do     
  
    ! set valid D
    D = 0d0
    do ii=1,N
      D(2*ii-1) = 1d0
    end do

    ! set valid C and B
    C = 0d0
    do ii=1,N
      C(3*ii) = -1d0
    end do
    B = -C
 
    ! set valid V
    V = cmplx(0d0,0d0,kind=8)
    do ii=1,N
      V(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do

    ! change Q(1)
    Q(1) = 0d0
    Q(2) = 1d0
    Q(3) = 0d0
    
    ! call deflation check
    call z_upr1fact_deflationcheck(.TRUE.,N,P,Q,D,C,B,N,V,ZERO)
    
    ! check ZERO
    if (ZERO.NE.1) then
      call u_test_failed(__LINE__)
    end if

    ! check Q
    if ( (Q(1).NE.1d0).OR.(Q(2).NE.0d0).OR.(Q(3).NE.0d0) ) then
      call u_test_failed(__LINE__)
    end if

    ! check D
    if ( (D(1).NE.0d0).OR.(D(2).NE.1d0).OR.(D(3).NE.0d0).OR.(D(4).NE.-1d0) ) then
      call u_test_failed(__LINE__)
    end if

    ! check C
    if ( (C(4).NE.0d0).OR.(C(5).NE.0d0).OR.(C(6).NE.-1d0) ) then
      call u_test_failed(__LINE__)
    end if

    ! check B
    if ( (B(4).NE.0d0).OR.(B(5).NE.0d0).OR.(B(6).NE.1d0) ) then
      call u_test_failed(__LINE__)
    end if

    ! check V
    if ( V(2,2).NE.cmplx(0d0,-1d0,kind=8) ) then
      call u_test_failed(__LINE__)
    end if

  ! end check 1)
  
  ! check 2)
    ! set P
    P = .TRUE.
    
    ! set valid Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii) = 1d0
    end do     
  
    ! set valid D
    D = 0d0
    do ii=1,N
      D(2*ii-1) = 1d0
    end do

    ! set valid C and B
    C = 0d0
    do ii=1,N
      C(3*ii) = -1d0
    end do
    B = -C
 
    ! set valid V
    V = cmplx(0d0,0d0,kind=8)
    do ii=1,N
      V(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do

    ! change Q(N-1)
    Q(3*N-5) = 0d0
    Q(3*N-4) = 1d0
    Q(3*N-3) = 0d0
    
    ! call deflation check
    call z_upr1fact_deflationcheck(.TRUE.,N,P,Q,D,C,B,N,V,ZERO)
    
    ! check ZERO
    if ( ZERO.NE.(N-1) ) then
      call u_test_failed(__LINE__)
    end if

    ! check Q
    ii = 3*(N-2)
    if ( (Q(ii+1).NE.1d0).OR.(Q(ii+2).NE.0d0).OR.(Q(ii+3).NE.0d0) ) then
      call u_test_failed(__LINE__)
    end if

    ! check D
    ii = 2*(N-2)
    if ( (D(ii+1).NE.0d0).OR.(D(ii+2).NE.1d0).OR. & 
         (D(ii+3).NE.0d0).OR.(D(ii+4).NE.-1d0) ) then
      call u_test_failed(__LINE__)
    end if

    ! check C
    ii = 3*(N-1)
    if ( (C(ii+1).NE.0d0).OR.(C(ii+2).NE.0d0).OR.(C(ii+3).NE.-1d0) ) then
      call u_test_failed(__LINE__)
    end if

    ! check B
    if ( (B(ii+1).NE.0d0).OR.(B(ii+2).NE.0d0).OR.(B(ii+3).NE.1d0) ) then
      call u_test_failed(__LINE__)
    end if

    ! check V
    if ( V(N-1,N-1).NE.cmplx(0d0,1d0,kind=8) ) then
      call u_test_failed(__LINE__)
    end if

  ! end check 2)
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_upr1fact_deflationcheck
