#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_buildbulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_buildbulge. 
! The following tests are run:
!
! 1) 1st factor in Q is [0, -1; 1, 0], both triangular parts are I,
!    SHFT = 1 and P(1) = FALSE
!
! 1) 1st factor in Q is [0, -1; 1, 0], both triangular parts are I,
!    SHFT = 1 and P(1) = TRUE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_buildbulge

  implicit none
  
  ! compute variables
  integer, parameter :: N = 3
  real(8), parameter :: tol = 1d1*epsilon(1d0)
  integer :: ii, INFO
  logical :: P(N-2)
  real(8) :: Q(3*(N-1)), D1(2*(N+1)), D2(2*(N+1))
  real(8) :: C1(3*N), B1(3*N), C2(3*N) ,B2(3*N), G(3)
  complex(8) :: SHFT
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! check 1)
    ! set INFO
    INFO = 0
    
    ! set P
    P = .FALSE.
    
    ! set valid Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii-2) = 1d0
    end do
    Q(1) = 0d0
    Q(2) = 0d0
    Q(3) = 1d0     
  
    ! set valid D1, D2
    D1 = 0d0
    do ii=1,(N+1)
      D1(2*ii-1) = 1d0
    end do
    D2 = D1

    ! set valid C1, B1, C2, B2
    C1 = 0d0
    B1 = 0d0
    do ii=1,N
      C1(3*ii) = -1d0
      B1(3*ii) = 1d0
    end do
    C2 = C1
    B2 = B1
    
    ! set valid SHFT
    SHFT = cmplx(1d0,0d0,kind=8)
    
    ! call build bulge
    call z_upr1fact_buildbulge('QZ',N,1,P,Q,D1,C1,B1,D2,C2,B2,SHFT,G,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if
    
    ! check results    
    if (abs(G(1)+1d0/sqrt(2d0)) > tol) then
      call u_test_failed(__LINE__)
    end if
    if (abs(G(2)) > tol) then
      call u_test_failed(__LINE__)
    end if
    if (abs(G(3)-1d0/sqrt(2d0)) > tol) then
      call u_test_failed(__LINE__)
    end if
  ! end check 1)
  
  ! check 2)
    ! set INFO
    INFO = 0
    
    ! set P
    P = .TRUE.
    
    ! set valid Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii-2) = 1d0
    end do
    Q(1) = 0d0
    Q(2) = 0d0
    Q(3) = 1d0     
  
    ! set valid D1, D2
    D1 = 0d0
    do ii=1,(N+1)
      D1(2*ii-1) = 1d0
    end do
    D2 = D1

    ! set valid C1, B1, C2, B2
    C1 = 0d0
    B1 = 0d0
    do ii=1,N
      C1(3*ii) = -1d0
      B1(3*ii) = 1d0
    end do
    C2 = C1
    B2 = B1
    
    ! set valid SHFT
    SHFT = cmplx(1d0,0d0,kind=8)
    
    ! call build bulge
    call z_upr1fact_buildbulge('QZ',N,1,P,Q,D1,C1,B1,D2,C2,B2,SHFT,G,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if
    
    ! check results    
    if (abs(G(1)-1d0/sqrt(2d0)) > tol) then
      call u_test_failed(__LINE__)
    end if
    if (abs(G(2)) > tol) then
      call u_test_failed(__LINE__)
    end if
    if (abs(G(3)-1d0/sqrt(2d0)) > tol) then
      call u_test_failed(__LINE__)
    end if
  ! end check 2)
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_upr1fact_buildbulge
