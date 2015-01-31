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
!    SHFT = 1 and P = FALSE
!
! 2) 1st factor in Q is [0, -1; 1, 0], both triangular parts are I,
!    SHFT = 1 and P = TRUE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_buildbulge

  implicit none
  
  ! compute variables
  real(8), parameter :: tol = 1d1*EISCOR_DBL_EPS
  logical :: P
  real(8) :: Q(6), D1(4), C1(6), B1(6)
  real(8) :: G(3), D2(4), C2(6), B2(6)
  complex(8) :: SHFT
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
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
    Q(3) = 1d0
    Q(4) = 1d0     
  
    ! set valid D1 and D2
    D1 = 0d0
    D1(1) = 1d0
    D1(3) = 1d0
    D2 = D1

    ! set valid C1, B1, C2, B2
    C1 = 0d0
    C1(3) = 1d0; C1(6) = 1d0
    C2 = C1; B1 = -C1; B2 = B1

    ! set valid SHFT
    SHFT = cmplx(1d0,0d0,kind=8)
    
    ! call build bulge
    call z_upr1fact_buildbulge(.TRUE.,P,Q,D1,C1,B1,D2,C2,B2,SHFT,G)
    
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
    ! set P
    P = .TRUE.
    
    ! set valid Q
    Q = 0d0
    Q(3) = 1d0
    Q(4) = 1d0     
  
    ! set valid D1 and D2
    D1 = 0d0
    D1(1) = 1d0
    D1(3) = 1d0
    D2 = D1

    ! set valid C1, B1, C2, B2
    C1 = 0d0
    C1(3) = 1d0; C1(6) = 1d0
    C2 = C1; B1 = -C1; B2 = B1

    ! set valid SHFT
    SHFT = cmplx(1d0,0d0,kind=8)
    
    ! call build bulge
    call z_upr1fact_buildbulge(.TRUE.,P,Q,D1,C1,B1,D2,C2,B2,SHFT,G)
    
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
