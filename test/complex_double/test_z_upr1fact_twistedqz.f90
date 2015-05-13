#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_twistedqz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_twistedqz. 
! The following tests are run:
!
! 1) check roots of unity with upper-hess QR
!
! 2) check roots of unity with inverse-hess QR
!
! 3) check roots of unity with cmv QR
!
! 4) check roots of unity with upper-hess QZ
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_twistedqz

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2**3
  real(8) :: tol
  integer :: ii, INFO, ITS(N-1)
  logical :: P(N-2)
  real(8) :: Q(3*(N-1)), D1(2*(N+1)), C1(3*N), B1(3*N)
  real(8) :: D2(2*(N+1)), C2(3*N), B2(3*N)
  complex(8) :: temp,V(N,N), W(N,N)
  interface
    function l_upr1fact_hess(m,flags)
      logical :: l_upr1fact_hess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_hess
  end interface
  interface
    function l_upr1fact_inversehess(m,flags)
      logical :: l_upr1fact_inversehess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_inversehess
  end interface
  interface
    function l_upr1fact_cmv(m,flags)
      logical :: l_upr1fact_cmv
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_cmv
  end interface
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! set tolerance
  tol = 1d1*dble(N)*EISCOR_DBL_EPS

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
      Q(3*ii) = 1d0
    end do     
  
    ! set valid D
    D1 = 0d0
    do ii=1,(N+1)
      D1(2*ii-1) = 1d0
    end do
    D1(2*N-1) = (-1d0)**(N-1)
    D2 = 0d0

    ! set valid C1 and B1
    C1 = 0d0
    do ii=1,N
      C1(3*ii) = -1d0
    end do
    B1 = -C1
    
    ! call twisted QZ
    call z_upr1fact_twistedqz(.FALSE.,.FALSE.,.FALSE.,l_upr1fact_hess &
    ,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

    do ii=1,(N)
      temp = -cmplx(D1(2*ii-1),D1(2*ii),kind=8)*B1(3*ii)/C1(3*ii)
      if (abs(temp**N-cmplx(1d0,0d0,kind=8)) >= tol) then
        call u_test_failed(__LINE__)
      end if
    end do

  ! end check 1)

  ! check 2)
    ! set INFO
    INFO = 0
    
    ! set P
    P = .TRUE.
    
    ! set valid Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii) = 1d0
    end do     
  
    ! set valid D
    D1 = 0d0
    do ii=1,(N+1)
      D1(2*ii-1) = 1d0
    end do
    D1(2*N-1) = (-1d0)**(N-1)
    D2 = 0d0

    ! set valid C1 and B1
    C1 = 0d0
    do ii=1,N
      C1(3*ii) = -1d0
    end do
    B1 = -C1
    
    ! call twisted QZ
    call z_upr1fact_twistedqz(.FALSE.,.FALSE.,.FALSE.,l_upr1fact_inversehess &
    ,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

    do ii=1,(N)
      temp = -cmplx(D1(2*ii-1),D1(2*ii),kind=8)*B1(3*ii)/C1(3*ii)
      if (abs(temp**N-cmplx(1d0,0d0,kind=8)) >= tol) then
        call u_test_failed(__LINE__)
      end if
    end do

  ! end check 2)

  ! check 3)
    ! set INFO
    INFO = 0
    
    ! set P
    do ii=1,(N-2)
      if (mod(ii,2).EQ.1) then 
        P(ii) = .FALSE.
      else 
        P(ii) = .TRUE.
      end if
    end do    

    ! set valid Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii) = 1d0
    end do     
  
    ! set valid D
    D1 = 0d0
    do ii=1,(N+1)
      D1(2*ii-1) = 1d0
    end do
    D1(2*N-1) = (-1d0)**(N-1)
    D2 = 0d0

    ! set valid C1 and B1
    C1 = 0d0
    do ii=1,N
      C1(3*ii) = -1d0
    end do
    B1 = -C1
    
    ! call twisted QZ
    call z_upr1fact_twistedqz(.FALSE.,.FALSE.,.FALSE.,l_upr1fact_cmv &
    ,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

    do ii=1,(N)
      temp = -cmplx(D1(2*ii-1),D1(2*ii),kind=8)*B1(3*ii)/C1(3*ii)
      if (abs(temp**N-cmplx(1d0,0d0,kind=8)) >= tol) then
        call u_test_failed(__LINE__)
      end if
    end do

  ! end check 3)

  ! check 4)
    ! set INFO
    INFO = 0
    
    ! set P
    P = .FALSE.
    
    ! set valid Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii) = 1d0
    end do     
  
    ! set valid D1 and D2
    D1 = 0d0
    do ii=1,(N+1)
      D1(2*ii-1) = 1d0
    end do
    D2 = D1
    D2(2*N-1) = (-1d0)**(N-1)

    ! set valid C1, B1, C2 and B2
    C1 = 0d0
    do ii=1,N
      C1(3*ii) = -1d0
    end do
    B1 = -C1
    C2 = C1
    B2 = B1
    
    ! call twisted QZ
    call z_upr1fact_twistedqz(.TRUE.,.FALSE.,.FALSE.,l_upr1fact_hess &
    ,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

    do ii=1,(N)
      temp = -cmplx(D1(2*ii-1),D1(2*ii),kind=8)*B1(3*ii)/C1(3*ii)
      temp = -temp/(cmplx(D2(2*ii-1),D2(2*ii),kind=8)*B2(3*ii)/C2(3*ii))
      if (abs(temp**N-cmplx(1d0,0d0,kind=8)) >= tol) then
        call u_test_failed(__LINE__)
      end if
    end do

  ! end check 4)

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_upr1fact_twistedqz
