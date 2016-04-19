#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_qr. 
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
! 5) check roots of unity with inverse-hess QZ
!
! 6) check roots of unity with cmv QZ
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_qr

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2**2
  real(8) :: tol
  integer :: ii, INFO, ITS(N-1)
  logical :: P(N-2)
  real(8) :: Q(3*(N-1)), D(2*N), C(3*N), B(3*N)
  complex(8) :: V(N,N)
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
  tol = 1d2*dble(N)*EISCOR_DBL_EPS

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
    D = 0d0
    do ii=1,N
      D(2*ii-1) = 1d0
    end do
    D(2*N-1) = (-1d0)**(N-1)

    ! set valid C and B
    C = 0d0
    do ii=1,N
      C(3*ii) = -1d0
    end do
    B = -C
    
    ! call twisted QZ
    call z_upr1fact_qr(.TRUE.,.TRUE.,l_upr1fact_hess &
    ,N,P,Q,D,C,B,N,V,ITS,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

    do ii=1,(N)
      temp = -cmplx(D(2*ii-1),D(2*ii),kind=8)*B(3*ii)/C(3*ii)
      if (abs(temp**N-cmplx(1d0,0d0,kind=8)) >= tol) then
        call u_test_failed(__LINE__)
      end if
    end do

  ! end check 1)

!  ! check 2)
!    ! set INFO
!    INFO = 0
!    
!    ! set P
!    P = .TRUE.
!    
!    ! set valid Q
!    Q = 0d0
!    do ii=1,(N-1)
!      Q(3*ii) = 1d0
!    end do     
!  
!    ! set valid D
!    D = 0d0
!    do ii=1,N
!      D(2*ii-1) = 1d0
!    end do
!    D(2*N-1) = (-1d0)**(N-1)
!    D = 0d0
!
!    ! set valid C and B
!    C = 0d0
!    do ii=1,N
!      C(3*ii) = -1d0
!    end do
!    B = -C
!    
!    ! call twisted QZ
!    call z_upr1fact_qr(.FALSE.,l_upr1fact_inversehess &
!    ,N,P,Q,D,C,B,N,V,ITS,INFO)
!    
!    ! check INFO
!    if (INFO.NE.0) then
!      call u_test_failed(__LINE__)
!    end if
!
!    do ii=1,(N)
!      temp = -cmplx(D(2*ii-1),D(2*ii),kind=8)*B(3*ii)/C(3*ii)
!      if (abs(temp**N-cmplx(1d0,0d0,kind=8)) >= tol) then
!        call u_test_failed(__LINE__)
!      end if
!    end do
!
!  ! end check 2)
!
!  ! check 3)
!    ! set INFO
!    INFO = 0
!    
!    ! set P
!    do ii=1,(N-2)
!      if (mod(ii,2).EQ.1) then 
!        P(ii) = .FALSE.
!      else 
!        P(ii) = .TRUE.
!      end if
!    end do    
!
!    ! set valid Q
!    Q = 0d0
!    do ii=1,(N-1)
!      Q(3*ii) = 1d0
!    end do     
!  
!    ! set valid D
!    D = 0d0
!    do ii=1,N
!      D(2*ii-1) = 1d0
!    end do
!    D(2*N-1) = (-1d0)**(N-1)
!
!    ! set valid C and B
!    C = 0d0
!    do ii=1,N
!      C(3*ii) = -1d0
!    end do
!    B = -C
!    
!    ! call twisted QZ
!    call z_upr1fact_qr(.FALSE.,l_upr1fact_cmv &
!    ,N,P,Q,D,C,B,N,V,ITS,INFO)
!    
!    ! check INFO
!    if (INFO.NE.0) then
!      call u_test_failed(__LINE__)
!    end if
!
!    do ii=1,(N)
!      temp = -cmplx(D(2*ii-1),D(2*ii),kind=8)*B(3*ii)/C(3*ii)
!      if (abs(temp**N-cmplx(1d0,0d0,kind=8)) >= tol) then
!        call u_test_failed(__LINE__)
!      end if
!    end do
!
!  ! end check 3)
!
!  ! check 4)
!    ! set INFO
!    INFO = 0
!    
!    ! set P
!    P = .FALSE.
!    
!    ! set valid Q
!    Q = 0d0
!    do ii=1,(N-1)
!      Q(3*ii) = 1d0
!    end do     
!  
!    ! set valid D and D
!    D = 0d0
!    do ii=1,N
!      D(2*ii-1) = 1d0
!    end do
!
!    ! set valid C, B, C2 and B2
!    C = 0d0
!    do ii=1,N
!      C(3*ii) = -1d0
!    end do
!    B = -C
!    
!    ! call twisted QZ
!    call z_upr1fact_qr(.TRUE.,.FALSE.,.FALSE.,l_upr1fact_hess &
!    ,N,P,Q,D,C,B,N,V,ITS,INFO)
!    
!    ! check INFO
!    if (INFO.NE.0) then
!      call u_test_failed(__LINE__)
!    end if
!
!    do ii=1,(N)
!      temp = -cmplx(D(2*ii-1),D(2*ii),kind=8)*B(3*ii)/C(3*ii)
!      temp = -temp/cmplx(D(2*ii-1),D(2*ii),kind=8)
!      if (abs(temp**N-cmplx(1d0,0d0,kind=8)) >= tol) then
!        call u_test_failed(__LINE__)
!      end if
!    end do
!
!  ! end check 4)
!  
!  ! check 5)
!    ! set INFO
!    INFO = 0
!    
!    ! set P
!    P = .TRUE.
!    
!    ! set valid Q
!    Q = 0d0
!    do ii=1,(N-1)
!      Q(3*ii) = 1d0
!    end do     
!  
!    ! set valid D 
!    D = 0d0
!    do ii=1,N
!      D(2*ii-1) = 1d0
!    end do
!
!    ! set valid C, B, C2 and B2
!    C = 0d0
!    do ii=1,N
!      C(3*ii) = -1d0
!    end do
!    B = -C
!    
!    ! call twisted QZ
!    call z_upr1fact_qr(.TRUE.,.FALSE.,.FALSE.,l_upr1fact_inversehess &
!    ,N,P,Q,D,C,B,N,V,ITS,INFO)
!    
!    ! check INFO
!    if (INFO.NE.0) then
!      call u_test_failed(__LINE__)
!    end if
!
!    do ii=1,(N)
!      temp = -cmplx(D(2*ii-1),D(2*ii),kind=8)*B(3*ii)/C(3*ii)
!      temp = -temp/cmplx(D(2*ii-1),D(2*ii),kind=8)
!      if (abs(temp**N-cmplx(1d0,0d0,kind=8)) >= tol) then
!        call u_test_failed(__LINE__)
!      end if
!    end do
!
!  ! end check 5)
!
!  ! check 6)
!    ! set INFO
!    INFO = 0
!    
!    ! set P
!    do ii=1,(N-2)
!      if (mod(ii,2).EQ.1) then 
!        P(ii) = .FALSE.
!      else 
!        P(ii) = .TRUE.
!      end if
!    end do    
!
!    ! set valid Q
!    Q = 0d0
!    do ii=1,(N-1)
!      Q(3*ii) = 1d0
!    end do     
!  
!    ! set valid D
!    D = 0d0
!    do ii=1,N
!      D(2*ii-1) = 1d0
!    end do
!    D(2*N-1) = (-1d0)**(N-1)
!
!    ! set valid C and B
!    C = 0d0
!    do ii=1,N
!      C(3*ii) = -1d0
!    end do
!    B = -C
!    
!    ! call twisted QZ
!    call z_upr1fact_qr(.TRUE.,.FALSE.,.FALSE.,l_upr1fact_cmv &
!    ,N,P,Q,D,C,B,N,V,ITS,INFO)
!    
!    ! check INFO
!    if (INFO.NE.0) then
!      call u_test_failed(__LINE__)
!    end if
!
!    do ii=1,(N)
!      temp = -cmplx(D(2*ii-1),D(2*ii),kind=8)*B(3*ii)/C(3*ii)
!      if (abs(temp**N-cmplx(1d0,0d0,kind=8)) >= tol) then
!        call u_test_failed(__LINE__)
!      end if
!    end do
!
!  ! end check 6)

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_upr1fact_qr
