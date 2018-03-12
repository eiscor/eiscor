#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_3x3qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_qr on 3x3 matrices. 
! The following tests are run:
!
! 1) check roots of unity with hess QR
!
! 2) check roots of unity with inversehess QR
!
! 3) check roots of unity with cmv QR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_3x3qr

  implicit none
  
  ! compute variables
  integer, parameter :: N = 3
  real(8) :: tol
  integer :: ii, INFO, ITS(N-1)
  logical :: P(N-2)
  real(8) :: Q(6), D(6), C(9), B(9)
  complex(8) :: V(3,3), T(3,3), H(3,3)
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
  tol = 1d2*EISCOR_DBL_EPS

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
    Q(3) = 1d0
    Q(6) = 1d0
  
    ! set valid D
    D = 0d0
    D(1) = 1d0
    D(3) = 1d0
    D(5) = 1d0

    ! set valid C and B
    C = 0d0
    C(3) = -1d0
    C(6) = -1d0
    C(9) = -1d0
    B = -C

    ! decompress
    call z_upr1utri_decompress(.FALSE.,N,D,C,B,T)
   
    ! set H
    H(1,1) = cmplx(Q(1),Q(2),kind=8)
    H(2,1) = cmplx(Q(3),0d0,kind=8)
    H(1,2) = -H(2,1)
    H(2,2) = conjg(H(1,1))
    T(2:3,:) = matmul(H(1:2,1:2),T(2:3,:))
    T(1:2,:) = matmul(H(1:2,1:2),T(1:2,:))
    H = T
    
    ! call twisted qr
    call z_upr1fact_qr(.TRUE.,.TRUE.,l_upr1fact_hess &
    ,N,P,Q,D,C,B,N,V,ITS,INFO)
    
    ! decompress
    call z_upr1utri_decompress(.FALSE.,N,D,C,B,T)

    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

    ! check residual
    H = abs(matmul(H,V)-matmul(V,T))
    if (maxval(dble(H)) >= tol) then   
      call u_test_failed(__LINE__)
    end if

    ! check orthogonality
    H = matmul(conjg(transpose(V)),V)
    H(1,1) = H(1,1) - cmplx(1d0,0d0,kind=8)
    H(2,2) = H(2,2) - cmplx(1d0,0d0,kind=8)
    H(3,3) = H(3,3) - cmplx(1d0,0d0,kind=8)
    H = abs(H)
    if (maxval(dble(H)) >= tol) then   
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
    Q(3) = 1d0
    Q(6) = 1d0
  
    ! set valid D
    D = 0d0
    D(1) = 1d0
    D(3) = 1d0
    D(5) = 1d0

    ! set valid C and B
    C = 0d0
    C(3) = -1d0
    C(6) = -1d0
    C(9) = -1d0
    B = -C

    ! decompress
    call z_upr1utri_decompress(.FALSE.,N,D,C,B,T)
   
    ! set H
    H(1,1) = cmplx(Q(1),Q(2),kind=8)
    H(2,1) = cmplx(Q(3),0d0,kind=8)
    H(1,2) = -H(2,1)
    H(2,2) = conjg(H(1,1))
    T(1:2,:) = matmul(H(1:2,1:2),T(1:2,:))
    T(2:3,:) = matmul(H(1:2,1:2),T(2:3,:))
    H = T
    
    ! call twisted qr
    call z_upr1fact_qr(.TRUE.,.TRUE.,l_upr1fact_inversehess &
    ,N,P,Q,D,C,B,N,V,ITS,INFO)
    
    ! decompress
    call z_upr1utri_decompress(.FALSE.,N,D,C,B,T)

    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

    ! check residual
    H = abs(matmul(H,V)-matmul(V,T))
    if (maxval(dble(H)) >= tol) then   
      call u_test_failed(__LINE__)
    end if

    ! check orthogonality
    H = matmul(conjg(transpose(V)),V)
    H(1,1) = H(1,1) - cmplx(1d0,0d0,kind=8)
    H(2,2) = H(2,2) - cmplx(1d0,0d0,kind=8)
    H(3,3) = H(3,3) - cmplx(1d0,0d0,kind=8)
    H = abs(H)
    if (maxval(dble(H)) >= tol) then   
      call u_test_failed(__LINE__)
    end if

  ! end check 2)





  ! check 3)
    ! set INFO
    INFO = 0
    
    ! set P
    P = .FALSE.
    
    ! set valid Q
    Q = 0d0
    Q(3) = 1d0
    Q(6) = 1d0
  
    ! set valid D
    D = 0d0
    D(1) = 1d0
    D(3) = 1d0
    D(5) = 1d0

    ! set valid C and B
    C = 0d0
    C(3) = -1d0
    C(6) = -1d0
    C(9) = -1d0
    B = -C

    ! decompress
    call z_upr1utri_decompress(.FALSE.,N,D,C,B,T)
   
    ! set H
    H(1,1) = cmplx(Q(1),Q(2),kind=8)
    H(2,1) = cmplx(Q(3),0d0,kind=8)
    H(1,2) = -H(2,1)
    H(2,2) = conjg(H(1,1))
    T(2:3,:) = matmul(H(1:2,1:2),T(2:3,:))
    T(1:2,:) = matmul(H(1:2,1:2),T(1:2,:))
    H = T
    
    ! call twisted qr
    call z_upr1fact_qr(.TRUE.,.TRUE.,l_upr1fact_cmv &
    ,N,P,Q,D,C,B,N,V,ITS,INFO)
    
    ! decompress
    call z_upr1utri_decompress(.FALSE.,N,D,C,B,T)

    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

    ! check residual
    H = abs(matmul(H,V)-matmul(V,T))
    if (maxval(dble(H)) >= tol) then   
      call u_test_failed(__LINE__)
    end if

    ! check orthogonality
    H = matmul(conjg(transpose(V)),V)
    H(1,1) = H(1,1) - cmplx(1d0,0d0,kind=8)
    H(2,2) = H(2,2) - cmplx(1d0,0d0,kind=8)
    H(3,3) = H(3,3) - cmplx(1d0,0d0,kind=8)
    H = abs(H)
    if (maxval(dble(H)) >= tol) then   
      call u_test_failed(__LINE__)
    end if

  ! end check 3)





  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_upr1fact_3x3qr
