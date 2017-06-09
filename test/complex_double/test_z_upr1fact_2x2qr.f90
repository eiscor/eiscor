#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_2x2qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_qr on 2x2 matrices. 
! The following tests are run:
!
! 1) check roots of unity with upper-hess QR
!
! 2) check H = [1,-1;1,1]/sqrt(2)
!
! 3) nonunitary example
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_2x2qr

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2
  real(8) :: tol
  integer :: ii, INFO, ITS(N-1)
  logical :: P(1)
  real(8) :: Q(3), D(4), C(6), B(6)
  complex(8) :: V(2,2), T(2,2), H(2,2)
  interface
    function l_upr1fact_hess(m,flags)
      logical :: l_upr1fact_hess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_hess
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
    Q(3) = 1d0
  
    ! set valid D
    D = 0d0
    D(1) = 1d0
    D(3) = -1d0

    ! set valid C and B
    C = 0d0
    C(3) = -1d0
    C(6) = -1d0
    B = -C

    ! decompress
    call z_upr1fact_decompress(N,P,Q,D,C,B,H)
   
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
    if (maxval(dble(abs(H))) >= tol) then   
      call u_test_failed(__LINE__)
    end if

  ! end check 1)




  ! check 2)
    ! set INFO
    INFO = 0
    
    ! set P
    P = .FALSE.
    
    ! set valid Q
    Q = 0d0
    Q(1) = 1d0/sqrt(2d0)
    Q(3) = 1d0/sqrt(2d0)
  
    ! set valid D
    D = 0d0
    D(1) = 1d0
    D(3) = 1d0

    ! set valid C and B
    C = 0d0
    C(3) = -1d0
    C(6) = -1d0
    B = -C

    ! decompress
    call z_upr1fact_decompress(N,P,Q,D,C,B,H)
   
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
    Q(1) = 1d0/sqrt(2d0)
    Q(3) = 1d0/sqrt(2d0)
  
    ! set valid D
    D = 0d0
    D(1) = 1d0
    D(3) = 1d0

    ! set valid C and B
    C(1) = 1d0/sqrt(2d0)
    C(2) = 0d0
    C(3) = 1d0/sqrt(2d0)
    C(4) = 0d0
    C(5) = 1d0/sqrt(2d0)
    C(6) = 1d0/sqrt(2d0)
    B = -C

    ! decompress
    call z_upr1fact_decompress(N,P,Q,D,C,B,H)
   
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
    H = abs(H)
    if (maxval(dble(H)) >= tol) then   
      call u_test_failed(__LINE__)
    end if

  ! end check 3)



  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_upr1fact_2x2qr
