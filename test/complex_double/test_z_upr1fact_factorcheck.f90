#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_factorcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_factorcheck. 
! The following tests are run:
!
! 1) check valid factorization
! 2) check N <= 2
! 3) check Q is NAN
! 4) check Q is INF
! 5) check Q is not orthogonal
! 6) check D is NAN
! 7) check D is INF
! 8) check D is not orthogonal
! 9) check R is NAN
! 10) check R is INF
! 11) check R is not orthogonal
! 12) check R for zero diagonal
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_factorcheck

  implicit none
  
  ! compute variables
  integer, parameter :: N = 3
  integer :: ii, INFO
  real(8) :: nan, inf
  real(8) :: Q(3*(N-1)), D(2*N), C(3*N), B(3*N)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! create inf
  inf = EISCOR_DBL_INF
  inf = 10d0*inf
  
  ! create nan
  nan = 0d0
  nan = 1d0/nan
  
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

  ! set valid C, B
  C = 0d0
  do ii=1,N
    C(3*ii) = 1d0
  end do
  B = C; 
  
  ! check 1)
    ! set INFO
    INFO = 0
  
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(N,Q,D,C,B,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if
    
  ! check 2)
    ! set INFO
    INFO = 0
    
    ! insert nan
    Q(1) = nan
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(N,Q,D,C,B,INFO)
    
    ! check INFO
    if (INFO.NE.-2) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset Q
    Q(1) = 0d0
    
  ! check 3)
    ! set INFO
    INFO = 0
    
    ! insert inf
    Q(1) = inf
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(N,Q,D,C,B,INFO)
    
    ! check INFO
    if (INFO.NE.-2) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset Q
    Q(1) = 0d0
    
  ! check 4)
    ! set INFO
    INFO = 0
    
    ! insert 1
    Q(1) = 1d0
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(N,Q,D,C,B,INFO)
    
    ! check INFO
    if (INFO.NE.-2) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset Q
    Q(1) = 0d0
    
  ! check 5)
    ! set INFO
    INFO = 0
    
    ! insert nan
    D(1) = nan
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(N,Q,D,C,B,INFO)
    
    ! check INFO
    if (INFO.NE.-3) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset D
    D(1) = 1d0
    
  ! check 6)
    ! set INFO
    INFO = 0
    
    ! insert inf
    D(1) = inf
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(N,Q,D,C,B,INFO)
    
    ! check INFO
    if (INFO.NE.-3) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset D
    D(1) = 1d0
    
  ! check 7)
    ! set INFO
    INFO = 0
    
    ! insert 1
    D(2) = 1d0
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(N,Q,D,C,B,INFO)
    
    ! check INFO
    if (INFO.NE.-3) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset D
    D(2) = 0d0
    
  ! check 8)
    ! set INFO
    INFO = 0
    
    ! insert nan
    C(1) = nan
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(N,Q,D,C,B,INFO)
    
    ! check INFO
    if (INFO.NE.-4) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset R
    C(1) = 1d0
    
  ! check 9)
    ! set INFO
    INFO = 0
    
    ! insert inf
    C(1) = inf
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(N,Q,D,C,B,INFO)
    
    ! check INFO
    if (INFO.NE.-4) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset R
    C(1) = 1d0
    
  ! check 10)
    ! set INFO
    INFO = 0
    
    ! insert 1
    C(1) = 1d0
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(N,Q,D,C,B,INFO)
    
    ! check INFO
    if (INFO.NE.-4) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset R
    C(1) = 0d0
    
  ! check 11)
    ! set INFO
    INFO = 0
    
    ! insert 1
    B(1) = 1d0
    B(3) = 0d0
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(N,Q,D,C,B,INFO)
    
    ! check INFO
    if (INFO.NE.-5) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset R
    B(1) = 0d0
    B(3) = 1d0
  
  ! check 12)
    ! set INFO
    INFO = 0
    
    ! insert 1
    C(1) = 1d0
    C(3) = 0d0
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(N,Q,D,C,B,INFO)
    
    ! check INFO
    if (INFO.NE.-4) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset R
    C(1) = 0d0
    C(3) = 1d0
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_upr1fact_factorcheck
