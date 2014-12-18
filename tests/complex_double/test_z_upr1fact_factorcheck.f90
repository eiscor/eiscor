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
! 2) check ALG is not QR or QZ
! 3) check N <= 2
! 4) check Q is NAN
! 5) check Q is INF
! 6) check Q is not orthogonal
! 7) check D is NAN
! 8) check D is INF
! 9) check D is not orthogonal
! 10) check R is NAN
! 11) check R is INF
! 12) check R is not orthogonal
! 13) check R is not upper-triangular
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_factorcheck

  implicit none
  
  ! compute variables
  integer, parameter :: N = 3
  integer :: ii, INFO
  character(2) :: ALG
  real(8) :: nan, inf
  real(8) :: Q(3*(N-1)), D(2,2*(N+1)), R(4,3*N)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! create inf
  inf = huge(1d0)
  inf = 10d0*inf
  
  ! create nan
  nan = 1d0
  nan = 1d0/(nan-1d0)
  
  ! set ALG
  ALG = 'QZ'
  
  ! set valid Q
  Q = 0d0
  do ii=1,(N-1)
    Q(3*ii) = 1d0
  end do     
  
  ! set valid D
  D = 0d0
  do ii=1,(N+1)
    D(:,2*ii-1) = 1d0
  end do

  ! set valid R
  R = 0d0
  do ii=1,N
    R(:,3*ii) = 1d0
  end do
  
  ! check 1)
    ! set INFO
    INFO = 0
  
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if
  
  ! check 2)
    ! set INFO
    INFO = 0
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck('NULL',N,Q,D,R,INFO)
    
    ! check INFO
    if (INFO.NE.-1) then
      call u_test_failed(__LINE__)
    end if

  ! check 3)
    ! set INFO
    INFO = 0
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(ALG,1,Q,D,R,INFO)
    
    ! check INFO
    if (INFO.NE.-2) then
      call u_test_failed(__LINE__)
    end if    
    
  ! check 4)
    ! set INFO
    INFO = 0
    
    ! insert nan
    Q(1) = nan
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    
    ! check INFO
    if (INFO.NE.-3) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset Q
    Q(1) = 0d0
    
  ! check 5)
    ! set INFO
    INFO = 0
    
    ! insert inf
    Q(1) = inf
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    
    ! check INFO
    if (INFO.NE.-3) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset Q
    Q(1) = 0d0
    
  ! check 6)
    ! set INFO
    INFO = 0
    
    ! insert 1
    Q(1) = 1d0
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    
    ! check INFO
    if (INFO.NE.-3) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset Q
    Q(1) = 0d0
    
  ! check 7)
    ! set INFO
    INFO = 0
    
    ! insert nan
    D(2,1) = nan
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    
    ! check INFO
    if (INFO.NE.-4) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset D
    D(2,1) = 1d0
    
  ! check 8)
    ! set INFO
    INFO = 0
    
    ! insert inf
    D(2,1) = inf
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    
    ! check INFO
    if (INFO.NE.-4) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset D
    D(2,1) = 1d0
    
  ! check 9)
    ! set INFO
    INFO = 0
    
    ! insert 1
    D(2,2) = 1d0
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    
    ! check INFO
    if (INFO.NE.-4) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset D
    D(2,2) = 0d0
    
  ! check 10)
    ! set INFO
    INFO = 0
    
    ! insert nan
    R(3,1) = nan
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    
    ! check INFO
    if (INFO.NE.-5) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset R
    R(3,1) = 1d0
    
  ! check 11)
    ! set INFO
    INFO = 0
    
    ! insert inf
    R(3,1) = inf
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    
    ! check INFO
    if (INFO.NE.-5) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset R
    R(3,1) = 1d0
    
  ! check 12)
    ! set INFO
    INFO = 0
    
    ! insert 1
    R(3,1) = 1d0
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    
    ! check INFO
    if (INFO.NE.-5) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset R
    R(3,1) = 0d0
    
  ! check 13)
    ! set INFO
    INFO = 0
    
    ! insert 1
    R(4,1) = 1d0
    R(4,3) = 0d0
    
    ! call z_upr1fact_factorcheck
    call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
    
    ! check INFO
    if (INFO.NE.-5) then
      call u_test_failed(__LINE__)
    end if 
    
    ! reset R
    R(4,1) = 0d0
    R(4,3) = 1d0
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_upr1fact_factorcheck
