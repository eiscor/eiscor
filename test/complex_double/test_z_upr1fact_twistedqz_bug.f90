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
! 1) check roots of unity with upperhess QR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_twistedqz

  implicit none
  
  ! compute variables
  integer, parameter :: N = 10
  integer :: ii, jj, INFO, ITS(N-1), IND(N)
  logical :: P(N-2)
  real(8) :: Q(3*(N-1)), D1(2*(N+1)), D2(2*(N+1))
  real(8) :: C1(3*N), B1(3*N), C2(3*N) ,B2(3*N)
  real(8) :: error, err, arg
  complex(8) :: V(N,N), W(N,N)
  interface
    function l_upr1fact_upperhess(m,flags)
      logical :: l_upr1fact_upperhess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_upperhess
  end interface
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! random seed
  call u_randomseed_initialize(INFO)
  
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

    ! set valid C1, B1
    C1 = 0d0
    B1 = 0d0
    do ii=1,N
      C1(3*ii) = -1d0
      B1(3*ii) = 1d0
    end do
    
    ! call twisted QZ
    call z_upr1fact_twistedqz(.FALSE.,.FALSE.,.FALSE.,l_upr1fact_upperhess,N,P &
    ,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
   
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

    ! check accuracy
    error = 0d0
    do ii=1,N
       call z_scalar_argument(D1(2*ii-1),D1(2*ii),arg,INFO)
       arg = arg * N / 2d0 / 3.14159265358979323846264338327950d0
       err = abs(idnint(arg) - arg)
       IND(ii) = idnint(arg)
       print*, ii, idnint(arg)
       do jj=1,ii-1
          if (IND(ii).EQ.IND(jj)) then
             print*, "found the same root twice in ", ii, " and ", jj, "."
             print*, D1(2*ii-1), D1(2*ii)
             print*, D1(2*jj-1), D1(2*jj)             
             call u_test_failed(__LINE__)
          end if
       end do
       if (err > error) then
          error = err
       end if
    end do

  ! end check 1)
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
     
end program test_z_upr1fact_twistedqz
