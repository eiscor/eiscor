#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1_factoridpspike
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1_factoridpspike. 
! The following tests are run:
!
! 1) last column
! 2) 4th column
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1_factoridpspike

  implicit none
  
  ! compute variables
  integer, parameter :: N = 6 ! at least 6
  real(8), parameter :: tol = 3d1*EISCOR_DBL_EPS
  integer :: ii, N2,N3, INFO
  real(8) :: D(2*N+2), C(3*N), B(3*N)
  complex(8) :: H1(N,N), H2(N,N), U(N)
  

  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! check 1)
  N2 = 6
  U(1) = cmplx(1d0,1d0,kind=8)
  U(2) = cmplx(1d-1,-1d-1,kind=8)
  U(3) = cmplx(4d-1,3d0,kind=8)
  U(4) = cmplx(-1d1,-1d0,kind=8)
  U(5) = cmplx(-5d-1,4d-1,kind=8)
  U(6) = cmplx(-1d0,-1d0,kind=8)

  H1 = cmplx(0d0,0d0,kind=8)
  do ii=1,N2
     H1(ii,ii)=cmplx(1d0,0d0,kind=8)
  end do
  H1(:,N2) = U

  call z_upr1_factoridpspike(.TRUE.,N2,N2,U,D,C,B,INFO)

  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if

  call z_upr1fact_extracttri(.FALSE.,N2,D,C,B,H2)

  if (DEBUG) then
     print*, ""
     do ii=1,N2
        print*, H1(ii,:)
     end do     
     print*, ""
     do ii=1,N2
        print*, H2(ii,:)
     end do
     print*, ""
     do ii=1,N2
        print*, H1(ii,:)-H2(ii,:)
     end do
  end if

  ! check difference
  if (maxval(abs(H1-H2)) > tol) then
     print*, maxval(abs(H1-H2))
     call u_test_failed(__LINE__)
  end if
  

  ! check 2)
  N2 = 6
  N3 = 4
  U(1) = cmplx(1d0,1d0,kind=8)
  U(2) = cmplx(1d-1,-1d-1,kind=8)
  U(3) = cmplx(4d-1,3d0,kind=8)
  U(4) = cmplx(-1d1,-1d0,kind=8)
  U(5) = cmplx(0d0,0d0,kind=8)
  U(6) = cmplx(0d0,0d0,kind=8)

  H1 = cmplx(0d0,0d0,kind=8)
  do ii=1,N2
     H1(ii,ii)=cmplx(1d0,0d0,kind=8)
  end do
  H1(:,N3) = U

  call z_upr1_factoridpspike(.TRUE.,N2,N3,U,D,C,B,INFO)

  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if

  call z_upr1fact_extracttri(.FALSE.,N2,D,C,B,H2)

  if (DEBUG) then
     print*, ""
     do ii=1,N2
        print*, abs(H1(ii,:))
     end do     
     print*, ""
     do ii=1,N2
        print*, abs(H2(ii,:))
     end do
     print*, ""
     do ii=1,N2
        print*, abs(H1(ii,:)-H2(ii,:))
     end do
     print*, H1(5,5), H2(5,5)
  end if

     

  ! check difference
  if (maxval(abs(H1-H2)) > tol) then
     print*, maxval(abs(H1-H2))
     call u_test_failed(__LINE__)
  end if

  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_upr1_factoridpspike

