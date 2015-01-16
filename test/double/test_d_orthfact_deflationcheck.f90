#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_orthfact_deflationcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_d_orthfact_deflationcheck.
! The following tests are run:
!
! 1) no deflation in 
!    Q1=[0.8, 0.6], Q2=[0.5,sqrt(3)/2], Q3=[-1/sqrt(2),1/sqrt(2)],
!    D = diag([1, -1, -1, 1])
!
! 2) deflation in position 3
!    Q1=[0.8, 0.6], Q2=[0.5,sqrt(3)/2], Q3=[-1,eps/2], 
!    Q4=[-1/sqrt(2),1/sqrt(2)],
!    D = diag([1, -1, -1, 1, 1])
!
! in DEBUG mode additionally the check of the input data is checked
! A) N = 1
! B) str = 0
! C) str = 4 (N=4, stp = 3)
! D) stp = 0
! E) stp = 4
! F) itcnt = -1
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_orthfact_deflationcheck

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*epsilon(1d0) ! accuracy (tolerance)
  
  ! compute variables
  integer :: info, str, stp, zero, itcnt, its(4)
  real(8) :: Q(8) ! N = 5
  real(8) :: D(10), nul = 0d0

  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 1)
  Q(1) = 8d-1
  Q(2) = 6d-1
  Q(3) = 5d-1
  Q(4) = -sqrt(3d0)/2d0
  Q(5) = -1d0/sqrt(2d0)
  Q(6) = 1d0/sqrt(2d0)
  
  D(1) = 1d0
  D(2) = 0d0
  D(3) = -1d0
  D(4) = 0d0
  D(5) = -1d0
  D(6) = 0d0
  D(7) = 1d0
  D(8) = 0d0
  
  str = 1
  stp = 3
  zero = 4
  itcnt = 4
  its(1) = 0
  its(2) = 0
  its(3) = 0
  
  call d_orthfact_deflationcheck(4,str,stp,zero,Q,D,itcnt,its,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  
  ! check Q
  if (abs(Q(1)-8d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-6d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-5d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(4)+sqrt(3d0)/2d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(5)+1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(6)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check D
  if (abs(D(1)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(2))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(3)+1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(4))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(5)+1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(6))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(7)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(8))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check integers
  if (str.NE.1) then
     call u_test_failed(__LINE__)
  end if
  if (stp.NE.3) then
     call u_test_failed(__LINE__)
  end if
  if (zero.NE.4) then
     call u_test_failed(__LINE__)
  end if
  if (itcnt.NE.4) then
     call u_test_failed(__LINE__)
  end if
  if (its(1).NE.0) then
     call u_test_failed(__LINE__)
  end if
  if (its(2).NE.0) then
     call u_test_failed(__LINE__)
  end if
  if (its(3).NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  Q(1) = 8d-1
  Q(2) = 6d-1
  Q(3) = 5d-1
  Q(4) = -sqrt(3d0)/2d0
  Q(5) = -1d0
  Q(6) = tol/2d1
  Q(7) = -1d0/sqrt(2d0)
  Q(8) = 1d0/sqrt(2d0)
  
  D(1) = 1d0
  D(2) = 0d0
  D(3) = -1d0
  D(4) = 0d0
  D(5) = -1d0
  D(6) = 0d0
  D(7) = 1d0
  D(8) = 0d0
  D(9) = 1d0
  D(10) = 0d0
  
  str = 1
  stp = 4
  zero = 5
  itcnt = 4
  its(1) = 0
  its(2) = 0
  its(3) = 0
  its(4) = 0

  call d_orthfact_deflationcheck(5,str,stp,zero,Q,D,itcnt,its,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if
  
  ! check Q
  if (abs(Q(1)-8d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-6d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-5d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(4)+sqrt(3d0)/2d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(5)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(6))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(7)+1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(8)+1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check D
  if (abs(D(1)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(2))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(3)+1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(4))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(5)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(6))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(7)+1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(8))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(9)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(10))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check integers
  if (str.NE.4) then
     call u_test_failed(__LINE__)
  end if
  if (stp.NE.4) then
     call u_test_failed(__LINE__)
  end if
  if (zero.NE.3) then
     call u_test_failed(__LINE__)
  end if
  if (itcnt.NE.0) then
     call u_test_failed(__LINE__)
  end if
  if (its(1).NE.0) then
     call u_test_failed(__LINE__)
  end if
  if (its(2).NE.0) then
     call u_test_failed(__LINE__)
  end if
  if (its(3).NE.4) then
     call u_test_failed(__LINE__)
  end if
  if (its(4).NE.0) then
     call u_test_failed(__LINE__)
  end if
  

  if (DEBUG) then
     !!!!!!!!!!!!!!!!!!!!
     ! check A)
     ! N = 1
     call d_orthfact_deflationcheck(1,str,stp,zero,Q,D,itcnt,its,INFO)

     ! check info
     if (info.NE.-1) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check B)
     ! str = 0
     str = 0
     call d_orthfact_deflationcheck(4,str,stp,zero,Q,D,itcnt,its,INFO)
     ! check info
     if (info.NE.-2) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check C)
     ! str = 4
     str = 4
     call d_orthfact_deflationcheck(4,str,stp,zero,Q,D,itcnt,its,INFO)
     ! check info
     if (info.NE.-2) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check D)
     ! stp = 0
     str = 1
     stp = 0
     call d_orthfact_deflationcheck(4,str,stp,zero,Q,D,itcnt,its,INFO)
     ! check info
     if (info.NE.-3) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check E)
     ! stp = 4
     stp = 4
     call d_orthfact_deflationcheck(4,str,stp,zero,Q,D,itcnt,its,INFO)
     ! check info
     if (info.NE.-3) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check G)
     ! Q contains INF
     Q(1) = 2d0/nul
     str = 2
     stp = 3
     call d_orthfact_deflationcheck(4,str,stp,zero,Q,D,itcnt,its,INFO)
     ! check info
     if (info.NE.-5) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check H)
     ! D contains INF
     D(1) = 2d0/nul
     Q(1) = 8d-1
     call d_orthfact_deflationcheck(4,str,stp,zero,Q,D,itcnt,its,INFO)
     ! check info
     if (info.NE.-6) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check I)
     ! itcnt = -1
     stp = 3
     itcnt = -1
     D(1) = 1d0
     call d_orthfact_deflationcheck(4,str,stp,zero,Q,D,itcnt,its,INFO)
     ! check info
     if (info.NE.-7) then
        print*, info
        call u_test_failed(__LINE__)
     end if

  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_d_orthfact_deflationcheck
