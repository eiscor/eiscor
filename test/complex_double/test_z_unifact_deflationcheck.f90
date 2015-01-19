#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_unifact_deflationcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_z_unifact_deflationcheck.
! The following tests are run:
!
! 1) Q1=[0.8, 0.6, 0.0], Q2=[0.3,0.4,-sqrt(3)/2], 
!    Q3=[-1/sqrt(3),1/sqrt(3), 1/sqrt(3)]
!    D = diag([1, 0.8+i0,6, 1/sqrt(2)+i/sqrt(2), 1])
!
! 1B) set STR = 2
!
! 2) set Q2 = [5/sqrt(60), -7/sqrt(60), 0], STR = 2
!
!
! in DEBUG mode additionally the check of the input data is checked
! A) N = 1
! B) STR = 0, STR = 4 
! C) STP = 0, STP = 4
! D) ZERO = -1, ZERO = 2
! E) ITCNT = -1
! F) Q invalid
! G) D invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_unifact_deflationcheck

  implicit none
  
  ! parameter
  real(8) :: tol = 1d1*epsilon(1d0) ! accuracy (tolerance)
  
  ! compute variables
  integer :: info, ITS(3), ITCNT, STR, STP, ZERO
  real(8) :: Q(9), Qs(9) ! N = 4
  real(8) :: D(8), nul = 0d0
  complex(8) :: H(2,2)
  
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
  Q(3) = 0d0
  Q(4) = 3d-1
  Q(5) = 4d-1
  Q(6) = -sqrt(3d0)/2d0
  Q(7) = -1d0/sqrt(3d0)
  Q(8) = 1d0/sqrt(3d0)
  Q(9) = 1d0/sqrt(3d0)
  
  D(1) = 1d0
  D(2) = 0d0
  D(3) = 0.8d0
  D(4) = 0.6d0
  D(5) = 1d0/sqrt(2d0)
  D(6) = 1d0/sqrt(2d0)
  D(7) = 1d0
  D(8) = 0d0
  
  ITS = 0
  STR = 1
  STP = 3
  ZERO = 0
  ITCNT = 42

  call z_unifact_deflationcheck(4,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)

  ! check info
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if

  ! check Q
  if (abs(Q(1)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(2))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(3))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(4)-48d-2)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(5)-14d-2)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(6)+sqrt(3d0)/2d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(7)+0.2d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(8)-1.4d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(9)-1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check D
  if (abs(D(1)-0.8d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(2)-0.6d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(3)-0.8d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(4)-0.6d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(5)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(6)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(7)-0.8d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(8)+0.6d0)>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check integers
  if (str.NE.2) then
     call u_test_failed(__LINE__)
  end if
  if (stp.NE.3) then
     call u_test_failed(__LINE__)
  end if
  if (zero.NE.1) then
     call u_test_failed(__LINE__)
  end if
  if (itcnt.NE.0) then
     call u_test_failed(__LINE__)
  end if
  if (its(1).NE.42) then
     call u_test_failed(__LINE__)
  end if
  if (its(2).NE.0) then
     call u_test_failed(__LINE__)
  end if
  if (its(3).NE.0) then
     call u_test_failed(__LINE__)
  end if

  !!!!!!!!!!!!!!!!!!!!
  ! check 1B)
  Q(1) = 8d-1
  Q(2) = 6d-1
  Q(3) = 0d0
  Q(4) = 3d-1
  Q(5) = 4d-1
  Q(6) = -sqrt(3d0)/2d0
  Q(7) = -1d0/sqrt(3d0)
  Q(8) = 1d0/sqrt(3d0)
  Q(9) = 1d0/sqrt(3d0)
  
  D(1) = 1d0
  D(2) = 0d0
  D(3) = 0.8d0
  D(4) = 0.6d0
  D(5) = 1d0/sqrt(2d0)
  D(6) = 1d0/sqrt(2d0)
  D(7) = 1d0
  D(8) = 0d0
  
  ITS = 0
  STR = 2
  STP = 3
  ZERO = 1
  ITCNT = 42

  call z_unifact_deflationcheck(4,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)

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
  if (abs(Q(3))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(4)-3d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(5)-4d-1)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(6)+sqrt(3d0)/2d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(7)+1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(8)-1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(9)-1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check D
  if (abs(D(1)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(2))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(3)-0.8d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(4)-0.6d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(5)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(6)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(7)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(8))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check integers
  if (str.NE.2) then
     call u_test_failed(__LINE__)
  end if
  if (stp.NE.3) then
     call u_test_failed(__LINE__)
  end if
  if (zero.NE.1) then
     call u_test_failed(__LINE__)
  end if
  if (itcnt.NE.42) then
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
  Q(3) = 0d0
  Q(4) = 5d0/sqrt(61d0)
  Q(5) = -6d0/sqrt(61d0)
  Q(6) = 0d0
  Q(7) = -1d0/sqrt(3d0)
  Q(8) = 1d0/sqrt(3d0)
  Q(9) = 1d0/sqrt(3d0)
  
  D(1) = 1d0
  D(2) = 0d0
  D(3) = 0.8d0
  D(4) = 0.6d0
  D(5) = 1d0/sqrt(2d0)
  D(6) = 1d0/sqrt(2d0)
  D(7) = 1d0
  D(8) = 0d0
  
  ITS = 0
  STR = 2
  STP = 3
  ZERO = 1
  ITCNT = 42

  call z_unifact_deflationcheck(4,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)

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
  if (abs(Q(3))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(4)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(5))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(6))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(7)+11d0/sqrt(183d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(8)+1d0/sqrt(183d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(9)-1d0/sqrt(3d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check D
  if (abs(D(1)-1d0)>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(2))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(3)-7.6d0/sqrt(61d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(4)+1.8d0/sqrt(61d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(5)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(6)-1d0/sqrt(2d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(7)-5d0/sqrt(61d0))>tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(D(8)-6d0/sqrt(61d0))>tol) then
     call u_test_failed(__LINE__)
  end if

  ! check integers
  if (str.NE.3) then
     call u_test_failed(__LINE__)
  end if
  if (stp.NE.3) then
     call u_test_failed(__LINE__)
  end if
  if (zero.NE.2) then
     call u_test_failed(__LINE__)
  end if
  if (itcnt.NE.0) then
     call u_test_failed(__LINE__)
  end if
  if (its(1).NE.0) then
     call u_test_failed(__LINE__)
  end if
  if (its(2).NE.42) then
     call u_test_failed(__LINE__)
  end if
  if (its(3).NE.0) then
     call u_test_failed(__LINE__)
  end if



  if (DEBUG) then
     !!!!!!!!!!!!!!!!!!!!
     ! check A)
     ! N = 1
     call z_unifact_deflationcheck(1,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)
     ! check info
     if (info.NE.-1) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check B)
     STR = 0
     call z_unifact_deflationcheck(4,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)
     ! check info
     if (info.NE.-2) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     STR = 4
     call z_unifact_deflationcheck(4,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)
     ! check info
     if (info.NE.-2) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check C)
     STR = 1
     STP = 0
     call z_unifact_deflationcheck(4,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)
     ! check info
     if (info.NE.-3) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     STP = 4
     call z_unifact_deflationcheck(4,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)
     ! check info
     if (info.NE.-3) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check D)
     STP = 3
     ZERO = -1
     call z_unifact_deflationcheck(4,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)
     ! check info
     if (info.NE.-4) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     ZERO = 2
     call z_unifact_deflationcheck(4,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)
     ! check info
     if (info.NE.-4) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check E)
     ZERO = 0
     ITCNT = -1
     call z_unifact_deflationcheck(4,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)
     ! check info
     if (info.NE.-7) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check F)
     ITCNT = 1
     Qs = Q
     Q(1) = 1d0/nul
     call z_unifact_deflationcheck(4,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)
     ! check info
     if (info.NE.-5) then
        print*, info
        call u_test_failed(__LINE__)
     end if

     !!!!!!!!!!!!!!!!!!!!
     ! check G)
     Q = Qs
     D(1) = 1d0/nul
     call z_unifact_deflationcheck(4,STR,STP,ZERO,Q,D,ITCNT,ITS,INFO)
     ! check info
     if (info.NE.-6) then
        print*, info
        call u_test_failed(__LINE__)
     end if


  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_unifact_deflationcheck
