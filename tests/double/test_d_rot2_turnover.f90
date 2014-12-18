#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_rot2_turnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine d_rot2_turnover (turnover). The following 
! tests are run:
!
! 1) three random rotations
! 2) one almost diagonal rotation
! 3) two almost diagonal rotations
! 4) three almost diagonal rotations
! 5) one almost anti-diagonal rotation
! 6) two almost anti-diagonal rotations
! 7) three almost anti-diagonal rotations
! 8) repeating 2), 3), 5), and 6) with one/two exact (anti-)diagonal rotations
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! If DEBUG mode is activated a histogram of the accuracy of the turnover is 
! printed. The program compares the histogram to a reference histogram. If the
! histogram is equal or better than the reference, then the test is passed.
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_rot2_turnover
  
  implicit none
  
  ! parameter
  integer, parameter :: nt = 120000 ! number of testcases
  integer :: accum = 100 ! number of successive turnovers
                         ! (deterioration)
  real(8) :: tol = 200d0*epsilon(1d0) ! accuracy of turnover

  ! compute variables
  logical :: pass_all = .TRUE.
  logical :: pass_cur
  integer :: ii, n, jj, INFO
  integer, allocatable :: seed(:)
  real(8) :: Q1(2),Q2(2),Q3(2)
  real(8) :: time
  double precision :: rp,nrm,pi = 3.141592653589793239d0

  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! debug
  integer :: histo(7),histo2(7,8), histot(7,8), h2,ht

  ! fix seed
  INFO = 0
  ! get size of see        
  call random_seed(size = n)
  ! allocate memory for seed
  allocate(seed(n))  
  ! check allocation
  if (allocated(seed).EQV..FALSE.) then
    INFO = 1  
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"Array allocation failed",INFO,INFO)
    end if    
    stop
  end if   
  ! store seeds    
  seed = 0
  seed(1) = 377679
  if (n>=12) then
     seed(2) = 154653 
     seed(3) = 331669 
     seed(4) = 194341
     seed(5) = 1451740
     seed(6) = 3974222
     seed(7) = 1274552
     seed(8) = 4130946
     seed(9) = 3816048
     seed(10) = 4989015
     seed(11) = 933389
     seed(12) = 4989015
  end if
  ! set the generator
  call random_seed(put = seed) 
  ! free memory        
  deallocate(seed)
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start) 
  ! print banner
  call u_test_banner(__FILE__)

  ! set accum to a even number
  if (mod(accum,2)==1) then
     accum=accum+1;
  end if

  ! Arbitrary random
  pass_cur = .TRUE.
  histo = 0
  do ii=1,nt
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  histo2(1:7,1)=histo(1:7)
  if (DEBUG) then
     if (.NOT.pass_cur) then
        write(*,*) ""
        write(*,*) "Arbitrary random rotation"
        write(*,*) ""
        write(*,*) "<1e-17",histo(1)
        write(*,*) "<1e-16",histo(2)
        write(*,*) "<1e-15",histo(3)
        write(*,*) "<1e-14",histo(4)
        write(*,*) "<1e-13",histo(5)
        write(*,*) "<1e-12",histo(6)
        write(*,*) ">1e-12",histo(7)
        write(*,*) "FAILED."
        pass_all = .FALSE.
     end if
  else
     if (.NOT. pass_cur) then
        call u_test_failed(__LINE__)
     end if
  end if

  ! one close to diagonal
  pass_cur = .TRUE.
  histo = 0
  do ii=1,nt/3
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp)*1e-18,Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/3
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp)*1e-18,Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/3
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp)*1e-18,Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  histo2(1:7,2)=histo(1:7)
  if (DEBUG) then
     if (.NOT.pass_cur) then
        write(*,*) ""
        write(*,*) "Arbitrary random rotation"
        write(*,*) ""
        write(*,*) "<1e-17",histo(1)
        write(*,*) "<1e-16",histo(2)
        write(*,*) "<1e-15",histo(3)
        write(*,*) "<1e-14",histo(4)
        write(*,*) "<1e-13",histo(5)
        write(*,*) "<1e-12",histo(6)
        write(*,*) ">1e-12",histo(7)
        write(*,*) "FAILED."
        pass_all = .FALSE.
     end if
  else
     if (.NOT. pass_cur) then
        call u_test_failed(__LINE__)
     end if
  end if

  ! two close to diagonal
  pass_cur = .TRUE.
  histo = 0
  do ii=1,nt/3
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp)*1e-18,Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp)*1e-18,Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/3
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp)*1e-18,Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp)*1e-18,Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/3
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp)*1e-18,Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp)*1e-18,Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  histo2(1:7,3)=histo(1:7)
  if (DEBUG) then
     if (.NOT.pass_cur) then
        write(*,*) ""
        write(*,*) "Arbitrary random rotation"
        write(*,*) ""
        write(*,*) "<1e-17",histo(1)
        write(*,*) "<1e-16",histo(2)
        write(*,*) "<1e-15",histo(3)
        write(*,*) "<1e-14",histo(4)
        write(*,*) "<1e-13",histo(5)
        write(*,*) "<1e-12",histo(6)
        write(*,*) ">1e-12",histo(7)
        write(*,*) "FAILED."
        pass_all = .FALSE.
     end if
  else
     if (.NOT. pass_cur) then
        call u_test_failed(__LINE__)
     end if
  end if

  ! three close to diagonal
  pass_cur = .TRUE.
  histo = 0
  do ii=1,nt
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp)*1e-18,Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp)*1e-18,Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp)*1e-18,Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do 
  histo2(1:7,4)=histo(1:7)
  if (DEBUG) then
     if (.NOT.pass_cur) then
        write(*,*) ""
        write(*,*) "Arbitrary random rotation"
        write(*,*) ""
        write(*,*) "<1e-17",histo(1)
        write(*,*) "<1e-16",histo(2)
        write(*,*) "<1e-15",histo(3)
        write(*,*) "<1e-14",histo(4)
        write(*,*) "<1e-13",histo(5)
        write(*,*) "<1e-12",histo(6)
        write(*,*) ">1e-12",histo(7)
        write(*,*) "FAILED."
        pass_all = .FALSE.
     end if
  else
     if (.NOT. pass_cur) then
        call u_test_failed(__LINE__)
     end if
  end if

  ! one close to anti-diagonal
  pass_cur = .TRUE.
  histo = 0
  do ii=1,nt/3
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp)*1e-18,sin(rp),Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/3
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp)*1e-18,sin(rp),Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/3
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp)*1e-18,sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  histo2(1:7,5)=histo(1:7)
  if (DEBUG) then
     if (.NOT.pass_cur) then
        write(*,*) ""
        write(*,*) "Arbitrary random rotation"
        write(*,*) ""
        write(*,*) "<1e-17",histo(1)
        write(*,*) "<1e-16",histo(2)
        write(*,*) "<1e-15",histo(3)
        write(*,*) "<1e-14",histo(4)
        write(*,*) "<1e-13",histo(5)
        write(*,*) "<1e-12",histo(6)
        write(*,*) ">1e-12",histo(7)
        write(*,*) "FAILED."
        pass_all = .FALSE.
     end if
  else
     if (.NOT. pass_cur) then
        call u_test_failed(__LINE__)
     end if
  end if

  ! two close to anti-diagonal
  pass_cur = .TRUE.
  histo = 0
  do ii=1,nt/3
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp)*1e-18,sin(rp),Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp)*1e-18,sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/3
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp)*1e-18,sin(rp),Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp)*1e-18,sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/3
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp)*1e-18,sin(rp),Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp)*1e-18,sin(rp),Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  histo2(1:7,6)=histo(1:7)
  if (DEBUG) then
     if (.NOT.pass_cur) then
        write(*,*) ""
        write(*,*) "Arbitrary random rotation"
        write(*,*) ""
        write(*,*) "<1e-17",histo(1)
        write(*,*) "<1e-16",histo(2)
        write(*,*) "<1e-15",histo(3)
        write(*,*) "<1e-14",histo(4)
        write(*,*) "<1e-13",histo(5)
        write(*,*) "<1e-12",histo(6)
        write(*,*) ">1e-12",histo(7)
        write(*,*) "FAILED."
        pass_all = .FALSE.
     end if
  else
     if (.NOT. pass_cur) then
        call u_test_failed(__LINE__)
     end if
  end if

  ! three close to anti-diagonal
  pass_cur = .TRUE.
  histo = 0
  do ii=1,nt
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp)*1e-18,sin(rp),Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp)*1e-18,sin(rp),Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp)*1e-18,sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  histo2(1:7,7)=histo(1:7)
  if (DEBUG) then
     if (.NOT.pass_cur) then
        write(*,*) ""
        write(*,*) "Arbitrary random rotation"
        write(*,*) ""
        write(*,*) "<1e-17",histo(1)
        write(*,*) "<1e-16",histo(2)
        write(*,*) "<1e-15",histo(3)
        write(*,*) "<1e-14",histo(4)
        write(*,*) "<1e-13",histo(5)
        write(*,*) "<1e-12",histo(6)
        write(*,*) ">1e-12",histo(7)
        write(*,*) "FAILED."
        pass_all = .FALSE.
     end if
  else
     if (.NOT. pass_cur) then
        call u_test_failed(__LINE__)
     end if
  end if

  ! some rotations exact diagonal / exact anti-diagonal
  pass_cur = .TRUE.
  histo = 0
  do ii=1,nt/12
     call d_rot2_vec2gen(1d0,0d0,Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/12
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q1(1),Q1(2),nrm,INFO)
     call d_rot2_vec2gen(1d0,0d0,Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/12
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q2(1),Q2(2),nrm,INFO)
     call d_rot2_vec2gen(1d0,0d0,Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/12
     call d_rot2_vec2gen(1d0,0d0,Q1(1),Q1(2),nrm,INFO)
     call d_rot2_vec2gen(1d0,0d0,Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/12
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q1(1),Q1(2),nrm,INFO)
     call d_rot2_vec2gen(1d0,0d0,Q2(1),Q2(2),nrm,INFO)
     call d_rot2_vec2gen(1d0,0d0,Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/12
     call d_rot2_vec2gen(1d0,0d0,Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q2(1),Q2(2),nrm,INFO)
     call d_rot2_vec2gen(1d0,0d0,Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/12
     call d_rot2_vec2gen(0d0,1d0,Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/12
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q1(1),Q1(2),nrm,INFO)
     call d_rot2_vec2gen(0d0,1d0,Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/12
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q2(1),Q2(2),nrm,INFO)
     call d_rot2_vec2gen(0d0,1d0,Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/12
     call d_rot2_vec2gen(0d0,1d0,Q1(1),Q1(2),nrm,INFO)
     call d_rot2_vec2gen(0d0,1d0,Q2(1),Q2(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/12
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q1(1),Q1(2),nrm,INFO)
     call d_rot2_vec2gen(0d0,1d0,Q2(1),Q2(2),nrm,INFO)
     call d_rot2_vec2gen(0d0,1d0,Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  do ii=1,nt/12
     call d_rot2_vec2gen(0d0,1d0,Q1(1),Q1(2),nrm,INFO)
     call random_number(rp)
     rp = 2d0*pi*rp
     call d_rot2_vec2gen(cos(rp),sin(rp),Q2(1),Q2(2),nrm,INFO)
     call d_rot2_vec2gen(0d0,1d0,Q3(1),Q3(2),nrm,INFO)
     call tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)
  end do
  histo2(1:7,8)=histo(1:7)
  if (DEBUG) then
     if (.NOT.pass_cur) then
        write(*,*) ""
        write(*,*) "Arbitrary random rotation"
        write(*,*) ""
        write(*,*) "<1e-17",histo(1)
        write(*,*) "<1e-16",histo(2)
        write(*,*) "<1e-15",histo(3)
        write(*,*) "<1e-14",histo(4)
        write(*,*) "<1e-13",histo(5)
        write(*,*) "<1e-12",histo(6)
        write(*,*) ">1e-12",histo(7)
        write(*,*) "FAILED."
        pass_all = .FALSE.
     end if
  else
     if (.NOT. pass_cur) then
        call u_test_failed(__LINE__)
     end if
  end if

  if (DEBUG) then
     ! print histogram
     write(*,*) ""
     write(*,*) "<1e-17",histo2(1,:)
     write(*,*) "<1e-16",histo2(2,:)
     write(*,*) "<1e-15",histo2(3,:)
     write(*,*) "<1e-14",histo2(4,:)
     write(*,*) "<1e-13",histo2(5,:)
     write(*,*) "<1e-12",histo2(6,:)
     write(*,*) ">1e-12",histo2(7,:)
  end if

  ! reference histogram, turnover passes test if histogram is better than this one
  histot(1,:) = (/          15,       16389,       68649,      120000,         103,      118387,      120000,       66152/)!
  histot(2,:) = (/         135,        9388,        2274,           0,         555,           0,           0,        2222/)!
  histot(3,:) = (/       67830,       73499,       44915,           0,       81487,        1613,           0,       37740/)!
  histot(4,:) = (/       51146,       20372,        4162,           0,       36978,           0,           0,       13595/)!
  histot(5,:) = (/         874,         352,           0,           0,         877,           0,           0,         291/)!
  histot(6,:) = (/           0,           0,           0,           0,           0,           0,           0,           0/)!
  histot(7,:) = (/           0,           0,           0,           0,           0,           0,           0,           0/)!

  ! compare histogram
  do jj=1,8
     h2 = histo2(7,jj)
     ht = histot(7,jj)
     do ii=7,1,-1
        if (ht>h2)then !*1.001) then
           if (DEBUG) then
              pass_all = .FALSE.
           else
     write(*,*) "<1e-17",histo2(1,:)
     write(*,*) "<1e-16",histo2(2,:)
     write(*,*) "<1e-15",histo2(3,:)
     write(*,*) "<1e-14",histo2(4,:)
     write(*,*) "<1e-13",histo2(5,:)
     write(*,*) "<1e-12",histo2(6,:)
     write(*,*) ">1e-12",histo2(7,:)
              call u_test_failed(__LINE__)           
           end if
        end if
        if (ii>1) then
           ht = ht + histot(ii-1,jj)
           h2 = h2 + histo2(ii-1,jj)
        end if
     end do
  end do
 
  if (DEBUG) then
     if (.NOT. pass_all) then
        write(*,*) "At least one turnover test FAILED."
     end if
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_d_rot2_turnover

subroutine tt(Q1,Q2,Q3,accum,tol,histo,pass_cur)

  implicit none

  real(8), intent(inout) :: Q1(2),Q2(2),Q3(2)
  integer,intent(inout) :: histo(7)
  integer, intent(in) :: accum
  logical, intent(inout) :: pass_cur

  ! compute variables
  integer :: jj, INFO
  real(8) ::B(2),H(3,3),Hs(3,3),A1(3,3),A2(3,3),A3(3,3)
  real(8) :: Q1s(2),Q2s(2),Q3s(2)
  real(8) :: tol, nrm
  
  Q1s = Q1
  Q2s = Q2
  Q3s = Q3
  A1=0d0
  A2=0d0
  A3=0d0
  A1(1,1)=Q1(1)
  A1(2,1)=Q1(2)
  A1(1,2)=-Q1(2)
  A1(2,2)=Q1(1)
  A1(3,3)=1d0
  A2(2,2)=Q2(1)
  A2(3,2)=Q2(2)
  A2(2,3)=-Q2(2)
  A2(3,3)=Q2(1)
  A2(1,1)=1d0
  A3(1,1)=Q3(1)
  A3(2,1)=Q3(2)
  A3(1,2)=-Q3(2)
  A3(2,2)=Q3(1)
  A3(3,3)=1d0
  Hs = matmul(A1,matmul(A2,A3))

  call d_rot2_turnover(Q1,Q2,Q3,INFO)

  B = Q1
  Q1 = Q3
  Q3 = Q2
  Q2 = B
  A1=0d0
  A2=0d0
  A3=0d0
  A1(2,2)=Q1(1)
  A1(3,2)=Q1(2)
  A1(2,3)=-Q1(2)
  A1(3,3)=Q1(1)
  A1(1,1)=1d0
  A2(1,1)=Q2(1)
  A2(2,1)=Q2(2)
  A2(1,2)=-Q2(2)
  A2(2,2)=Q2(1)
  A2(3,3)=1d0
  A3(2,2)=Q3(1)
  A3(3,2)=Q3(2)
  A3(2,3)=-Q3(2)
  A3(3,3)=Q3(1)
  A3(1,1)=1d0
  H = matmul(A1,matmul(A2,A3))
  H = H-Hs
  nrm = sqrt(H(1,1)*H(1,1) + H(1,2)*H(1,2) + H(1,3)*H(1,3) +&
       &H(2,1)*H(2,1) + H(2,2)*H(2,2) + H(2,3)*H(2,3) +&
       H(3,1)*H(3,1) + H(3,2)*H(3,2) + H(3,3)*H(3,3))
  
  do jj=2,accum

     call d_rot2_turnover(Q1,Q2,Q3,INFO)
     
     B = Q1
     Q1 = Q3
     Q3 = Q2
     Q2 = B
  end do
  
  A1=0d0
  A2=0d0
  A3=0d0
  A1(1,1)=Q1(1)
  A1(2,1)=Q1(2)
  A1(1,2)=-Q1(2)
  A1(2,2)=Q1(1)
  A1(3,3)=1d0
  A2(2,2)=Q2(1)
  A2(3,2)=Q2(2)
  A2(2,3)=-Q2(2)
  A2(3,3)=Q2(1)
  A2(1,1)=1d0
  A3(1,1)=Q3(1)
  A3(2,1)=Q3(2)
  A3(1,2)=-Q3(2)
  A3(2,2)=Q3(1)
  A3(3,3)=1d0
  H = matmul(A1,matmul(A2,A3))
  H = H-Hs
  nrm = nrm + sqrt(H(1,1)*H(1,1) + H(1,2)*H(1,2) + H(1,3)*H(1,3) +&
       &H(2,1)*H(2,1) + H(2,2)*H(2,2) + H(2,3)*H(2,3) +&
       H(3,1)*H(3,1) + H(3,2)*H(3,2) + H(3,3)*H(3,3))

  if (nrm>tol) then
     pass_cur = .FALSE.
  end if
  if (nrm<1e-17) then
     histo(1) = histo(1) + 1
  else if (nrm<1e-16) then
     histo(2) = histo(2) + 1
  else if (nrm<1e-15) then
     histo(3) = histo(3) + 1
  else if (nrm<1e-14) then
     histo(4) = histo(4) + 1
  else if (nrm<1e-13) then
     histo(5) = histo(5) + 1
  else if (nrm<1e-12) then
     histo(6) = histo(6) + 1
  else
     histo(7) = histo(7) + 1
     write(*,*) Q1(2), Q1s(2), Q2(2), Q2s(2), Q3(2), Q3s(2)
     write(*,*) ((Q1(2)>0).AND.(Q2(2)>0).AND.(Q3(2)>0)), ((Q1s(2)>0).AND.(Q2s(2)>0).AND.(Q3s(2)>0))
  end if

end subroutine tt
