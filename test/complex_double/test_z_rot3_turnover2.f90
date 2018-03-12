#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_rot3_turnover2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_rot3_turnover2 (turnover). The following 
! tests are run:
!
! 1) a random test
! 2) a particular test which failed earlier 
!  2A) with double input (d0)
!  2B) with sinlge input (e0)
!  2C) with the input from 2B in double precision
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rot3_turnover2
  
  implicit none
  
  ! parameter
  integer, parameter :: nt = 120000 ! number of testcases
  integer :: accum = 100 ! number of successive turnovers
                         ! (deterioration)
  real(8) :: tol

  ! compute variables
  logical :: pass_all = .TRUE.
  logical :: pass_cur
  integer :: ii, n, jj
  integer, allocatable :: seed(:)
  real(8) :: Q1(3), Q2(3), Q3(3)
  real(8) :: time
  real(8) :: rp, rm, nrm, pi = 3.141592653589793239d0

  real(8) :: B(3)
  complex(8) :: H(3,3), Hs(3,3), Hd(3,3), A1(3,3), A2(3,3), A3(3,3)
  real(8) :: Q1s(3), Q2s(3), Q3s(3)


  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! VERBOSE
  integer :: histo(7), histo2(7,8), histot(7,8), h2, ht

  ! tol depending on accum
  tol = 2d0*accum*EISCOR_DBL_EPS ! accuracy of turnover

  ! fix seed
  
  ! get size of see        
  call random_seed(size = n)
  ! allocate memory for seed
  allocate(seed(n))  
  ! check allocation
  if (allocated(seed).EQV..FALSE.) then
    call u_test_failed(__LINE__)
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

  ! fixed 1 (test case random chosen
  !print*, "===== 1  ====="
  call random_number(rp)
  rp = 2d0*pi*rp
  call random_number(rm)
  rm = 2d0*pi*rm
  call z_rot3_vec3gen(cos(rp)*cos(rm),sin(rp)*cos(rm),sin(rm),Q1(1),Q1(2),Q1(3),nrm)
  call random_number(rp)
  rp = 2d0*pi*rp
  call random_number(rm)
  rm = 2d0*pi*rm
  call z_rot3_vec3gen(cos(rp)*cos(rm),sin(rp)*cos(rm),sin(rm),Q2(1),Q2(2),Q2(3),nrm)
  call random_number(rp)
  rp = 2d0*pi*rp
  call random_number(rm)
  rm = 2d0*pi*rm
  call z_rot3_vec3gen(cos(rp)*cos(rm),sin(rp)*cos(rm),sin(rm),Q3(1),Q3(2),Q3(3),nrm)

  ! store Q1, Q2, Q3
  Q1s = Q1
  Q2s = Q2
  Q3s = Q3
  ! compute Hs
  A1=cmplx(0d0,0d0,kind=8)
  A2=cmplx(0d0,0d0,kind=8)
  A3=cmplx(0d0,0d0,kind=8)
  A1(1,1)=cmplx(Q1(1),Q1(2),kind=8)
  A1(2,1)=Q1(3)
  A1(1,2)=-Q1(3)
  A1(2,2)=cmplx(Q1(1),-Q1(2),kind=8)
  A1(3,3)=cmplx(1d0,0d0,kind=8)
  A2(2,2)=cmplx(Q2(1),Q2(2),kind=8)
  A2(3,2)=Q2(3)
  A2(2,3)=-Q2(3)
  A2(3,3)=cmplx(Q2(1),-Q2(2),kind=8)
  A2(1,1)=cmplx(1d0,0d0,kind=8)
  A3(1,1)=cmplx(Q3(1),Q3(2),kind=8)
  A3(2,1)=Q3(3)
  A3(1,2)=-Q3(3)
  A3(2,2)=cmplx(Q3(1),-Q3(2),kind=8)
  A3(3,3)=cmplx(1d0,0d0,kind=8)
  Hs = matmul(A1,matmul(A2,A3))

  ! first turnover
  call z_rot3_turnover(Q1,Q2,Q3)

  ! switch position of rotations
  B = Q1
  Q1 = Q3
  Q3 = Q2
  Q2 = B

  ! compute H
  A1=cmplx(0d0,0d0,kind=8)
  A2=cmplx(0d0,0d0,kind=8)
  A3=cmplx(0d0,0d0,kind=8)
  A1(1,1)=cmplx(1d0,0d0,kind=8)
  A1(2,2)=cmplx(Q1(1),Q1(2),kind=8)
  A1(3,2)=Q1(3)
  A1(2,3)=-Q1(3)
  A1(3,3)=cmplx(Q1(1),-Q1(2),kind=8)
  A2(3,3)=cmplx(1d0,0d0,kind=8)
  A2(1,1)=cmplx(Q2(1),Q2(2),kind=8)
  A2(2,1)=Q2(3)
  A2(1,2)=-Q2(3)
  A2(2,2)=cmplx(Q2(1),-Q2(2),kind=8)
  A3(2,2)=cmplx(Q3(1),Q3(2),kind=8)
  A3(3,2)=Q3(3)
  A3(2,3)=-Q3(3)
  A3(3,3)=cmplx(Q3(1),-Q3(2),kind=8)
  A3(1,1)=cmplx(1d0,0d0,kind=8)

  ! first part of nrm
  H = matmul(A1,matmul(A2,A3))
  Hd = H-Hs
  nrm = sqrt(Hd(1,1)*conjg(Hd(1,1)) + Hd(1,2)*conjg(Hd(1,2)) + Hd(1,3)*conjg(Hd(1,3)) +&
       &Hd(2,1)*conjg(Hd(2,1)) + Hd(2,2)*conjg(Hd(2,2)) + Hd(2,3)*conjg(Hd(2,3)) +&
       Hd(3,1)*conjg(Hd(3,1)) + Hd(3,2)*conjg(Hd(3,2)) + Hd(3,3)*conjg(Hd(3,3)))

  
  if (nrm>1e-15) then 
     print*, ""
     print*, "||Hd||_F",nrm
     print*, "Hs"
     print*, Hs(1,:)
     print*, Hs(2,:)
     print*, Hs(3,:)
     print*, "H"
     print*, H(1,:)
     print*, H(2,:)
     print*, H(3,:)
     print*, "Hd"
     print*, Hd(1,:)
     print*, Hd(2,:)
     print*, Hd(3,:)
     call u_test_failed(__LINE__)           
  end if


  ! fixed 2A (double input, the less problematic case)
  !print*, ""
  !print*, "===== 2A ====="
  Q1(1)=0.27782287709264308d0
  Q1(2)=0.92445277493108780d0    
  Q1(3)=-0.26115419944197238d0     
  Q2(1)=0.35035947196078610d0
  Q2(2)=0.70026206456609341d0       
  Q2(3)=0.62199781457573600d0     
  Q3(1)=0.86471669054433509d0       
  Q3(2)=0.13828034916775414d0       
  Q3(3)=0.48284944871884922d0     


  ! store Q1, Q2, Q3
  Q1s = Q1
  Q2s = Q2
  Q3s = Q3
  ! compute Hs
  A1=cmplx(0d0,0d0,kind=8)
  A2=cmplx(0d0,0d0,kind=8)
  A3=cmplx(0d0,0d0,kind=8)
  A1(1,1)=cmplx(Q1(1),Q1(2),kind=8)
  A1(2,1)=Q1(3)
  A1(1,2)=-Q1(3)
  A1(2,2)=cmplx(Q1(1),-Q1(2),kind=8)
  A1(3,3)=cmplx(1d0,0d0,kind=8)
  A2(2,2)=cmplx(Q2(1),Q2(2),kind=8)
  A2(3,2)=Q2(3)
  A2(2,3)=-Q2(3)
  A2(3,3)=cmplx(Q2(1),-Q2(2),kind=8)
  A2(1,1)=cmplx(1d0,0d0,kind=8)
  A3(1,1)=cmplx(Q3(1),Q3(2),kind=8)
  A3(2,1)=Q3(3)
  A3(1,2)=-Q3(3)
  A3(2,2)=cmplx(Q3(1),-Q3(2),kind=8)
  A3(3,3)=cmplx(1d0,0d0,kind=8)
  Hs = matmul(A1,matmul(A2,A3))

  ! first turnover
  call z_rot3_turnover(Q1,Q2,Q3)

  ! switch position of rotations
  B = Q1
  Q1 = Q3
  Q3 = Q2
  Q2 = B

  ! compute H
  A1=cmplx(0d0,0d0,kind=8)
  A2=cmplx(0d0,0d0,kind=8)
  A3=cmplx(0d0,0d0,kind=8)
  A1(1,1)=cmplx(1d0,0d0,kind=8)
  A1(2,2)=cmplx(Q1(1),Q1(2),kind=8)
  A1(3,2)=Q1(3)
  A1(2,3)=-Q1(3)
  A1(3,3)=cmplx(Q1(1),-Q1(2),kind=8)
  A2(3,3)=cmplx(1d0,0d0,kind=8)
  A2(1,1)=cmplx(Q2(1),Q2(2),kind=8)
  A2(2,1)=Q2(3)
  A2(1,2)=-Q2(3)
  A2(2,2)=cmplx(Q2(1),-Q2(2),kind=8)
  A3(2,2)=cmplx(Q3(1),Q3(2),kind=8)
  A3(3,2)=Q3(3)
  A3(2,3)=-Q3(3)
  A3(3,3)=cmplx(Q3(1),-Q3(2),kind=8)
  A3(1,1)=cmplx(1d0,0d0,kind=8)

  ! first part of nrm
  H = matmul(A1,matmul(A2,A3))
  Hd = H-Hs
  nrm = sqrt(Hd(1,1)*conjg(Hd(1,1)) + Hd(1,2)*conjg(Hd(1,2)) + Hd(1,3)*conjg(Hd(1,3)) +&
       &Hd(2,1)*conjg(Hd(2,1)) + Hd(2,2)*conjg(Hd(2,2)) + Hd(2,3)*conjg(Hd(2,3)) +&
       Hd(3,1)*conjg(Hd(3,1)) + Hd(3,2)*conjg(Hd(3,2)) + Hd(3,3)*conjg(Hd(3,3)))

  
  if (nrm>1e-15) then
     print*, ""
     print*, "||Hd||_F",nrm
     print*, "Hs"
     print*, Hs(1,:)
     print*, Hs(2,:)
     print*, Hs(3,:)
     print*, "H"
     print*, H(1,:)
     print*, H(2,:)
     print*, H(3,:)
     print*, "Hd"
     print*, Hd(1,:)
     print*, Hd(2,:)
     print*, Hd(3,:)
     call u_test_failed(__LINE__)           
  end if


  ! fixed 2B (single input, the problematic case) 
  !print*, ""
  !print*, "===== 2B ====="
  Q1(1)=0.27782287709264308e0
  Q1(2)=0.92445277493108780e0    
  Q1(3)=-0.26115419944197238e0     
  Q2(1)=0.35035947196078610e0
  Q2(2)=0.70026206456609341e0       
  Q2(3)=0.62199781457573600e0     
  Q3(1)=0.86471669054433509e0       
  Q3(2)=0.13828034916775414e0       
  Q3(3)=0.48284944871884922e0     

  ! the input are not rotations (not normalized)
  call z_rot3_vec3gen(Q1(1),Q1(2),Q1(3),Q1(1),Q1(2),Q1(3),nrm)
  call z_rot3_vec3gen(Q2(1),Q2(2),Q2(3),Q2(1),Q2(2),Q2(3),nrm)
  call z_rot3_vec3gen(Q3(1),Q3(2),Q3(3),Q3(1),Q3(2),Q3(3),nrm)

  ! store Q1, Q2, Q3
  Q1s = Q1
  Q2s = Q2
  Q3s = Q3

  ! compute Hs
  A1=cmplx(0d0,0d0,kind=8)
  A2=cmplx(0d0,0d0,kind=8)
  A3=cmplx(0d0,0d0,kind=8)
  A1(1,1)=cmplx(Q1(1),Q1(2),kind=8)
  A1(2,1)=Q1(3)
  A1(1,2)=-Q1(3)
  A1(2,2)=cmplx(Q1(1),-Q1(2),kind=8)
  A1(3,3)=cmplx(1d0,0d0,kind=8)
  A2(2,2)=cmplx(Q2(1),Q2(2),kind=8)
  A2(3,2)=Q2(3)
  A2(2,3)=-Q2(3)
  A2(3,3)=cmplx(Q2(1),-Q2(2),kind=8)
  A2(1,1)=cmplx(1d0,0d0,kind=8)
  A3(1,1)=cmplx(Q3(1),Q3(2),kind=8)
  A3(2,1)=Q3(3)
  A3(1,2)=-Q3(3)
  A3(2,2)=cmplx(Q3(1),-Q3(2),kind=8)
  A3(3,3)=cmplx(1d0,0d0,kind=8)
  Hs = matmul(A1,matmul(A2,A3))

  ! first turnover
  call z_rot3_turnover(Q1,Q2,Q3)

  ! switch position of rotations
  B = Q1
  Q1 = Q3
  Q3 = Q2
  Q2 = B

  ! compute H
  A1=cmplx(0d0,0d0,kind=8)
  A2=cmplx(0d0,0d0,kind=8)
  A3=cmplx(0d0,0d0,kind=8)
  A1(1,1)=cmplx(1d0,0d0,kind=8)
  A1(2,2)=cmplx(Q1(1),Q1(2),kind=8)
  A1(3,2)=Q1(3)
  A1(2,3)=-Q1(3)
  A1(3,3)=cmplx(Q1(1),-Q1(2),kind=8)
  A2(3,3)=cmplx(1d0,0d0,kind=8)
  A2(1,1)=cmplx(Q2(1),Q2(2),kind=8)
  A2(2,1)=Q2(3)
  A2(1,2)=-Q2(3)
  A2(2,2)=cmplx(Q2(1),-Q2(2),kind=8)
  A3(2,2)=cmplx(Q3(1),Q3(2),kind=8)
  A3(3,2)=Q3(3)
  A3(2,3)=-Q3(3)
  A3(3,3)=cmplx(Q3(1),-Q3(2),kind=8)
  A3(1,1)=cmplx(1d0,0d0,kind=8)

  ! first part of nrm
  H = matmul(A1,matmul(A2,A3))
  Hd = H-Hs
  nrm = sqrt(Hd(1,1)*conjg(Hd(1,1)) + Hd(1,2)*conjg(Hd(1,2)) + Hd(1,3)*conjg(Hd(1,3)) +&
       &Hd(2,1)*conjg(Hd(2,1)) + Hd(2,2)*conjg(Hd(2,2)) + Hd(2,3)*conjg(Hd(2,3)) +&
       Hd(3,1)*conjg(Hd(3,1)) + Hd(3,2)*conjg(Hd(3,2)) + Hd(3,3)*conjg(Hd(3,3)))

  
  if (nrm>1e-15) then
     print*, ""
     print*, "||Hd||_F",nrm
     print*, "Hs"
     print*, Hs(1,:)
     print*, Hs(2,:)
     print*, Hs(3,:)
     print*, "H"
     print*, H(1,:)
     print*, H(2,:)
     print*, H(3,:)
     print*, "Hd"
     print*, Hd(1,:)
     print*, Hd(2,:)
     print*, Hd(3,:)
     !call u_test_failed(__LINE__)           
  end if


  ! fixed 2C (single input, the problematic case)
  !print*, ""
  !print*, "===== 2C ====="
  Q1(1) = 0.27782288193702698d0
  Q1(2) = 0.92445278167724609d0      
  Q1(3) = -0.26115420460700989d0    
  Q2(1) = 0.35035946965217590d0     
  Q2(2) = 0.70026206970214844d0     
  Q2(3) = 0.62199783325195312d0   
  Q3(1) = 0.86471670866012573d0
  Q3(2) = 0.13828034698963165d0
  Q3(3) = 0.48284944891929626d0 

  call z_rot3_vec3gen(Q1(1),Q1(2),Q1(3),Q1(1),Q1(2),Q1(3),nrm)
  call z_rot3_vec3gen(Q2(1),Q2(2),Q2(3),Q2(1),Q2(2),Q2(3),nrm)
  call z_rot3_vec3gen(Q3(1),Q3(2),Q3(3),Q3(1),Q3(2),Q3(3),nrm)

  ! store Q1, Q2, Q3
  Q1s = Q1
  Q2s = Q2
  Q3s = Q3

  ! compute Hs
  A1=cmplx(0d0,0d0,kind=8)
  A2=cmplx(0d0,0d0,kind=8)
  A3=cmplx(0d0,0d0,kind=8)
  A1(1,1)=cmplx(Q1(1),Q1(2),kind=8)
  A1(2,1)=Q1(3)
  A1(1,2)=-Q1(3)
  A1(2,2)=cmplx(Q1(1),-Q1(2),kind=8)
  A1(3,3)=cmplx(1d0,0d0,kind=8)
  A2(2,2)=cmplx(Q2(1),Q2(2),kind=8)
  A2(3,2)=Q2(3)
  A2(2,3)=-Q2(3)
  A2(3,3)=cmplx(Q2(1),-Q2(2),kind=8)
  A2(1,1)=cmplx(1d0,0d0,kind=8)
  A3(1,1)=cmplx(Q3(1),Q3(2),kind=8)
  A3(2,1)=Q3(3)
  A3(1,2)=-Q3(3)
  A3(2,2)=cmplx(Q3(1),-Q3(2),kind=8)
  A3(3,3)=cmplx(1d0,0d0,kind=8)
  Hs = matmul(A1,matmul(A2,A3))

  ! first turnover
  call z_rot3_turnover(Q1,Q2,Q3)

  ! switch position of rotations
  B = Q1
  Q1 = Q3
  Q3 = Q2
  Q2 = B

  ! compute H
  A1=cmplx(0d0,0d0,kind=8)
  A2=cmplx(0d0,0d0,kind=8)
  A3=cmplx(0d0,0d0,kind=8)
  A1(1,1)=cmplx(1d0,0d0,kind=8)
  A1(2,2)=cmplx(Q1(1),Q1(2),kind=8)
  A1(3,2)=Q1(3)
  A1(2,3)=-Q1(3)
  A1(3,3)=cmplx(Q1(1),-Q1(2),kind=8)
  A2(3,3)=cmplx(1d0,0d0,kind=8)
  A2(1,1)=cmplx(Q2(1),Q2(2),kind=8)
  A2(2,1)=Q2(3)
  A2(1,2)=-Q2(3)
  A2(2,2)=cmplx(Q2(1),-Q2(2),kind=8)
  A3(2,2)=cmplx(Q3(1),Q3(2),kind=8)
  A3(3,2)=Q3(3)
  A3(2,3)=-Q3(3)
  A3(3,3)=cmplx(Q3(1),-Q3(2),kind=8)
  A3(1,1)=cmplx(1d0,0d0,kind=8)

  ! first part of nrm
  H = matmul(A1,matmul(A2,A3))
  Hd = H-Hs
  nrm = sqrt(Hd(1,1)*conjg(Hd(1,1)) + Hd(1,2)*conjg(Hd(1,2)) + Hd(1,3)*conjg(Hd(1,3)) +&
       &Hd(2,1)*conjg(Hd(2,1)) + Hd(2,2)*conjg(Hd(2,2)) + Hd(2,3)*conjg(Hd(2,3)) +&
       Hd(3,1)*conjg(Hd(3,1)) + Hd(3,2)*conjg(Hd(3,2)) + Hd(3,3)*conjg(Hd(3,3)))

  
  if (nrm>1e-15) then
     print*, ""
     print*, "||Hd||_F",nrm
     print*, "Hs"
     print*, Hs(1,:)
     print*, Hs(2,:)
     print*, Hs(3,:)
     print*, "H"
     print*, H(1,:)
     print*, H(2,:)
     print*, H(3,:)
     print*, "Hd"
     print*, Hd(1,:)
     print*, Hd(2,:)
     print*, Hd(3,:)
     call u_test_failed(__LINE__)           
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_rot3_turnover2

