#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_2x2array_eig
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine test_z_2x2array_eig. 
! The following tests are run:
!
! 0) A = [0, 1; 1, 0], B = [1, 0; 0, 1], FLAG = .FALSE.
! 1) A = [0, 1; 1, 0], B = [1, 0; 0, 1], FLAG = .TRUE.
! 2) A = [1, 0; 0, 1], B = [0, 1; 1, 0], FLAG = .TRUE.
! 3) A = [1, 2i; 3, 4i], B = [1, 0; 0, 1], FLAG = .TRUE.
! 4) A = [1, 2i; 3, 4i], B = [1, 0; 0, 0], FLAG = .TRUE.
! 5) A = [1, 2i; 2, 4i], B = [1, 0; 0, 1], FLAG = .TRUE.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_2x2array_eig

  implicit none
  
  ! parameter
  real(8), parameter :: tol = 1d1*EISCOR_DBL_EPS 
  
  ! compute variables
  complex(8) :: A(2,2), B(2,2), Q(2,2), Z(2,2)
  complex(8) :: T(2,2), S(2,2), eye(2,2)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! set eye
  eye(1,1) = cmplx(1d0,0d0,kind=8)
  eye(2,1) = cmplx(0d0,0d0,kind=8)
  eye(1,2) = cmplx(0d0,0d0,kind=8)
  eye(2,2) = cmplx(1d0,0d0,kind=8)
  
  !!!!!!!!!!!!
  ! check 0)

  T(1,1) = cmplx(0d0,0d0,kind=8)
  T(2,1) = cmplx(1d0,0d0,kind=8)
  T(1,2) = cmplx(-1d0,0d0,kind=8)
  T(2,2) = cmplx(0d0,0d0,kind=8)
  
  A = T
  
  call z_2x2array_eig(.FALSE.,A,B,Q,Z)
  
  ! check results
  T = matmul(T,Q)
  T = matmul(conjg(transpose(Q)),T)
  
  if (maxval(abs(A-T)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(matmul(conjg(transpose(Q)),Q)-eye)) > tol) then
    call u_test_failed(__LINE__)
  end if
  
  !!!!!!!!!!!!
  ! check 1)

  T(1,1) = cmplx(0d0,0d0,kind=8)
  T(2,1) = cmplx(1d0,0d0,kind=8)
  T(1,2) = cmplx(-1d0,0d0,kind=8)
  T(2,2) = cmplx(0d0,0d0,kind=8)
  
  S(1,1) = cmplx(1d0,0d0,kind=8)
  S(2,1) = cmplx(0d0,0d0,kind=8)
  S(1,2) = cmplx(0d0,0d0,kind=8)
  S(2,2) = cmplx(1d0,0d0,kind=8)
  
  A = T
  B = S
  
  call z_2x2array_eig(.TRUE.,A,B,Q,Z)
  
  ! check results
  T = matmul(T,Z)
  T = matmul(conjg(transpose(Q)),T)
  S = matmul(S,Z)
  S = matmul(conjg(transpose(Q)),S)
  
  if (maxval(abs(A-T)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(B-S)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(matmul(conjg(transpose(Z)),Z)-eye)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(matmul(conjg(transpose(Q)),Q)-eye)) > tol) then
    call u_test_failed(__LINE__)
  end if
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 2)

  S(1,1) = cmplx(0d0,0d0,kind=8)
  S(2,1) = cmplx(1d0,0d0,kind=8)
  S(1,2) = cmplx(-1d0,0d0,kind=8)
  S(2,2) = cmplx(0d0,0d0,kind=8)
  
  T(1,1) = cmplx(1d0,0d0,kind=8)
  T(2,1) = cmplx(0d0,0d0,kind=8)
  T(1,2) = cmplx(0d0,0d0,kind=8)
  T(2,2) = cmplx(1d0,0d0,kind=8)
  
  A = T
  B = S
  
  call z_2x2array_eig(.TRUE.,A,B,Q,Z)
  
  ! check results
  T = matmul(T,Z)
  T = matmul(conjg(transpose(Q)),T)
  S = matmul(S,Z)
  S = matmul(conjg(transpose(Q)),S)
  
  if (maxval(abs(A-T)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(B-S)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(matmul(conjg(transpose(Z)),Z)-eye)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(matmul(conjg(transpose(Q)),Q)-eye)) > tol) then
    call u_test_failed(__LINE__)
  end if
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 3)

  T(1,1) = cmplx(1d0,0d0,kind=8)
  T(1,2) = cmplx(0d0,2d0,kind=8)
  T(2,1) = cmplx(3d0,0d0,kind=8)
  T(2,2) = cmplx(0d0,4d0,kind=8)
  
  S(1,1) = cmplx(1d0,0d0,kind=8)
  S(2,1) = cmplx(0d0,0d0,kind=8)
  S(1,2) = cmplx(0d0,0d0,kind=8)
  S(2,2) = cmplx(1d0,0d0,kind=8)
  
  A = T
  B = S
  
  call z_2x2array_eig(.TRUE.,A,B,Q,Z)
  
  ! check results
  T = matmul(T,Z)
  T = matmul(conjg(transpose(Q)),T)
  S = matmul(S,Z)
  S = matmul(conjg(transpose(Q)),S)
  
  if (maxval(abs(A-T)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(B-S)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(matmul(conjg(transpose(Z)),Z)-eye)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(matmul(conjg(transpose(Q)),Q)-eye)) > tol) then
    call u_test_failed(__LINE__)
  end if
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 4)

  T(1,1) = cmplx(1d0,0d0,kind=8)
  T(1,2) = cmplx(0d0,2d0,kind=8)
  T(2,1) = cmplx(3d0,0d0,kind=8)
  T(2,2) = cmplx(0d0,4d0,kind=8)
  
  S(1,1) = cmplx(1d0,0d0,kind=8)
  S(2,1) = cmplx(0d0,0d0,kind=8)
  S(1,2) = cmplx(0d0,0d0,kind=8)
  S(2,2) = cmplx(0d0,0d0,kind=8)
  
  A = T
  B = S
  
  call z_2x2array_eig(.TRUE.,A,B,Q,Z)
  
  ! check results
  T = matmul(T,Z)
  T = matmul(conjg(transpose(Q)),T)
  S = matmul(S,Z)
  S = matmul(conjg(transpose(Q)),S)
  
  if (maxval(abs(A-T)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(B-S)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(matmul(conjg(transpose(Z)),Z)-eye)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(matmul(conjg(transpose(Q)),Q)-eye)) > tol) then
    call u_test_failed(__LINE__)
  end if
  
  !!!!!!!!!!!!!!!!!!!!
  ! check 5)

  T(1,1) = cmplx(1d0,0d0,kind=8)
  T(1,2) = cmplx(0d0,2d0,kind=8)
  T(2,1) = cmplx(2d0,0d0,kind=8)
  T(2,2) = cmplx(0d0,4d0,kind=8)
  
  S(1,1) = cmplx(1d0,0d0,kind=8)
  S(2,1) = cmplx(0d0,0d0,kind=8)
  S(1,2) = cmplx(0d0,0d0,kind=8)
  S(2,2) = cmplx(1d0,0d0,kind=8)
  
  A = T
  B = S
  
  call z_2x2array_eig(.TRUE.,A,B,Q,Z)
  
  ! check results
  T = matmul(T,Z)
  T = matmul(conjg(transpose(Q)),T)
  S = matmul(S,Z)
  S = matmul(conjg(transpose(Q)),S)
  
  if (maxval(abs(A-T)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(B-S)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(matmul(conjg(transpose(Z)),Z)-eye)) > tol) then
    call u_test_failed(__LINE__)
  end if
  if (maxval(abs(matmul(conjg(transpose(Q)),Q)-eye)) > tol) then
    call u_test_failed(__LINE__)
  end if
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
     
end program test_z_2x2array_eig
