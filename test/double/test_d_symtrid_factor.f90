#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_symtrid_factor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the factorization of a tridiagonal 
! matrix and compares it with a precomputed solution.  
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_symtrid_factor

  implicit none
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! parameters
  integer, parameter :: N = 4
  logical, parameter :: sca = .FALSE.
  real(8), parameter :: tol = 2d0*EISCOR_DBL_EPS ! accuracy (tolerance)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute variables
  integer :: M
  integer :: ii, jj, INFO
  real(8) :: Q(3*N-3), Dq(2*N), D(N), E(N-1), scale
  complex(8) :: Z(N,N)
  logical :: flag
    
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)

  ! check zero matrix
  D(1) = 0d0
  D(2) = 0d0
  D(3) = 0d0
  D(4) = 0d0
  E(1) = 0d0
  E(2) = 0d0
  E(3) = 0d0

  call d_symtrid_factor(.TRUE.,.TRUE.,sca,N,D,E,Q,Dq,scale,N,Z,INFO)
  if (INFO.NE.-56) then
     call u_test_failed(__LINE__)
  end if
  
  ! check one tridiagonal matrix
  D(1) = 1d0
  D(2) = 2d0
  D(3) = 3d0
  D(4) = 4d0
  E(1) = 1d0
  E(2) = 2d0
  E(3) = 3d0

  call d_symtrid_factor(.TRUE.,.TRUE.,sca,N,D,E,Q,Dq,scale,N,Z,INFO)

  ! scale
  if ((abs(scale-1d0).GT.tol).OR.(scale.NE.scale)) then
     call u_test_failed(__LINE__)
  end if
  
  ! Z
  call z_2Darray_check(4,4,Z,flag)
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if
  
  if (abs(Z(1,1)-cmplx(1d0,0d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Z(1,2)-cmplx(0d0,0d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Z(1,3)-cmplx(0d0,0d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Z(1,4)-cmplx(0d0,0d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if

  if (abs(Z(2,1)-cmplx(0d0,0d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Z(2,2)-cmplx(0.15961737689352434d0,0.55866081912733567d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Z(2,3)-cmplx(0.78502228754006542d0,0.12221903877869293d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Z(2,4)-cmplx(0.17670001921612805d0,0d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if

  if (abs(Z(3,1)-cmplx(0d0,0d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Z(3,2)-cmplx(-0.63846950757409759d0,-0.15961737689352445d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Z(3,3)-cmplx(0.15982489686444459d0,-0.35725565181464070d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Z(3,4)-cmplx(0.61845006725644835d0,-0.17670001921612821d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if

  if (abs(Z(4,1)-cmplx(0d0,0d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Z(4,2)-cmplx(0.47885213068057331d0,0d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Z(4,3)-cmplx(-0.31024832920745099d0,0.34550382116284334d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Z(4,4)-cmplx(0.70680007686451241d0,-0.23560002562150431d0,kind=8)).GT.tol) then
     call u_test_failed(__LINE__)
  end if
     

  ! Q
  call d_1Darray_check(3*(N-1),Q,flag)
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(1)-0.66951402318200481d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(2)-0.25106775869325187d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(3)-0.69908222213656168d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(4)-0.44829720819750452d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(5)-0.22494818294809418d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(6)-0.86511729153374006d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(7)+0.10626285016678665d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(8)+0.99219791040296457d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Q(9)-6.5203629242714667d-2).GT.tol) then
     call u_test_failed(__LINE__)
  end if


  ! Dq
  call d_1Darray_check(2*N,Dq,flag)
  if (.NOT.flag) then
     call u_test_failed(__LINE__)
  end if

  if (abs(Dq(1)-0.39054984685616950d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Dq(2)-0.92058178187525663d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Dq(3)-0d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Dq(4)+1d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Dq(5)+1.4432899320127035d-15).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Dq(6)+1d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Dq(7)-0.39054984685616800d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if
  if (abs(Dq(8)-0.92058178187525719d0).GT.tol) then
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_d_symtrid_factor
