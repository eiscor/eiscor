#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_rfr3_symturnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks limiting and phase-sensitive cases in the root-free
! symmetric turnover.
!
! Tested subroutine interface:
!
!     call z_rfr3_symturnover(NORMALIZE,OMEGA,CC,SS,U,VV,RHO)
!
! where
!
!     OMEGA = GAMMA*SIGMA**2,
!     CC    = C**2,
!     SS    = S**2,
!     VV    = V**2.
!
! This test uses
!
!     NORMALIZE = 2,
!
! i.e., renormalize both
!
!     |U|^2 + VV = 1
!
! and
!
!     |OMEGA| = 1.
!
! The phase convention checked here is
!
!     OMEGAhat = -RHO*OMEGA*Z**2/abs(Z)**2,
!
! where
!
!     T = conjg(OMEGA)*U,
!
! and
!
!     Z = 1 + T,                                           if real(T) >= 0,
!     Z = (VV + 2*i*aimag(T))/(1 - conjg(T)),              otherwise.
!
! The expected outputs below were computed a priori in extended precision from
! the root-free formulas and then rounded to double precision constants.
!
! The tests are:
!
! 1) Exact singular case:
!      U = 1, VV = 0, OMEGA = -1.
!
! 2) Nearly singular real case:
!      U = 1-eps, VV = eps*(2-eps), OMEGA = -1.
!
! 3) Nearly singular complex case:
!      U = 1-eps, VV = eps*(2-eps), OMEGA close to -1.
!
! 4) Generic phase-sensitive case with non-real RHO and OMEGA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rfr3_symturnover

  implicit none

  ! parameters
  integer, parameter :: NORMALIZE = 2
  real(8), parameter :: eps = (EISCOR_DBL_EPS)
  real(8), parameter :: tol = 1000d0*(EISCOR_DBL_EPS)

  ! input variables
  real(8) :: vv, cc, ss
  complex(8) :: u, omega, rho

  ! expected output variables
  real(8) :: vv_out, cc_out, ss_out
  complex(8) :: u_out, omega_out

  ! auxiliary variables
  real(8) :: r, theta

  ! timing variables
  integer :: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)




  !!!!!!!!!!!!!!!!!!!!
  ! check 1)
  !
  ! Exact singular case:
  !
  !     omega + u = 0,
  !     vv = 0.
  !
  ! This forces the m = 0 branch.  The implementation chooses OMEGAhat = 1.

  u = cmplx(1d0,0d0,kind=8)
  vv = 0d0

  omega = cmplx(-1d0,0d0,kind=8)

  cc = 5d-1
  ss = 5d-1

  rho = cmplx(1d0,0d0,kind=8)

  omega_out = cmplx(1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  cc_out    =       0.000000000000000d+00
  ss_out    =       1.000000000000000d+00
  u_out     = cmplx(1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  vv_out    =       0.000000000000000d+00

  call check_case(NORMALIZE,omega,cc,ss,u,vv,rho, &
       omega_out,cc_out,ss_out,u_out,vv_out,tol,__LINE__)




  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  !
  ! Nearly singular real case:
  !
  !     u = 1 - eps,
  !     vv = 1 - |u|^2 = eps*(2-eps),
  !     omega = -1.
  !
  ! This tests the real cancellation branch.

  u = cmplx(1d0-eps,0d0,kind=8)
  vv = eps*(2d0-eps)

  omega = cmplx(-1d0,0d0,kind=8)

  cc = 5d-1
  ss = 5d-1

  rho = cmplx(1d0,0d0,kind=8)

  omega_out = cmplx(1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  cc_out    =       5.551115123125783d-17
  ss_out    =       9.999999999999999d-01
  u_out     = cmplx(9.999999999999999d-01, 0.000000000000000d+00, kind=8)
  vv_out    =       2.220446049250313d-16

  call check_case(NORMALIZE,omega,cc,ss,u,vv,rho, &
       omega_out,cc_out,ss_out,u_out,vv_out,tol,__LINE__)




  !!!!!!!!!!!!!!!!!!!!
  ! check 3)
  !
  ! Nearly singular complex case:
  !
  !     u = 1 - eps,
  !     vv = 1 - |u|^2 = eps*(2-eps),
  !     omega approximately -1 + 2*i*eps.
  !
  ! This tests the complex cancellation branch.

  u = cmplx(1d0-eps,0d0,kind=8)
  vv = eps*(2d0-eps)

  r = sqrt(1d0 + (2d0*eps)**2)
  omega = cmplx(-1d0/r,2d0*eps/r,kind=8)

  cc = 5d-1
  ss = 5d-1

  rho = cmplx(1d0,0d0,kind=8)

  omega_out = cmplx(-5.999999999999999d-01, -8.000000000000001d-01, kind=8)
  cc_out    =       2.775557561562890d-16
  ss_out    =       9.999999999999997d-01
  u_out     = cmplx(9.999999999999999d-01, -2.220446049250313d-16, kind=8)
  vv_out    =       2.220446049250313d-16

  call check_case(NORMALIZE,omega,cc,ss,u,vv,rho, &
       omega_out,cc_out,ss_out,u_out,vv_out,tol,__LINE__)




  !!!!!!!!!!!!!!!!!!!!
  ! check 4)
  !
  ! Generic phase-sensitive case.
  !
  ! This catches the convention
  !
  !     OMEGAhat = -RHO*OMEGA*Z**2/abs(Z)**2
  !
  ! rather than its conjugated-RHO alternative.

  theta = 7d-1
  omega = cmplx(cos(theta),sin(theta),kind=8)

  u = cmplx(4d-1,3d-1,kind=8)
  vv = 1d0 - (dble(u)**2 + aimag(u)**2)

  cc = 3d-1
  ss = 7d-1

  theta = -9d-1
  rho = cmplx(cos(theta),sin(theta),kind=8)

  omega_out = cmplx(-9.718911923360310d-01, 2.354304786123664d-01, kind=8)
  cc_out    =       4.735073491387741d-01
  ss_out    =       5.264926508612259d-01
  u_out     = cmplx(1.831199678440116d-02, 4.999754712008334d-02, kind=8)
  vv_out    =       9.971649160557431d-01

  call check_case(NORMALIZE,omega,cc,ss,u,vv,rho, &
       omega_out,cc_out,ss_out,u_out,vv_out,tol,__LINE__)




  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! check_case
  !
  ! Calls z_rfr3_symturnover and compares the result against hard-coded
  ! expected output constants.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_case(NORMALIZE,OMEGA0,CC0,SS0,U0,VV0,RHO, &
       OMEGAe,CCe,SSe,Ue,VVe,tol,line)

    implicit none

    ! input variables
    integer, intent(in) :: NORMALIZE
    complex(8), intent(in) :: OMEGA0, U0, RHO
    real(8), intent(in) :: CC0, SS0, VV0

    ! expected output variables
    complex(8), intent(in) :: OMEGAe, Ue
    real(8), intent(in) :: CCe, SSe, VVe

    ! tolerance and line
    real(8), intent(in) :: tol
    integer, intent(in) :: line

    ! computed output
    complex(8) :: OMEGA, U
    real(8) :: CC, SS, VV

    ! error
    real(8) :: nrm, err

    ! copy input
    OMEGA = OMEGA0
    CC = CC0
    SS = SS0
    U = U0
    VV = VV0

    ! perform turnover
    call z_rfr3_symturnover(NORMALIZE,OMEGA,CC,SS,U,VV,RHO)

    ! compare outputs
    nrm = abs(OMEGA-OMEGAe) + abs(CC-CCe) + abs(SS-SSe) &
        + abs(U-Ue) + abs(VV-VVe)

    if (nrm > tol) then
       call u_test_failed(line)
    end if

    ! check root-free normalization conditions
    err = abs(dble(OMEGA)**2 + aimag(OMEGA)**2 - 1d0)
    err = max(err,abs(CC + SS - 1d0))
    err = max(err,abs(dble(U)**2 + aimag(U)**2 + VV - 1d0))

    if (err > tol) then
       call u_test_failed(line)
    end if

  end subroutine check_case

end program test_z_rfr3_symturnover

