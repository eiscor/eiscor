#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_urffact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine diagonalizes a unitary upper Hessenberg matrix stored in
! square-root-free symmetric factored form.
!
! The matrix is represented by
!
!   | nu  0 | | u1       -v1 |
!   |  0  1 | | v1  conj(u1) | | u2       -v2 |
!                              | v2  conj(u2) | | u3       -v3 | | 1   0 |
!                                               | v3  conj(u3) | | 0  u4 |
!
! but the square-root-free algorithm stores vi^2 rather than vi.
!
! The arrays U and VV contain
!
!      U(i) = ui,
!     VV(i) = vi^2.
!
! The input must satisfy
!
!     |U(i)|^2 + VV(i) = 1,   i = 1,...,N,
!                  VV(N) = 0.
!
! On successful output, U contains the eigenvalues.
!
! NORMALIZE is passed to z_rfr3_symturnover through z_urffact_singlestep:
!
!     NORMALIZE = 0: no renormalization
!     NORMALIZE = 1: renormalize U and VV only
!     NORMALIZE = 2: renormalize U, VV, and OMEGA
!
! Any other value defaults to NORMALIZE = 2 inside z_rfr3_symturnover.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  NORMALIZE       INTEGER
!                    controls turnover renormalization
!
!  N               INTEGER
!                    dimension of matrix
!
!  U               COMPLEX(8) array of dimension N
!                    array of complex generators for core transformations
!                    on output contains eigenvalues
!
!  VV              REAL(8) array of dimension N
!                    array of squared real generators vi^2
!
! OUTPUT VARIABLES:
!
!  ITS             INTEGER array of dimension N-1
!                    contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO =  1 implies no convergence
!                    INFO =  0 implies successful computation
!                    INFO = -1 implies N is invalid
!                    INFO = -2 implies U or VV is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_urffact_qr(NORMALIZE,N,U,VV,ITS,INFO)

  implicit none

  ! input/output variables
  integer, intent(in) :: NORMALIZE
  integer, intent(in) :: N
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: VV(N)
  integer, intent(inout) :: INFO, ITS(N-1)

  ! compute variables
  integer :: ii, kk
  integer :: STR, STP, ZERO, ITMAX, ITCNT
  real(8) :: xx
  complex(8) :: nu

  ! initialize info
  INFO = 0

  ! check N
  if (N < 2) then
    INFO = -1

    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N must be >= 2.",INFO,INFO)
    end if

    return
  end if

  ! check U and VV
  do ii = 1,N

    if (abs(abs(U(ii))**2 + VV(ii) - 1d0) > 10d0*EISCOR_DBL_EPS) then
      INFO = -2

      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"U or VV is invalid.",INFO,INFO)
      end if

      return
    end if

  end do

  ! initialize storage
  ITS = 0

  ! initialize indices
  STR = 1
  STP = N-1
  ZERO = 0
  ITMAX = 20*N
  ITCNT = 0

  ! iteration loop
  do kk = 1,ITMAX

    ! check for completion
    if (STP <= 0) then

      ! store eigenvalues in U
      do ii = 1,N-1

        U(N+1-ii) = conjg(U(N-ii))*U(N+1-ii)

        ! first-order renormalization
        xx = dble(U(N+1-ii))**2 + aimag(U(N+1-ii))**2
        U(N+1-ii) = 5d-1*U(N+1-ii)*(3d0-xx)

      end do

      exit

    end if

    ! check for deflation
    call z_urffact_deflationcheck(STP-STR+1,U(STR:STP),VV(STR:STP),ZERO)

    if (ZERO.GT.0) then
      ITS(STR+ZERO-1) = ITS(STR+ZERO-1) + ITCNT
      ITCNT = 0
    end if

    ! if 1x1 block, remove and check again
    if (STP == (STR+ZERO-1)) then

      ! update indices
      STP = STP - 1
      ZERO = 0
      STR = 1

    ! if greater than 1x1, chase a bulge
    else

      ! check ZERO
      if (ZERO.GT.0) then
        STR = STR + ZERO
      end if

      ! set nu for top deflations
      if (STR > 1) then
        nu = conjg(U(STR-1))
      else
        nu = cmplx(1d0,0d0,kind=8)
      end if

      ! perform single-shift iteration
      call z_urffact_singlestep( &
           NORMALIZE, &
           STP-STR+2, &
           U(STR:STP+1), &
           VV(STR:STP+1), &
           nu, &
           ITCNT)

      ! update iteration counter
      ITCNT = ITCNT + 1

    end if

    ! if ITMAX hit
    if (kk == ITMAX) then
      INFO = 1
      ITS(STR+STP-1) = ITCNT
    end if

  end do

end subroutine z_urffact_qr
