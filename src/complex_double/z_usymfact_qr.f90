#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_usymfact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine diagonalizes a unitary upper Hessenberg matrix stored in
! symmetric factored form as a product of N Givens rotations.
!
! NORMALIZE is passed to z_rot3_symturnover through z_usymfact_singlestep:
!
!     NORMALIZE(1) = .TRUE.  : renormalize U and V
!     NORMALIZE(2) = .TRUE.  : renormalize SIGMA
!     NORMALIZE(3) = .TRUE.  : renormalize GAMMA
!
! The three flags are independent.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_usymfact_qr(NORMALIZE,VEC,ID,N,U,V,M,Z,ITS,INFO)

  implicit none

  ! input/output variables
  logical, intent(in) :: NORMALIZE(3)
  logical, intent(in) :: VEC, ID
  integer, intent(in) :: N, M
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: V(N)
  complex(8), intent(inout) :: Z(M,N)
  integer, intent(inout) :: INFO, ITS(N-1)

  ! compute variables
  logical :: flg
  integer :: ii, jj, kk
  integer :: STR, STP, ZERO, ITMAX, ITCNT
  real(8) :: xx
  complex(8) :: nu

  ! initialize info
  INFO = 0

  ! check N
  if (N < 2) then
    INFO = -1

    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N must be >= 2.",INFO,INFO)
    end if

    return
  end if

  ! check U and V
  do ii = 1,N
    if (abs(abs(U(ii))**2 + V(ii)**2 - 1d0) > 10d0*EISCOR_DBL_EPS) then
      INFO = -2

      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"U or V is invalid.",INFO,INFO)
      end if

      return
    end if
  end do

  ! check M
  if (VEC .AND. (M < 1)) then
    INFO = -3

    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"M must be at least 1.",INFO,INFO)
    end if

    return
  end if

  ! check Z
  if (VEC .AND. .NOT.ID) then
    call z_2Darray_check(M,N,Z,flg)

    if (.NOT.flg) then
      INFO = -4

      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"Z is invalid.",INFO,INFO)
      end if

      return
    end if
  end if

  ! initialize storage
  ITS = 0

  ! initialize Schur/eigenvector accumulator
  if (VEC .AND. ID) then

    Z = cmplx(0d0,0d0,kind=8)

    do ii = 1,min(M,N)
      Z(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do

  end if

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
        xx = 5d-1*(3d0-xx)
        U(N+1-ii) = U(N+1-ii)*xx

      end do

      ! renormalize columns of Z after all eigenvalues have been computed
      if (VEC) then

        do jj = 1,N

          xx = 0d0

          do ii = 1,M
            xx = xx + dble(Z(ii,jj))**2 + aimag(Z(ii,jj))**2
          end do

          ! first-order renormalization factor for 1/sqrt(xx)
          xx = 5d-1*(3d0-xx)

          do ii = 1,M
            Z(ii,jj) = Z(ii,jj)*xx
          end do

        end do

      end if

      ! exit qr algorithm
      exit

    end if

    ! check for deflation
    call z_usymfact_deflationcheck(STP-STR+1,U(STR:STP),V(STR:STP),ZERO)

    if (ZERO.GT.0) then
      ITS(STR+ZERO-1) = ITS(STR+ZERO-1) + ITCNT
      ITCNT = 0
    end if

    ! if 1x1 block, remove and check again
    if (STP == (STR+ZERO-1)) then

      STP = STP - 1
      ZERO = 0
      STR = 1

    ! if greater than 1x1, chase a bulge
    else

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
      call z_usymfact_singlestep( &
           NORMALIZE, &
           VEC, &
           STP-STR+2, &
           U(STR:STP+1), &
           V(STR:STP+1), &
           nu, &
           M, &
           Z(:,STR:STP+1))

      ! update iteration counter
      ITCNT = ITCNT + 1

    end if

    ! if ITMAX hit
    if (kk == ITMAX) then
      INFO = 1
      ITS(STR+STP-1) = ITCNT
    end if

  end do

end subroutine z_usymfact_qr
