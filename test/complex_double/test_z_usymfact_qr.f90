#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_usymfact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_usymfact_qr with the eigenvector interface
!
!     call z_usymfact_qr(VEC,ID,N,U,V,M,Z,ITS,INFO)
!
! The following tests are run:
!
! 1) Compute roots of unity and check forward error for powers of 2.
!    Also check the eigenvector residual
!
!        A*Z - Z*diag(U).
!
! 2) Compute the 2 x 2 example corresponding to eigenvalues i and -i.
!    In symmetric factored form we use
!
!        U(1) = 0, V(1) = 1, U(2) = 1, V(2) = 0,
!
!    which represents
!
!        [ 0  -1 ]
!        [ 1   0 ],
!
!    a real-positive-subdiagonal representative with eigenvalues i and -i.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_usymfact_qr

  implicit none

  ! compute variables
  integer, parameter :: MPOW = 10
  integer, parameter :: NMAX = 2**MPOW
  real(8), parameter :: twopi = 2d0*EISCOR_DBL_PI

  logical :: VEC, ID
  integer :: ii, INFO, jj, kk, M, idx
  integer :: ITS(NMAX-1), Aidx(NMAX)

  complex(8) :: U(NMAX), E(NMAX), swap
  real(8) :: V(NMAX)
  complex(8) :: H(NMAX,NMAX), Z(NMAX,NMAX)

  real(8) :: tol, small, err1, err2

  ! timing variables
  integer :: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)

  ! use eigenvector accumulation
  VEC = .TRUE.
  ID  = .TRUE.




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Check 1)
  !
  ! Roots of unity for powers of 2.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do kk = 1,MPOW

    ! set current degree
    M = 2**kk

    ! initialize U and V for the cyclic permutation representative
    U = cmplx(0d0,0d0,kind=8)
    V = 1d0

    U(M) = cmplx(sign(1d0,(-1d0)**(M-1)),0d0,kind=8)
    V(M) = 0d0

    ! explicitly form the input matrix before QR overwrites U
    call form_usymfact_matrix(M,U(1:M),V(1:M),H(1:M,1:M))

    ! call QR
    call z_usymfact_qr(VEC,ID,M,U(1:M),V(1:M),M,Z(1:M,1:M),ITS(1:M-1),INFO)

    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if

    ! check eigenvector residual
    !
    !     H*Z - Z*diag(U).
    !
    tol = 100d0*max(10d0,dble(M))*EISCOR_DBL_EPS
    call check_eigen_residual(M,H(1:M,1:M),Z(1:M,1:M),U(1:M),tol,__LINE__)

    ! compute argument indices for sorting
    do ii = 1,M
      Aidx(ii) = nint(dble(M)*(aimag(log(U(ii)))/twopi))
    end do

    ! sort by argument index
    do ii = 1,M
      small = 2d0*dble(M) + 1d0
      idx = ii

      do jj = ii,M
        if (dble(Aidx(jj)) < small) then
          idx = jj
          small = dble(Aidx(idx))
        end if
      end do

      Aidx(idx) = Aidx(ii)
      Aidx(ii) = int(small)

      swap = U(idx)
      U(idx) = U(ii)
      U(ii) = swap
    end do

    ! true eigenvalues
    do ii = 1,M
      small = twopi*dble(Aidx(ii))/dble(M)
      E(ii) = cmplx(cos(small),sin(small),kind=8)
    end do

    ! compute maximum forward error
    small = 0d0
    do ii = 1,M
      if (abs(U(ii)-E(ii)) > small) then
        small = abs(U(ii)-E(ii))
      end if
    end do

    ! check maximum forward error
    tol = dble(M)*EISCOR_DBL_EPS
    if (small >= tol) then
      call u_test_failed(__LINE__)
    end if

  end do




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Check 2)
  !
  ! 2 x 2 matrix with eigenvalues i and -i.
  !
  ! Symmetric factored representative:
  !
  !     U(1) = 0, V(1) = 1,
  !     U(2) = 1, V(2) = 0.
  !
  ! This represents
  !
  !     [ 0  -1 ]
  !     [ 1   0 ].
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  M = 2

  ! initialize U and V
  U = cmplx(0d0,0d0,kind=8)
  V = 0d0

  U(1) = cmplx(0d0,0d0,kind=8)
  V(1) = 1d0

  U(2) = cmplx(1d0,0d0,kind=8)
  V(2) = 0d0

  ! explicitly form the input matrix before QR overwrites U
  call form_usymfact_matrix(M,U(1:M),V(1:M),H(1:M,1:M))

  ! call QR
  call z_usymfact_qr(VEC,ID,M,U(1:M),V(1:M),M,Z(1:M,1:M),ITS(1:M-1),INFO)

  ! check INFO
  if (INFO.NE.0) then
    call u_test_failed(__LINE__)
  end if

  ! check eigenvector residual
  tol = 100d0*EISCOR_DBL_EPS
  call check_eigen_residual(M,H(1:M,1:M),Z(1:M,1:M),U(1:M),tol,__LINE__)

  ! check eigenvalues against i and -i, allowing either ordering
  err1 = abs(U(1)-cmplx(0d0, 1d0,kind=8)) &
       + abs(U(2)-cmplx(0d0,-1d0,kind=8))

  err2 = abs(U(1)-cmplx(0d0,-1d0,kind=8)) &
       + abs(U(2)-cmplx(0d0, 1d0,kind=8))

  if (min(err1,err2) >= 100d0*EISCOR_DBL_EPS) then
    call u_test_failed(__LINE__)
  end if




  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! form_usymfact_matrix
  !
  ! Forms the dense matrix represented by the symmetric unitary factorization
  !
  !     A = G1 * G2 * ... * G_{N-1} * diag(1,...,1,U(N)),
  !
  ! where
  !
  !     G_i(i:i+1,i:i+1) =
  !
  !       [ U(i)       -V(i)        ]
  !       [ V(i)   conjg(U(i))      ].
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine form_usymfact_matrix(N,U,V,A)

    implicit none

    integer, intent(in) :: N
    complex(8), intent(in) :: U(N)
    real(8), intent(in) :: V(N)
    complex(8), intent(out) :: A(N,N)

    complex(8) :: zero, one
    complex(8) :: a1, a2
    integer :: ii, kk

    zero = cmplx(0d0,0d0,kind=8)
    one  = cmplx(1d0,0d0,kind=8)

    ! A = identity
    A = zero
    do ii = 1,N
      A(ii,ii) = one
    end do

    ! multiply on the right by G_1, ..., G_{N-1}
    do kk = 1,N-1

      do ii = 1,N

        a1 = A(ii,kk)
        a2 = A(ii,kk+1)

        A(ii,kk)   = a1*U(kk) + a2*V(kk)
        A(ii,kk+1) = -a1*V(kk) + a2*conjg(U(kk))

      end do

    end do

    ! right phase diag(1,...,1,U(N))
    do ii = 1,N
      A(ii,N) = A(ii,N)*U(N)
    end do

  end subroutine form_usymfact_matrix




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! check_eigen_residual
  !
  ! Checks
  !
  !     H*Z - Z*diag(LAMBDA)
  !
  ! entrywise.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_eigen_residual(N,H,Z,LAMBDA,tol,line)

    implicit none

    integer, intent(in) :: N
    complex(8), intent(in) :: H(N,N)
    complex(8), intent(in) :: Z(N,N)
    complex(8), intent(in) :: LAMBDA(N)
    real(8), intent(in) :: tol
    integer, intent(in) :: line

    complex(8) :: R(N,N)
    real(8) :: err
    integer :: ii

    R = matmul(H,Z)

    do ii = 1,N
      R(:,ii) = R(:,ii) - Z(:,ii)*LAMBDA(ii)
    end do

    err = maxval(abs(R))

    if (err >= tol) then
      call u_test_failed(line)
    end if

  end subroutine check_eigen_residual

end program test_z_usymfact_qr
