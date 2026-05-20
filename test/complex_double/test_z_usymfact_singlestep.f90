#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_usymfact_singlestep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine tests one Francis single-shift step for a unitary upper
! Hessenberg matrix stored in symmetric factored form.
!
! Tested subroutine interface:
!
!     call z_usymfact_singlestep(NORMALIZE,VEC,N,U,V,NU,M,Z,ITCNT)
!
! Check 1:
!
!     N  = 256,
!     ui = 0, vi = 1,  i = 1,...,N-1,
!     uN = 1, vN = 0,
!     nu = 1.
!
! Check 2:
!
!     N = 32,
!     random normalized unitary factors.
!
! With VEC = .TRUE., M = N, and Z initialized to the identity, the routine
! accumulates the unitary similarity transformation used in the QR step.
! Therefore the output matrix should satisfy
!
!     A1 = Z^* A0 Z.
!
! The characteristic-polynomial check has been dropped.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_usymfact_singlestep

  implicit none

  ! parameters
  integer, parameter :: N = 256
  integer, parameter :: M = N
  integer, parameter :: NR = 32
  integer, parameter :: MR = NR

  real(8), parameter :: eps = (EISCOR_DBL_EPS)
  real(8), parameter :: tol = 100d0*(EISCOR_DBL_EPS)

  ! deterministic N = 256 test variables
  complex(8) :: U(N), U0(N)
  real(8) :: V(N), V0(N)
  complex(8) :: NU
  integer :: ITCNT
  complex(8) :: Z(M,N)
  complex(8) :: A0(N,N), A1(N,N)

  ! random N = 32 test variables
  complex(8) :: UR(NR), UR0(NR)
  real(8) :: VR(NR), VR0(NR)
  complex(8) :: NUR
  integer :: ITCNTR
  complex(8) :: ZR(MR,NR)
  complex(8) :: AR0(NR,NR), AR1(NR,NR)

  ! normalization and eigenvector/similarity flags
  logical, parameter :: NORMALIZE(3) = (/ .TRUE., .TRUE., .TRUE. /)
  logical :: VEC

  ! timing variables
  integer :: c_start, c_stop, c_rate

  ! loop index
  integer :: ii

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! check 1)
  !
  ! N = 256,
  !
  !     ui = 0, vi = 1,  i = 1,...,N-1,
  !     uN = 1, vN = 0,
  !     nu = 1.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! initialize factors
  NU = cmplx(1d0,0d0,kind=8)

  do ii = 1,N-1
    U(ii) = cmplx(0d0,0d0,kind=8)
    V(ii) = 1d0
  end do

  U(N) = cmplx(1d0,0d0,kind=8)
  V(N) = 0d0

  ITCNT = 0

  ! save input factors
  U0 = U
  V0 = V

  ! explicitly form the input matrix
  call form_usymfact_matrix(N,U0,V0,NU,A0)

  ! initialize eigenvector/similarity accumulator
  VEC = .TRUE.
  call set_identity(M,N,Z)

  ! check input normalization
  call check_factor_normalization(N,U,V,NU,tol,__LINE__)

  ! check input unitarity
  call check_unitary_matrix(N,A0,tol,__LINE__)

  ! perform one QR step and accumulate the similarity transformation in Z
  call z_usymfact_singlestep(NORMALIZE,VEC,N,U,V,NU,M,Z,ITCNT)

  ! explicitly form the output matrix
  call form_usymfact_matrix(N,U,V,NU,A1)

  ! check output normalization
  call check_factor_normalization(N,U,V,NU,tol,__LINE__)

  ! check output unitarity
  call check_unitary_matrix(N,A1,tol,__LINE__)

  ! check that Z is unitary
  call check_unitary_matrix(N,Z,tol,__LINE__)

  ! check the explicit similarity relation
  !
  !     A1 = Z^* A0 Z.
  !
  call check_similarity_matrix(N,A0,A1,Z,tol,__LINE__)




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! check 2)
  !
  ! Random normalized unitary factors, N = 32.
  !
  !     |u_i|^2 + v_i^2 = 1,  i = 1,...,NR-1,
  !     |u_NR| = 1,
  !     v_NR = 0,
  !     |nu| = 1.
  !
  ! We initialize ZR = I and check
  !
  !     AR1 = ZR^* AR0 ZR.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! deterministic random seed
  call set_random_seed()

  ! initialize random normalized factors
  call init_random_usymfact(NR,UR,VR,NUR)

  ITCNTR = 0

  ! save input factors
  UR0 = UR
  VR0 = VR

  ! explicitly form the input matrix
  call form_usymfact_matrix(NR,UR0,VR0,NUR,AR0)

  ! initialize eigenvector/similarity accumulator
  VEC = .TRUE.
  call set_identity(MR,NR,ZR)

  ! check input normalization
  call check_factor_normalization(NR,UR,VR,NUR,tol,__LINE__)

  ! check input unitarity
  call check_unitary_matrix(NR,AR0,tol,__LINE__)

  ! perform one QR step and accumulate the similarity transformation in ZR
  call z_usymfact_singlestep(NORMALIZE,VEC,NR,UR,VR,NUR,MR,ZR,ITCNTR)

  ! explicitly form the output matrix
  call form_usymfact_matrix(NR,UR,VR,NUR,AR1)

  ! check output normalization
  call check_factor_normalization(NR,UR,VR,NUR,tol,__LINE__)

  ! check output unitarity
  call check_unitary_matrix(NR,AR1,tol,__LINE__)

  ! check that ZR is unitary
  call check_unitary_matrix(NR,ZR,tol,__LINE__)

  ! check the explicit similarity relation
  !
  !     AR1 = ZR^* AR0 ZR.
  !
  call check_similarity_matrix(NR,AR0,AR1,ZR,tol,__LINE__)




  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! set_identity
  !
  ! Sets Z = I_N inside an M x N array.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_identity(M,N,Z)

    implicit none

    integer, intent(in) :: M, N
    complex(8), intent(out) :: Z(M,N)

    complex(8) :: zero, one
    integer :: ii

    zero = cmplx(0d0,0d0,kind=8)
    one  = cmplx(1d0,0d0,kind=8)

    Z = zero

    do ii = 1,N
      Z(ii,ii) = one
    end do

  end subroutine set_identity




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! set_random_seed
  !
  ! Sets a deterministic random seed so that the randomized test is reproducible.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine set_random_seed()

    implicit none

    integer :: nseed, ii
    integer, allocatable :: seed(:)

    call random_seed(size=nseed)
    allocate(seed(nseed))

    do ii = 1,nseed
      seed(ii) = 7919 + 104729*ii
    end do

    call random_seed(put=seed)

    deallocate(seed)

  end subroutine set_random_seed




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! init_random_usymfact
  !
  ! Randomly initializes a normalized unitary upper-Hessenberg factorization:
  !
  !     |U(i)|^2 + V(i)^2 = 1,  i = 1,...,N-1,
  !     |U(N)| = 1,
  !     V(N) = 0,
  !     |NU| = 1.
  !
  ! The variables V(i) are nonnegative.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_random_usymfact(N,U,V,NU)

    implicit none

    integer, intent(in) :: N
    complex(8), intent(out) :: U(N)
    real(8), intent(out) :: V(N)
    complex(8), intent(out) :: NU

    real(8) :: r1, r2, theta, phi
    real(8) :: pi
    integer :: ii

    pi = 4d0*atan(1d0)

    ! leading phase NU
    call random_number(r1)
    phi = 2d0*pi*r1
    NU = cmplx(cos(phi),sin(phi),kind=8)

    ! proper core transformations
    do ii = 1,N-1

      call random_number(r1)
      call random_number(r2)

      ! theta in [0,pi/2], so V(ii) >= 0
      theta = 5d-1*pi*r1
      phi   = 2d0*pi*r2

      U(ii) = cos(theta)*cmplx(cos(phi),sin(phi),kind=8)
      V(ii) = sin(theta)

    end do

    ! trailing diagonal phase
    call random_number(r1)
    phi = 2d0*pi*r1

    U(N) = cmplx(cos(phi),sin(phi),kind=8)
    V(N) = 0d0

  end subroutine init_random_usymfact




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! form_usymfact_matrix
  !
  ! Forms the dense matrix represented by the symmetric unitary factorization
  !
  !     A = diag(NU,1,...,1) * G1 * G2 * ... * G_{N-1}
  !         * diag(1,...,1,U(N)),
  !
  ! where
  !
  !     G_i(i:i+1,i:i+1) =
  !
  !       [ U(i)       -V(i)        ]
  !       [ V(i)   conjg(U(i))      ].
  !
  ! This implementation applies each core transformation by updating the
  ! affected columns of A. It avoids constructing a dense G and calling
  ! MATMUL inside the loop.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine form_usymfact_matrix(N,U,V,NU,A)

    implicit none

    integer, intent(in) :: N
    complex(8), intent(in) :: U(N)
    real(8), intent(in) :: V(N)
    complex(8), intent(in) :: NU
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

    ! left phase diag(NU,1,...,1)
    A(1,1) = NU

    ! multiply on the right by G_1, ..., G_{N-1}
    !
    ! Columns kk and kk+1 are updated by
    !
    !     [ col_kk  col_kk+1 ] <- [ col_kk  col_kk+1 ]
    !
    !       [ U(kk)       -V(kk)       ]
    !       [ V(kk)    conjg(U(kk))   ].
    !
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
  ! check_factor_normalization
  !
  ! Checks
  !
  !     |U(i)|^2 + V(i)^2 = 1,   i = 1,...,N-1,
  !     V(N) = 0,
  !     |U(N)| = 1,
  !     |NU| = 1.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_factor_normalization(N,U,V,NU,tol,line)

    implicit none

    integer, intent(in) :: N
    complex(8), intent(in) :: U(N)
    real(8), intent(in) :: V(N)
    complex(8), intent(in) :: NU
    real(8), intent(in) :: tol
    integer, intent(in) :: line

    real(8) :: err
    integer :: ii

    err = abs(dble(NU)**2 + aimag(NU)**2 - 1d0)

    do ii = 1,N-1
      err = max(err,abs(dble(U(ii))**2 + aimag(U(ii))**2 + V(ii)**2 - 1d0))
    end do

    err = max(err,abs(V(N)))
    err = max(err,abs(dble(U(N))**2 + aimag(U(N))**2 - 1d0))

    if (err > tol) then
       call u_test_failed(line)
    end if

  end subroutine check_factor_normalization




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! check_unitary_matrix
  !
  ! Checks
  !
  !     A^* A = I
  !
  ! entrywise.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_unitary_matrix(N,A,tol,line)

    implicit none

    integer, intent(in) :: N
    complex(8), intent(in) :: A(N,N)
    real(8), intent(in) :: tol
    integer, intent(in) :: line

    complex(8) :: B(N,N)
    complex(8) :: zero, one
    real(8) :: err
    integer :: ii, jj

    zero = cmplx(0d0,0d0,kind=8)
    one  = cmplx(1d0,0d0,kind=8)

    B = matmul(conjg(transpose(A)),A)

    err = 0d0
    do jj = 1,N
      do ii = 1,N
        if (ii == jj) then
          err = max(err,abs(B(ii,jj)-one))
        else
          err = max(err,abs(B(ii,jj)-zero))
        end if
      end do
    end do

    if (err > tol) then
       call u_test_failed(line)
    end if

  end subroutine check_unitary_matrix




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! check_similarity_matrix
  !
  ! Checks
  !
  !     A1 = Z^* A0 Z.
  !
  ! This is the dense form of the accumulated QR-step similarity.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_similarity_matrix(N,A0,A1,Z,tol,line)

    implicit none

    integer, intent(in) :: N
    complex(8), intent(in) :: A0(N,N)
    complex(8), intent(in) :: A1(N,N)
    complex(8), intent(in) :: Z(N,N)
    real(8), intent(in) :: tol
    integer, intent(in) :: line

    complex(8) :: AZ(N,N), S(N,N)
    real(8) :: err
    integer :: ii, jj

    AZ = matmul(A0,Z)
    S  = matmul(conjg(transpose(Z)),AZ)

    err = 0d0
    do jj = 1,N
      do ii = 1,N
        err = max(err,abs(A1(ii,jj)-S(ii,jj)))
      end do
    end do

    if (err > tol) then
       call u_test_failed(line)
    end if

  end subroutine check_similarity_matrix

end program test_z_usymfact_singlestep
