#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_rot3_symturnover
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks limiting and near-unitary cases in the standard
! symmetric turnover.
!
! Tested subroutine interface:
!
!     call z_rot3_symturnover(GAMMA,SIGMA,C,S,U,V,RHO)
!
! where
!
!     GAMMA, SIGMA, U, RHO are complex,
!     C, S, V are real,
!
! and the active phase product is
!
!     GAMMA*SIGMA**2.
!
! For every scalar-output test, this program also explicitly forms
!
!     H    = R * tildeQ * U * Q
!     Hhat = Qhat * Uhat * Rhat * tildeQhat
!
! from the input and computed output parameters and checks that
! max(abs(H-Hhat)) is small.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_rot3_symturnover

  implicit none

  ! compute variables
  real(8), parameter :: eps = (EISCOR_DBL_EPS)
  real(8), parameter :: tol = 100d0*(EISCOR_DBL_EPS)
  real(8) :: nrm, C, S, V
  complex(8) :: GAMMA, SIGMA, U, RHO

  real(8) :: C_in, S_in, V_in
  complex(8) :: GAMMA_in, SIGMA_in, U_in

  real(8) :: C_out, S_out, V_out
  complex(8) :: GAMMA_out, SIGMA_out, U_out

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
  ! Non-trivial, no deflation.

  ! known inputs
  GAMMA_in = cmplx(7.071067811865476d-01, 7.071067811865476d-01, kind=8)
  SIGMA_in = cmplx(7.071067811865475d-01, 7.071067811865475d-01, kind=8)
  C_in     =       8.164965809277260d-01
  S_in     =       5.773502691896257d-01
  U_in     = cmplx(5.773502691896257d-01, 5.773502691896257d-01, kind=8)
  V_in     =       5.773502691896257d-01
  RHO      = cmplx(8.944271909999159d-01, 4.472135954999579d-01, kind=8)

  ! known outputs
  GAMMA_out = cmplx(-9.486832980505138d-01, -3.162277660168379d-01, kind=8)
  SIGMA_out = cmplx(-8.449489742783178d-01, -5.348469228349536d-01, kind=8)
  C_out     =        8.770580193070294d-01
  S_out     =        4.803844614152615d-01
  U_out     = cmplx( 4.690174004250532d-01, -5.463892354512888d-01, kind=8)
  V_out     =        6.938886664887108d-01

  ! store inputs for comparison
  GAMMA = GAMMA_in
  SIGMA = SIGMA_in
  C = C_in
  S = S_in
  U = U_in
  V = V_in

  ! perform symmetric turnover
  call z_rot3_symturnover(GAMMA,SIGMA,C,S,U,V,RHO)

  ! check explicit 3x3 products
  call check_symturnover_product( &
       GAMMA_in,SIGMA_in,C_in,S_in,U_in,V_in,RHO, &
       GAMMA,SIGMA,C,S,U,V, &
       tol,__LINE__)

  ! check error against known output
  nrm = abs(GAMMA-GAMMA_out) + abs(SIGMA-SIGMA_out) &
      + abs(C-C_out) + abs(S-S_out) &
      + abs(U-U_out) + abs(V-V_out)

  if (nrm > tol) then
     call u_test_failed(__LINE__)
  end if




  !!!!!!!!!!!!!!!!!!!!
  ! check 2)
  !
  ! Diagonal example. Here C = 0, so SIGMA is arbitrary on input.
  ! We choose SIGMA = 1.

  ! known inputs
  GAMMA_in = cmplx(6.000000000000000d-01, -8.000000000000000d-01, kind=8)
  SIGMA_in = cmplx(1.000000000000000d+00,  0.000000000000000d+00, kind=8)
  C_in     =       0.000000000000000d+00
  S_in     =       1.000000000000000d+00
  U_in     = cmplx(4.472135954999579d-01, -8.944271909999159d-01, kind=8)
  V_in     =       0.000000000000000d+00
  RHO      = cmplx(8.944271909999159d-01,  4.472135954999579d-01, kind=8)

  ! known outputs
  GAMMA_out = cmplx(-1.788854381999831d-01,  9.838699100999074d-01, kind=8)
  SIGMA_out = cmplx( 1.000000000000000d+00,  0.000000000000000d+00, kind=8)
  C_out     =        0.000000000000000d+00
  S_out     =        1.000000000000000d+00
  U_out     = cmplx( 0.000000000000000d+00, -9.999999999999999d-01, kind=8)
  V_out     =        0.000000000000000d+00

  ! store inputs for comparison
  GAMMA = GAMMA_in
  SIGMA = SIGMA_in
  C = C_in
  S = S_in
  U = U_in
  V = V_in

  ! perform symmetric turnover
  call z_rot3_symturnover(GAMMA,SIGMA,C,S,U,V,RHO)

  ! check explicit 3x3 products
  call check_symturnover_product( &
       GAMMA_in,SIGMA_in,C_in,S_in,U_in,V_in,RHO, &
       GAMMA,SIGMA,C,S,U,V, &
       tol,__LINE__)

  ! check error against known output
  nrm = abs(GAMMA-GAMMA_out) + abs(SIGMA-SIGMA_out) &
      + abs(C-C_out) + abs(S-S_out) &
      + abs(U-U_out) + abs(V-V_out)

  if (nrm > tol) then
     call u_test_failed(__LINE__)
  end if




  !!!!!!!!!!!!!!!!!!!!
  ! check 3)
  !
  ! Non-trivial with deflation / singular z branch.
  ! The output has C = 0, so SIGMA is not uniquely determined.
  ! The implementation chooses SIGMA = 1.

  ! known inputs
  GAMMA_in = cmplx(1.000000000000000d+00,  0.000000000000000d+00, kind=8)
  SIGMA_in = cmplx(0.000000000000000d+00, -1.000000000000000d+00, kind=8)
  C_in     =       7.071067811865476d-01
  S_in     =       7.071067811865476d-01
  U_in     = cmplx(1.000000000000000d+00,  0.000000000000000d+00, kind=8)
  V_in     =       0.000000000000000d+00
  RHO      = cmplx(0.000000000000000d+00,  1.000000000000000d+00, kind=8)

  ! known outputs
  GAMMA_out = cmplx(0.000000000000000d+00,  1.000000000000000d+00, kind=8)
  SIGMA_out = cmplx(1.000000000000000d+00,  0.000000000000000d+00, kind=8)
  C_out     =       0.000000000000000d+00
  S_out     =       1.000000000000000d+00
  U_out     = cmplx(0.000000000000000d+00, -1.000000000000000d+00, kind=8)
  V_out     =       0.000000000000000d+00

  ! store inputs for comparison
  GAMMA = GAMMA_in
  SIGMA = SIGMA_in
  C = C_in
  S = S_in
  U = U_in
  V = V_in

  ! perform symmetric turnover
  call z_rot3_symturnover(GAMMA,SIGMA,C,S,U,V,RHO)

  ! check explicit 3x3 products
  call check_symturnover_product( &
       GAMMA_in,SIGMA_in,C_in,S_in,U_in,V_in,RHO, &
       GAMMA,SIGMA,C,S,U,V, &
       tol,__LINE__)

  ! check error against known output
  nrm = abs(GAMMA-GAMMA_out) + abs(SIGMA-SIGMA_out) &
      + abs(C-C_out) + abs(S-S_out) &
      + abs(U-U_out) + abs(V-V_out)

  if (nrm > tol) then
     call u_test_failed(__LINE__)
  end if




  !!!!!!!!!!!!!!!!!!!!
  ! check 4)
  !
  ! U = 1, V = sqrt(eps), SIGMA = 1, GAMMA = 1, RHO = -1,
  ! C = S = 1/sqrt(2).
  ! This tests the Re(T) >= 0 branch near tiny V.

  ! known inputs
  GAMMA_in = cmplx(1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  SIGMA_in = cmplx(1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  C_in     =       7.071067811865476d-01
  S_in     =       7.071067811865476d-01
  U_in     = cmplx(1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  V_in     =       1.490116119384766d-08
  RHO      = cmplx(-1.000000000000000d+00, 0.000000000000000d+00, kind=8)

  ! known outputs
  GAMMA_out = cmplx(1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  SIGMA_out = cmplx(1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  C_out     =       1.000000000000000d+00
  S_out     =       1.053671212772351d-08
  U_out     = cmplx(0.000000000000000d+00, 0.000000000000000d+00, kind=8)
  V_out     =       1.000000000000000d+00

  ! store inputs for comparison
  GAMMA = GAMMA_in
  SIGMA = SIGMA_in
  C = C_in
  S = S_in
  U = U_in
  V = V_in

  ! perform symmetric turnover
  call z_rot3_symturnover(GAMMA,SIGMA,C,S,U,V,RHO)

  ! check explicit 3x3 products
  call check_symturnover_product( &
       GAMMA_in,SIGMA_in,C_in,S_in,U_in,V_in,RHO, &
       GAMMA,SIGMA,C,S,U,V, &
       tol,__LINE__)

  ! check error against known output
  nrm = abs(GAMMA-GAMMA_out) + abs(SIGMA-SIGMA_out) &
      + abs(C-C_out) + abs(S-S_out) &
      + abs(U-U_out) + abs(V-V_out)

  if (nrm > tol) then
     call u_test_failed(__LINE__)
  end if




  !!!!!!!!!!!!!!!!!!!!
  ! check 5)
  !
  ! Random nearly unitary inputs with S = 0.
  ! The constants below are deterministic random-looking values.

  ! known inputs
  GAMMA_in = cmplx( 9.930091528749599d-01, -1.180373767353156d-01, kind=8)
  SIGMA_in = cmplx(-9.930213065497917d-01,  1.179350869679778d-01, kind=8)
  C_in     =        1.000000000000001d+00
  S_in     =        0.000000000000000d+00
  U_in     = cmplx(-9.213754130678973d-01,  3.886738326591608d-01, kind=8)
  V_in     =        3.332000937312528d-08
  RHO      = cmplx( 9.786370151869102d-01, -2.055957015748513d-01, kind=8)

  ! known outputs
  GAMMA_out = cmplx(-9.960634907047851d-01, -8.864266740683405d-02, kind=8)
  SIGMA_out = cmplx( 3.405592881001832d-01,  9.402230433725267d-01, kind=8)
  C_out     =        9.999999999997192d-01
  S_out     =        7.495297542322489d-07
  U_out     = cmplx(-9.891178742925306d-01,  1.471252213423164d-01, kind=8)
  V_out     =        0.000000000000000d+00

  ! store inputs for comparison
  GAMMA = GAMMA_in
  SIGMA = SIGMA_in
  C = C_in
  S = S_in
  U = U_in
  V = V_in

  ! perform symmetric turnover
  call z_rot3_symturnover(GAMMA,SIGMA,C,S,U,V,RHO)

  ! check explicit 3x3 products
  call check_symturnover_product( &
       GAMMA_in,SIGMA_in,C_in,S_in,U_in,V_in,RHO, &
       GAMMA,SIGMA,C,S,U,V, &
       tol,__LINE__)

  ! check error against known output
  nrm = abs(GAMMA-GAMMA_out) + abs(SIGMA-SIGMA_out) &
      + abs(C-C_out) + abs(S-S_out) &
      + abs(U-U_out) + abs(V-V_out)

  if (nrm > tol) then
     call u_test_failed(__LINE__)
  end if




  !!!!!!!!!!!!!!!!!!!!
  ! check 6)
  !
  ! U = -1, V = sqrt(eps), SIGMA = 1, GAMMA = 1, RHO = -1,
  ! C = S = 1/sqrt(2).
  ! This tests the Re(T) < 0 branch near cancellation.

  ! known inputs
  GAMMA_in = cmplx(1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  SIGMA_in = cmplx(1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  C_in     =       7.071067811865476d-01
  S_in     =       7.071067811865476d-01
  U_in     = cmplx(-1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  V_in     =       1.490116119384766d-08
  RHO      = cmplx(-1.000000000000000d+00, 0.000000000000000d+00, kind=8)

  ! known outputs
  GAMMA_out = cmplx(1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  SIGMA_out = cmplx(1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  C_out     =       5.268356063861754d-09
  S_out     =       1.000000000000000d+00
  U_out     = cmplx(1.000000000000000d+00, 0.000000000000000d+00, kind=8)
  V_out     =       1.053671212772351d-08

  ! store inputs for comparison
  GAMMA = GAMMA_in
  SIGMA = SIGMA_in
  C = C_in
  S = S_in
  U = U_in
  V = V_in

  ! perform symmetric turnover
  call z_rot3_symturnover(GAMMA,SIGMA,C,S,U,V,RHO)

  ! check explicit 3x3 products
  call check_symturnover_product( &
       GAMMA_in,SIGMA_in,C_in,S_in,U_in,V_in,RHO, &
       GAMMA,SIGMA,C,S,U,V, &
       tol,__LINE__)

  ! check error against known output
  nrm = abs(GAMMA-GAMMA_out) + abs(SIGMA-SIGMA_out) &
      + abs(C-C_out) + abs(S-S_out) &
      + abs(U-U_out) + abs(V-V_out)

  if (nrm > tol) then
     call u_test_failed(__LINE__)
  end if




  ! stop timer
  call system_clock(count=c_stop)

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! check_symturnover_product
  !
  ! Explicitly forms
  !
  !     H    = R * tildeQ * U * Q
  !     Hhat = Qhat * Uhat * Rhat * tildeQhat
  !
  ! and checks max(abs(H-Hhat)) <= tol.
  !
  ! Pre-turnover variables:
  !
  !     GAMMA0, SIGMA0, C0, S0, U0, V0, RHO
  !
  ! Post-turnover variables:
  !
  !     GAMMA1, SIGMA1, C1, S1, U1, V1
  !
  ! The convention used here is
  !
  !     Q =
  !       [ C*SIGMA       -S              ]
  !       [ S              C*conjg(SIGMA) ]
  !
  !     tildeQ =
  !       [ C*GAMMA*SIGMA        -S                         ]
  !       [ S                     C*conjg(GAMMA*SIGMA)       ]
  !
  ! with the obvious 3 x 3 embeddings.
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine check_symturnover_product( &
       GAMMA0,SIGMA0,C0,S0,U0,V0,RHO, &
       GAMMA1,SIGMA1,C1,S1,U1,V1, &
       tol,line)

    implicit none

    ! input variables
    real(8), intent(in) :: C0, S0, V0
    real(8), intent(in) :: C1, S1, V1
    real(8), intent(in) :: tol
    integer, intent(in) :: line

    complex(8), intent(in) :: GAMMA0, SIGMA0, U0
    complex(8), intent(in) :: GAMMA1, SIGMA1, U1
    complex(8), intent(in) :: RHO

    ! local matrices
    complex(8) :: R(3,3), TQ(3,3), UM(3,3), Q(3,3), H(3,3)
    complex(8) :: QH(3,3), UH(3,3), RH(3,3), TQH(3,3), HH(3,3)
    complex(8) :: zero, one
    complex(8) :: gs0, gs1
    real(8) :: herr
    integer :: ii, jj

    zero = cmplx(0d0,0d0,kind=8)
    one  = cmplx(1d0,0d0,kind=8)

    gs0 = GAMMA0*SIGMA0
    gs1 = GAMMA1*SIGMA1

    ! initialize matrices to zero
    R   = zero
    TQ  = zero
    UM  = zero
    Q   = zero
    H   = zero

    QH  = zero
    UH  = zero
    RH  = zero
    TQH = zero
    HH  = zero

    ! -------------------------------------------------------------------------
    ! R_i, acting on rows/columns 1 and 2:
    !
    !     [ -conjg(rho)       0      ]
    !     [      0          -rho     ]
    !
    ! embedded in 3 x 3.
    ! -------------------------------------------------------------------------
    R(1,1) = -conjg(RHO)
    R(2,2) = -RHO
    R(3,3) = one

    ! -------------------------------------------------------------------------
    ! tilde Q_i, acting on rows/columns 1 and 2:
    !
    !     [ c*gamma*sigma          -s                    ]
    !     [ s                       c*conjg(gamma*sigma) ]
    ! -------------------------------------------------------------------------
    TQ(1,1) = C0*gs0
    TQ(1,2) = -S0
    TQ(2,1) = S0
    TQ(2,2) = C0*conjg(gs0)
    TQ(3,3) = one

    ! -------------------------------------------------------------------------
    ! U_{i+1}, acting on rows/columns 2 and 3:
    !
    !     [ 1       0          0        ]
    !     [ 0       u         -v        ]
    !     [ 0       v       conjg(u)    ]
    ! -------------------------------------------------------------------------
    UM(1,1) = one
    UM(2,2) = U0
    UM(2,3) = -V0
    UM(3,2) = V0
    UM(3,3) = conjg(U0)

    ! -------------------------------------------------------------------------
    ! Q_i, acting on rows/columns 1 and 2:
    !
    !     [ c*sigma          -s              ]
    !     [ s                 c*conjg(sigma)]
    ! -------------------------------------------------------------------------
    Q(1,1) = C0*SIGMA0
    Q(1,2) = -S0
    Q(2,1) = S0
    Q(2,2) = C0*conjg(SIGMA0)
    Q(3,3) = one

    ! H = R * tildeQ * U * Q
    H = matmul(R,matmul(TQ,matmul(UM,Q)))

    ! -------------------------------------------------------------------------
    ! Q_{i+1}, acting on rows/columns 2 and 3:
    !
    !     [ 1       0              0              ]
    !     [ 0       c*sigma       -s             ]
    !     [ 0       s              c*conjg(sigma)]
    ! -------------------------------------------------------------------------
    QH(1,1) = one
    QH(2,2) = C1*SIGMA1
    QH(2,3) = -S1
    QH(3,2) = S1
    QH(3,3) = C1*conjg(SIGMA1)

    ! -------------------------------------------------------------------------
    ! Uhat_i, acting on rows/columns 1 and 2:
    !
    !     [ uhat      -vhat        0 ]
    !     [ vhat    conjg(uhat)    0 ]
    !     [ 0          0           1 ]
    ! -------------------------------------------------------------------------
    UH(1,1) = U1
    UH(1,2) = -V1
    UH(2,1) = V1
    UH(2,2) = conjg(U1)
    UH(3,3) = one

    ! -------------------------------------------------------------------------
    ! R_{i+1}, acting on rows/columns 2 and 3:
    !
    !     [ 1        0              0   ]
    !     [ 0      -conjg(rho)      0   ]
    !     [ 0        0            -rho  ]
    ! -------------------------------------------------------------------------
    RH(1,1) = one
    RH(2,2) = -conjg(RHO)
    RH(3,3) = -RHO

    ! -------------------------------------------------------------------------
    ! tilde Q_{i+1}, acting on rows/columns 2 and 3:
    !
    !     [ 1       0                       0                   ]
    !     [ 0       c*gamma*sigma          -s                   ]
    !     [ 0       s                       c*conjg(gamma*sigma)]
    ! -------------------------------------------------------------------------
    TQH(1,1) = one
    TQH(2,2) = C1*gs1
    TQH(2,3) = -S1
    TQH(3,2) = S1
    TQH(3,3) = C1*conjg(gs1)

    ! Hhat = Qhat * Uhat * Rhat * tildeQhat
    HH = matmul(QH,matmul(UH,matmul(RH,TQH)))

    ! entrywise maximum absolute error
    herr = 0d0
    do jj = 1,3
      do ii = 1,3
        herr = max(herr,abs(H(ii,jj)-HH(ii,jj)))
      end do
    end do

    if (herr > tol) then
       call u_test_failed(line)
    end if

  end subroutine check_symturnover_product

end program test_z_rot3_symturnover
