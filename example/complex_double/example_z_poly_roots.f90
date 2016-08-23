#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_z_poly_roots
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program runs the examples from the paper
! Jared L. Aurentz, Thomas Mach, Raf Vandebril, and David S. Watkins. 
! [_Fast and backward stable computation of roots of polynomials._]
! (http://epubs.siam.org/doi/10.1137/140983434) 
! SIAM Journal on Matrix Analysis and Applications. Vol. 36, No. 3, 
! pp. 942-973. 2015.
!
! 1) 48 special polynomials, backward errors and forward errors as
!    in Table 8 and 10
!
!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_z_poly_roots

  implicit none
  
  ! compute variables
  real(8), parameter :: pi = 3.14159265358979323846264338327950d0

  integer :: N 
  real(8) :: res, forw, err
  integer :: ii, jj, kk, ll, INFO
  integer, allocatable :: ARGS(:)
  real(8) :: TABLE_F(48), TABLE_B(48), normofp, lambda
  real(8), allocatable :: RESIDUALS(:), FORWARD(:)
  complex(8), allocatable :: COEFFS(:), EROOTS(:), ROOTS(:)
  logical :: eigsknown, EK(48)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print banner
  print*,""
  print*,"example_z_poly_roots"
  print*,""

  TABLE_B = -1d0
  TABLE_F = -1d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 1) loop overall 48 special polynomials
  do kk = 1,48 
     select case (kk)
     case (1) ! Wilkinson Polynomial 10
        N = 10
     case (2) ! Wilkinson Polynomial 15
        N = 15
     case (3) ! Wilkinson Polynomial 20
        N = 20
     case (4)  ! monic polynomial roots [-2.1:0.2:1.7]
        N = 20
     case (5) ! inverse Wilkinson Polynomial 10
        N = 10
     case (6) ! inverse Wilkinson Polynomial 15
        N = 15
     case (7) ! inverse Wilkinson Polynomial 20
        N = 20
     case (8)  ! univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9
        N = 20   
     case (9)  ! univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9
        N = 20
     case (10) ! Chebyshev polynomial of degree 20 
        N = 20

     case (11) ! z^20 + z^19 + ... + z + 1
        N = 20
     case (12) ! C. Traverso 24 MPSolve"
        N = 24
     case (13) ! Mandelbort 31 MPSolve
        N = 31
     case (14) ! Mandelbort 63 MPSolve
        N = 63
     case (15) ! Vanni Noferini's example
        N = 12
     case (16) ! Vanni Noferini's example
        N = 35
     case (17) ! cubic rcoeffsnomial small roots
        N = 3
     case (18) ! cubic rcoeffsnomial very small roots
        N = 3
     case (19) ! cubic rcoeffsnomial large roots
        N = 3
     case (20) ! cubic rcoeffsnomial very large roots
        N = 3

     case (21) ! underflow test
        N = 10
     case (22) ! underflow test
        N = 20
     case (23) ! deflation stability test
        N = 3
     case (24) ! deflation stability test
        N = 3
     case (25) ! deflation stability test
        N = 3
     case (26) ! deflation stability test M = 15
        N = 60
     case (27)  ! Bernoulli rcoeffsnomial of degree 20  
        N = 20
     case (28)   ! p(z) = (20!) sum_{k=0}^{20} z^k/k!
        N = 20

     case (29) ! Bevilacqua, Del Corso, Gemignani P1
        N = 20
     case (30) ! Bevilacqua, Del Corso, Gemignani P1
        N = 40
     case (31) ! Bevilacqua, Del Corso, Gemignani P1
        N = 60
     case (32) ! Bevilacqua, Del Corso, Gemignani P1
        N = 512
     case (33) ! Bevilacqua, Del Corso, Gemignani P1
        N = 1024

     case (34) ! Bevilacqua, Del Corso, Gemignani P2
        N = 20
     case (35) ! Bevilacqua, Del Corso, Gemignani P2
        N = 40
     case (36) ! Bevilacqua, Del Corso, Gemignani P2
        N = 60
     case (37) ! Bevilacqua, Del Corso, Gemignani P2
        N = 512
     case (38) ! Bevilacqua, Del Corso, Gemignani P2
        N = 1024

     case (39) ! Bevilacqua, Del Corso, Gemignani P3, lambda = 0.9
        N = 20
     case (40) ! Bevilacqua, Del Corso, Gemignani P3, lambda = 0.9
        N = 40
     case (41) ! Bevilacqua, Del Corso, Gemignani P3, lambda = 0.9
        N = 60
     case (42) ! Bevilacqua, Del Corso, Gemignani P3, lambda = 0.9
        N = 512
     case (43) ! Bevilacqua, Del Corso, Gemignani P3, lambda = 0.9
        N = 1024

     case (44) ! Bevilacqua, Del Corso, Gemignani P3, lambda = 0.999
        N = 20
     case (45) ! Bevilacqua, Del Corso, Gemignani P3, lambda = 0.999
        N = 40
     case (46) ! Bevilacqua, Del Corso, Gemignani P3, lambda = 0.999
        N = 60
     case (47) ! Bevilacqua, Del Corso, Gemignani P3, lambda = 0.999
        N = 512
     case (48) ! Bevilacqua, Del Corso, Gemignani P3, lambda = 0.999
        N = 1024

     end select

     allocate(RESIDUALS(N),FORWARD(N),COEFFS(N+1),EROOTS(N),ROOTS(N),ARGS(N))

     eigsknown = .FALSE.
     select case (kk)
     case (1) ! Wilkinson Polynomial !check
        print*, "Wilkinson Polynomial degree 10"
        eigsknown = .TRUE.
        EROOTS(1)  = cmplx(1d0,0d0,kind=8)
        EROOTS(2)  = cmplx(2d0,0d0,kind=8)
        EROOTS(3)  = cmplx(3d0,0d0,kind=8)
        EROOTS(4)  = cmplx(4d0,0d0,kind=8)
        EROOTS(5)  = cmplx(5d0,0d0,kind=8)
        EROOTS(6)  = cmplx(6d0,0d0,kind=8)
        EROOTS(7)  = cmplx(7d0,0d0,kind=8)
        EROOTS(8)  = cmplx(8d0,0d0,kind=8)
        EROOTS(9)  = cmplx(9d0,0d0,kind=8)
        EROOTS(10) = cmplx(10d0,0d0,kind=8)
        COEFFS(  1) = cmplx( 1.000000000000000000d+00, 0.000000000000000000d+00,kind=8)
        COEFFS(  2) = cmplx(-5.500000000000000000d+01, 0.000000000000000000d+00,kind=8)
        COEFFS(  3) = cmplx( 1.320000000000000000d+03, 0.000000000000000000d+00,kind=8)
        COEFFS(  4) = cmplx(-1.815000000000000000d+04, 0.000000000000000000d+00,kind=8)
        COEFFS(  5) = cmplx( 1.577730000000000000d+05, 0.000000000000000000d+00,kind=8)
        COEFFS(  6) = cmplx(-9.020550000000000000d+05, 0.000000000000000000d+00,kind=8)
        COEFFS(  7) = cmplx( 3.416930000000000000d+06, 0.000000000000000000d+00,kind=8)
        COEFFS(  8) = cmplx(-8.409500000000000000d+06, 0.000000000000000000d+00,kind=8)
        COEFFS(  9) = cmplx( 1.275357600000000000d+07, 0.000000000000000000d+00,kind=8)
        COEFFS( 10) = cmplx(-1.062864000000000000d+07, 0.000000000000000000d+00,kind=8)
        COEFFS( 11) = cmplx( 3.628800000000000000d+06, 0.000000000000000000d+00,kind=8)
        
     case (2) ! Wilkinson Polynomial !check
        print*, "Wilkinson Polynomial degree 15"
        eigsknown = .TRUE.
        EROOTS(1)  = cmplx(1d0,0d0,kind=8)
        EROOTS(2)  = cmplx(2d0,0d0,kind=8)
        EROOTS(3)  = cmplx(3d0,0d0,kind=8)
        EROOTS(4)  = cmplx(4d0,0d0,kind=8)
        EROOTS(5)  = cmplx(5d0,0d0,kind=8)
        EROOTS(6)  = cmplx(6d0,0d0,kind=8)
        EROOTS(7)  = cmplx(7d0,0d0,kind=8)
        EROOTS(8)  = cmplx(8d0,0d0,kind=8)
        EROOTS(9)  = cmplx(9d0,0d0,kind=8)
        EROOTS(10) = cmplx(10d0,0d0,kind=8)
        EROOTS(11) = cmplx(11d0,0d0,kind=8)
        EROOTS(12) = cmplx(12d0,0d0,kind=8)
        EROOTS(13) = cmplx(13d0,0d0,kind=8)
        EROOTS(14) = cmplx(14d0,0d0,kind=8)
        EROOTS(15) = cmplx(15d0,0d0,kind=8)
        COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  2) = cmplx(-1.200000000000000000D+02, 0.000000000000000000D+00,kind=8)
        COEFFS(  3) = cmplx( 6.580000000000000000D+03, 0.000000000000000000D+00,kind=8)
        COEFFS(  4) = cmplx(-2.184000000000000000D+05, 0.000000000000000000D+00,kind=8)
        COEFFS(  5) = cmplx( 4.899622000000000000D+06, 0.000000000000000000D+00,kind=8)
        COEFFS(  6) = cmplx(-7.855848000000000000D+07, 0.000000000000000000D+00,kind=8)
        COEFFS(  7) = cmplx( 9.280957400000000000D+08, 0.000000000000000000D+00,kind=8)
        COEFFS(  8) = cmplx(-8.207628000000000000D+09, 0.000000000000000000D+00,kind=8)
        COEFFS(  9) = cmplx( 5.463112955300000000D+10, 0.000000000000000000D+00,kind=8)
        COEFFS( 10) = cmplx(-2.728032106800000000D+11, 0.000000000000000000D+00,kind=8)
        COEFFS( 11) = cmplx( 1.009672107080000000D+12, 0.000000000000000000D+00,kind=8)
        COEFFS( 12) = cmplx(-2.706813345600000000D+12, 0.000000000000000000D+00,kind=8)
        COEFFS( 13) = cmplx( 5.056995703824000000D+12, 0.000000000000000000D+00,kind=8)
        COEFFS( 14) = cmplx(-6.165817614720000000D+12, 0.000000000000000000D+00,kind=8)
        COEFFS( 15) = cmplx( 4.339163001600000000D+12, 0.000000000000000000D+00,kind=8)
        COEFFS( 16) = cmplx(-1.307674368000000000D+12, 0.000000000000000000D+00,kind=8)
        
     case (3) ! Wilkinson Polynomial
        print*, "Wilkinson Polynomial"
        eigsknown = .TRUE.
        EROOTS(1)  = cmplx(1d0,0d0,kind=8)
        EROOTS(2)  = cmplx(2d0,0d0,kind=8)
        EROOTS(3)  = cmplx(3d0,0d0,kind=8)
        EROOTS(4)  = cmplx(4d0,0d0,kind=8)
        EROOTS(5)  = cmplx(5d0,0d0,kind=8)
        EROOTS(6)  = cmplx(6d0,0d0,kind=8)
        EROOTS(7)  = cmplx(7d0,0d0,kind=8)
        EROOTS(8)  = cmplx(8d0,0d0,kind=8)
        EROOTS(9)  = cmplx(9d0,0d0,kind=8)
        EROOTS(10) = cmplx(10d0,0d0,kind=8)
        EROOTS(11) = cmplx(11d0,0d0,kind=8)
        EROOTS(12) = cmplx(12d0,0d0,kind=8)
        EROOTS(13) = cmplx(13d0,0d0,kind=8)
        EROOTS(14) = cmplx(14d0,0d0,kind=8)
        EROOTS(15) = cmplx(15d0,0d0,kind=8)
        EROOTS(16) = cmplx(16d0,0d0,kind=8)
        EROOTS(17) = cmplx(17d0,0d0,kind=8)
        EROOTS(18) = cmplx(18d0,0d0,kind=8)
        EROOTS(19) = cmplx(19d0,0d0,kind=8)
        EROOTS(20) = cmplx(20d0,0d0,kind=8)
        COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  2) = cmplx(-2.100000000000000000D+02, 0.000000000000000000D+00,kind=8)
        COEFFS(  3) = cmplx( 2.061500000000000000D+04, 0.000000000000000000D+00,kind=8)
        COEFFS(  4) = cmplx(-1.256850000000000000D+06, 0.000000000000000000D+00,kind=8)
        COEFFS(  5) = cmplx( 5.332794600000000000D+07, 0.000000000000000000D+00,kind=8)
        COEFFS(  6) = cmplx(-1.672280820000000000D+09, 0.000000000000000000D+00,kind=8)
        COEFFS(  7) = cmplx( 4.017177163000000000D+10, 0.000000000000000000D+00,kind=8)
        COEFFS(  8) = cmplx(-7.561111845000000000D+11, 0.000000000000000000D+00,kind=8)
        COEFFS(  9) = cmplx( 1.131027699538100000D+13, 0.000000000000000000D+00,kind=8)
        COEFFS( 10) = cmplx(-1.355851828995300000D+14, 0.000000000000000000D+00,kind=8)
        COEFFS( 11) = cmplx( 1.307535010540395000D+15, 0.000000000000000000D+00,kind=8)
        COEFFS( 12) = cmplx(-1.014229986551145000D+16, 0.000000000000000000D+00,kind=8)
        COEFFS( 13) = cmplx( 6.303081209929489600D+16, 0.000000000000000000D+00,kind=8)
        COEFFS( 14) = cmplx(-3.113336431613906560D+17, 0.000000000000000000D+00,kind=8)
        COEFFS( 15) = cmplx( 1.206647803780373248D+18, 0.000000000000000000D+00,kind=8)
        COEFFS( 16) = cmplx(-3.599979517947607040D+18, 0.000000000000000000D+00,kind=8)
        COEFFS( 17) = cmplx( 8.037811822645051392D+18, 0.000000000000000000D+00,kind=8)
        COEFFS( 18) = cmplx(-1.287093124515098829D+19, 0.000000000000000000D+00,kind=8)
        COEFFS( 19) = cmplx( 1.380375975364070400D+19, 0.000000000000000000D+00,kind=8)
        COEFFS( 20) = cmplx(-8.752948036761600000D+18, 0.000000000000000000D+00,kind=8)
        COEFFS( 21) = cmplx( 2.432902008176640000D+18, 0.000000000000000000D+00,kind=8)

     case (4)  ! monic polynomial roots [-2.1:0.2:1.7]
        print*, "monic polynomial roots [-2.1:0.2:1.7]"
        eigsknown = .TRUE.           
        EROOTS(1)  = cmplx(-2.1d0,0d0,kind=8)
        EROOTS(2)  = cmplx(-1.9d0,0d0,kind=8)
        EROOTS(3)  = cmplx(-1.7d0,0d0,kind=8)
        EROOTS(4)  = cmplx(-1.5d0,0d0,kind=8)
        EROOTS(5)  = cmplx(-1.3d0,0d0,kind=8)
        EROOTS(6)  = cmplx(-1.1d0,0d0,kind=8)
        EROOTS(7)  = cmplx(-0.9d0,0d0,kind=8)
        EROOTS(8)  = cmplx(-0.7d0,0d0,kind=8)
        EROOTS(9)  = cmplx(-0.5d0,0d0,kind=8)
        EROOTS(10) = cmplx(-0.3d0,0d0,kind=8)
        EROOTS(11) = cmplx(-0.1d0,0d0,kind=8)
        EROOTS(12) = cmplx(0.1d0,0d0,kind=8)
        EROOTS(13) = cmplx(0.3d0,0d0,kind=8)
        EROOTS(14) = cmplx(0.5d0,0d0,kind=8)
        EROOTS(15) = cmplx(0.7d0,0d0,kind=8)
        EROOTS(16) = cmplx(0.9d0,0d0,kind=8)
        EROOTS(17) = cmplx(1.1d0,0d0,kind=8)
        EROOTS(18) = cmplx(1.3d0,0d0,kind=8)
        EROOTS(19) = cmplx(1.5d0,0d0,kind=8)
        EROOTS(20) = cmplx(1.7d0,0d0,kind=8)
        
        COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  2) = cmplx( 4.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  3) = cmplx(-5.700000000000000178D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  4) = cmplx(-3.875999999999999801D+01, 0.000000000000000000D+00,kind=8)
        COEFFS(  5) = cmplx(-1.065899999999999181D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  6) = cmplx( 1.503888000000000034D+02, 0.000000000000000000D+00,kind=8)
        COEFFS(  7) = cmplx( 7.471635999999999456D+01, 0.000000000000000000D+00,kind=8)
        COEFFS(  8) = cmplx(-3.011858720000000176D+02, 0.000000000000000000D+00,kind=8)
        COEFFS(  9) = cmplx(-2.166913265400000057D+02, 0.000000000000000000D+00,kind=8)
        COEFFS( 10) = cmplx( 3.349663231200000268D+02, 0.000000000000000000D+00,kind=8)
        COEFFS( 11) = cmplx( 2.822715856739999936D+02, 0.000000000000000000D+00,kind=8)
        COEFFS( 12) = cmplx(-2.074292865528000220D+02, 0.000000000000000000D+00,kind=8)
        COEFFS( 13) = cmplx(-1.899148773100940275D+02, 0.000000000000000000D+00,kind=8)
        COEFFS( 14) = cmplx( 6.798334410529599836D+01, 0.000000000000000000D+00,kind=8)
        COEFFS( 15) = cmplx( 6.520295820747719517D+01, 0.000000000000000000D+00,kind=8)
        COEFFS( 16) = cmplx(-1.044171015022224047D+01, 0.000000000000000000D+00,kind=8)
        COEFFS( 17) = cmplx(-1.027240495854272240D+01, 0.000000000000000000D+00,kind=8)
        COEFFS( 18) = cmplx( 5.728036652158500219D-01, 0.000000000000000000D+00,kind=8)
        COEFFS( 19) = cmplx( 5.701842040814797397D-01, 0.000000000000000000D+00,kind=8)
        COEFFS( 20) = cmplx(-4.749807885322501061D-03, 0.000000000000000000D+00,kind=8)
        COEFFS( 21) = cmplx(-4.737933365609194557D-03, 0.000000000000000000D+00,kind=8)

     case (5)  ! inverse Wilkinson polynomial degree 10
        print*, "inverse Wilkinson polynomial degree", 10
        do ii=1,n
           EROOTS(ii) = cmplx(1d0/ii,0d0,kind=8)
        end do
        eigsknown = .TRUE.

        COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  2) = cmplx(-2.928968253968253777D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  3) = cmplx( 3.514543650793650720D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  4) = cmplx(-2.317432760141093340D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  5) = cmplx( 9.416143077601411004D-01, 0.000000000000000000D+00,kind=8)
        COEFFS(  6) = cmplx(-2.485821759259259078D-01, 0.000000000000000000D+00,kind=8)
        COEFFS(  7) = cmplx( 4.347800925925925791D-02, 0.000000000000000000D+00,kind=8)
        COEFFS(  8) = cmplx(-5.001653439153438442D-03, 0.000000000000000000D+00,kind=8)
        COEFFS(  9) = cmplx( 3.637566137566137473D-04, 0.000000000000000000D+00,kind=8)
        COEFFS( 10) = cmplx(-1.515652557319223834D-05, 0.000000000000000000D+00,kind=8)
        COEFFS( 11) = cmplx( 2.755731922398588828D-07, 0.000000000000000000D+00,kind=8)

     case (6)  ! inverse Wilkinson polynomial degree 15
        print*, "inverse Wilkinson polynomial degree", 15
        do ii=1,n
           EROOTS(ii) = cmplx(1d0/ii,0d0,kind=8)
        end do
        eigsknown = .TRUE.
        
        COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  2) = cmplx(-3.318228993228993229D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  3) = cmplx( 4.715101684030255313D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  4) = cmplx(-3.867167413825151723D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  5) = cmplx( 2.069944484527817874D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  6) = cmplx(-7.721127918292269854D-01, 0.000000000000000000D+00,kind=8)
        COEFFS(  7) = cmplx( 2.086170818636096380D-01, 0.000000000000000000D+00,kind=8)
        COEFFS(  8) = cmplx(-4.177731925460512891D-02, 0.000000000000000000D+00,kind=8)
        COEFFS(  9) = cmplx( 6.276507516586881524D-03, 0.000000000000000000D+00,kind=8)
        COEFFS( 10) = cmplx(-7.097300082584473077D-04, 0.000000000000000000D+00,kind=8)
        COEFFS( 11) = cmplx( 6.007495590828923856D-05, 0.000000000000000000D+00,kind=8)
        COEFFS( 12) = cmplx(-3.746821165802646595D-06, 0.000000000000000000D+00,kind=8)
        COEFFS( 13) = cmplx( 1.670140559029447750D-07, 0.000000000000000000D+00,kind=8)
        COEFFS( 14) = cmplx(-5.031833735537438450D-09, 0.000000000000000000D+00,kind=8)
        COEFFS( 15) = cmplx( 9.176596478183778313D-11, 0.000000000000000000D+00,kind=8)
        COEFFS( 16) = cmplx(-7.647163731819815396D-13, 0.000000000000000000D+00,kind=8)
        
     case (7)  ! inverse Wilkinson polynomial degree 20
        print*, "inverse Wilkinson polynomial degree", 20
        do ii=1,n
           EROOTS(ii) = cmplx(1d0/ii,0d0,kind=8)
        end do
        eigsknown = .TRUE.
        
        COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  2) = cmplx(-3.597739657143681935D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  3) = cmplx( 5.673783698335657100D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  4) = cmplx(-5.290361552538329626D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  5) = cmplx( 3.303795958748482864D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  6) = cmplx(-1.479705925618288065D+00, 0.000000000000000000D+00,kind=8)
        COEFFS(  7) = cmplx( 4.959705733009387441D-01, 0.000000000000000000D+00,kind=8)
        COEFFS(  8) = cmplx(-1.279680160216244700D-01, 0.000000000000000000D+00,kind=8)
        COEFFS(  9) = cmplx( 2.590766577834094017D-02, 0.000000000000000000D+00,kind=8)
        COEFFS( 10) = cmplx(-4.168807387812830528D-03, 0.000000000000000000D+00,kind=8)
        COEFFS( 11) = cmplx( 5.374384196921842176D-04, 0.000000000000000000D+00,kind=8)
        COEFFS( 12) = cmplx(-5.572981667319412582D-05, 0.000000000000000000D+00,kind=8)
        COEFFS( 13) = cmplx( 4.648883085865667695D-06, 0.000000000000000000D+00,kind=8)
        COEFFS( 14) = cmplx(-3.107857126833785075D-07, 0.000000000000000000D+00,kind=8)
        COEFFS( 15) = cmplx( 1.651187408904606119D-08, 0.000000000000000000D+00,kind=8)
        COEFFS( 16) = cmplx(-6.873605325572917044D-10, 0.000000000000000000D+00,kind=8)
        COEFFS( 17) = cmplx( 2.191947962588394225D-11, 0.000000000000000000D+00,kind=8)
        COEFFS( 18) = cmplx(-5.166052704859893225D-13, 0.000000000000000000D+00,kind=8)
        COEFFS( 19) = cmplx( 8.473419780458025638D-15, 0.000000000000000000D+00,kind=8)
        COEFFS( 20) = cmplx(-8.631667008955544507D-17, 0.000000000000000000D+00,kind=8)
        COEFFS( 21) = cmplx( 4.110317623312163881D-19, 0.000000000000000000D+00,kind=8)

        case (8)  ! univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9
           print*, "univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9"         
           eigsknown = .TRUE.
           EROOTS(1)  = cmplx(1d0/1024d0,0d0,kind=8)
           EROOTS(2)  = cmplx(1d0/512d0,0d0,kind=8)
           EROOTS(3)  = cmplx(1d0/256d0,0d0,kind=8)
           EROOTS(4)  = cmplx(1d0/128d0,0d0,kind=8)
           EROOTS(5)  = cmplx(1d0/64d0,0d0,kind=8)
           EROOTS(6)  = cmplx(1d0/32d0,0d0,kind=8)
           EROOTS(7)  = cmplx(0.0625d0,0d0,kind=8)
           EROOTS(8)  = cmplx(0.125d0,0d0,kind=8)
           EROOTS(9)  = cmplx(0.25d0,0d0,kind=8)
           EROOTS(10) = cmplx(0.5d0,0d0,kind=8)
           EROOTS(11) = cmplx(1d0,0d0,kind=8)
           EROOTS(12) = cmplx(2d0,0d0,kind=8)
           EROOTS(13) = cmplx(4d0,0d0,kind=8)
           EROOTS(14) = cmplx(8d0,0d0,kind=8)
           EROOTS(15) = cmplx(16d0,0d0,kind=8)
           EROOTS(16) = cmplx(32d0,0d0,kind=8)
           EROOTS(17) = cmplx(64d0,0d0,kind=8)
           EROOTS(18) = cmplx(128d0,0d0,kind=8)
           EROOTS(19) = cmplx(256d0,0d0,kind=8)
           EROOTS(20) = cmplx(512d0,0d0,kind=8)

           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx(-1.023999023437500000D+03, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx( 3.495243333339691162D+05, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx(-5.113022171493675560D+07, 0.000000000000000000D+00,kind=8)
           COEFFS(  5) = cmplx( 3.490463172082539082D+09, 0.000000000000000000D+00,kind=8)
           COEFFS(  6) = cmplx(-1.152961209588856964D+11, 0.000000000000000000D+00,kind=8)
           COEFFS(  7) = cmplx( 1.873962299335221924D+12, 0.000000000000000000D+00,kind=8)
           COEFFS(  8) = cmplx(-1.510882103839022656D+13, 0.000000000000000000D+00,kind=8)
           COEFFS(  9) = cmplx( 6.066487898306585938D+13, 0.000000000000000000D+00,kind=8)
           COEFFS( 10) = cmplx(-1.215375144010052344D+14, 0.000000000000000000D+00,kind=8)
           COEFFS( 11) = cmplx( 1.215969169007124688D+14, 0.000000000000000000D+00,kind=8)
           COEFFS( 12) = cmplx(-6.076875720050261719D+13, 0.000000000000000000D+00,kind=8)
           COEFFS( 13) = cmplx( 1.516621974576646484D+13, 0.000000000000000000D+00,kind=8)
           COEFFS( 14) = cmplx(-1.888602629798778320D+12, 0.000000000000000000D+00,kind=8)
           COEFFS( 15) = cmplx( 1.171226437084513702D+11, 0.000000000000000000D+00,kind=8)
           COEFFS( 16) = cmplx(-3.603003779965178013D+09, 0.000000000000000000D+00,kind=8)
           COEFFS( 17) = cmplx( 5.453848706378967315D+07, 0.000000000000000000D+00,kind=8)
           COEFFS( 18) = cmplx(-3.994548571479434031D+05, 0.000000000000000000D+00,kind=8)
           COEFFS( 19) = cmplx( 1.365329427085816860D+03, 0.000000000000000000D+00,kind=8)
           COEFFS( 20) = cmplx(-1.999998092651367188D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 21) = cmplx( 9.765625000000000000D-04, 0.000000000000000000D+00,kind=8)
           
        case (9)  ! univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9
           print*, "univariate polynomial with zeros 2^-10, 2^-9, ..., 2^8, 2^9"
           eigsknown = .TRUE.
           EROOTS(1)  = cmplx(1d0/1024d0-3d0,0d0,kind=8)
           EROOTS(2)  = cmplx(1d0/512d0-3d0,0d0,kind=8)
           EROOTS(3)  = cmplx(1d0/256d0-3d0,0d0,kind=8)
           EROOTS(4)  = cmplx(1d0/128d0-3d0,0d0,kind=8)
           EROOTS(5)  = cmplx(1d0/64d0-3d0,0d0,kind=8)
           EROOTS(6)  = cmplx(1d0/32d0-3d0,0d0,kind=8)
           EROOTS(7)  = cmplx(0.0625d0-3d0,0d0,kind=8)
           EROOTS(8)  = cmplx(0.125d0-3d0,0d0,kind=8)
           EROOTS(9)  = cmplx(0.25d0-3d0,0d0,kind=8)
           EROOTS(10) = cmplx(0.5d0-3d0,0d0,kind=8)
           EROOTS(11) = cmplx(1d0-3d0,0d0,kind=8)
           EROOTS(12) = cmplx(2d0-3d0,0d0,kind=8)
           EROOTS(13) = cmplx(4d0-3d0,0d0,kind=8)
           EROOTS(14) = cmplx(8d0-3d0,0d0,kind=8)
           EROOTS(15) = cmplx(16d0-3d0,0d0,kind=8)
           EROOTS(16) = cmplx(32d0-3d0,0d0,kind=8)
           EROOTS(17) = cmplx(64d0-3d0,0d0,kind=8)
           EROOTS(18) = cmplx(128d0-3d0,0d0,kind=8)
           EROOTS(19) = cmplx(256d0-3d0,0d0,kind=8)
           EROOTS(20) = cmplx(512d0-3d0,0d0,kind=8)

           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx(-9.639990234375000000D+02, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx( 2.928663889980316162D+05, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx(-3.380106221197273582D+07, 0.000000000000000000D+00,kind=8)
           COEFFS(  5) = cmplx( 1.337718430171444893D+09, 0.000000000000000000D+00,kind=8)
           COEFFS(  6) = cmplx(-2.954283131394730568D+09, 0.000000000000000000D+00,kind=8)
           COEFFS(  7) = cmplx(-3.996458470881940918D+11, 0.000000000000000000D+00,kind=8)
           COEFFS(  8) = cmplx(-1.730686076106421143D+12, 0.000000000000000000D+00,kind=8)
           COEFFS(  9) = cmplx( 3.208167078493442578D+13, 0.000000000000000000D+00,kind=8)
           COEFFS( 10) = cmplx( 3.929729145248658125D+14, 0.000000000000000000D+00,kind=8)
           COEFFS( 11) = cmplx( 1.583818591945143250D+15, 0.000000000000000000D+00,kind=8)
           COEFFS( 12) = cmplx(-7.674659101259701250D+14, 0.000000000000000000D+00,kind=8)
           COEFFS( 13) = cmplx(-3.560718163405644800D+16, 0.000000000000000000D+00,kind=8)
           COEFFS( 14) = cmplx(-1.723861401101336320D+17, 0.000000000000000000D+00,kind=8)
           COEFFS( 15) = cmplx(-4.426679893231359360D+17, 0.000000000000000000D+00,kind=8)
           COEFFS( 16) = cmplx(-6.494223815728798720D+17, 0.000000000000000000D+00,kind=8)
           COEFFS( 17) = cmplx(-4.124865865242539520D+17, 0.000000000000000000D+00,kind=8)
           COEFFS( 18) = cmplx( 2.669806373742526720D+17, 0.000000000000000000D+00,kind=8)
           COEFFS( 19) = cmplx( 7.356501255213547520D+17, 0.000000000000000000D+00,kind=8)
           COEFFS( 20) = cmplx( 5.552041109785232640D+17, 0.000000000000000000D+00,kind=8)
           COEFFS( 21) = cmplx( 1.534961300051973760D+17, 0.000000000000000000D+00,kind=8)

        case (10) ! Chebyshev polynomial of degree 20 
           print*, "Chebyshev polynomial of degree 20"
           eigsknown = .TRUE.
           EROOTS(1)  = cmplx(cos(1d0/40d0*pi),0d0,kind=8)
           EROOTS(2)  = cmplx(cos(3d0/40d0*pi),0d0,kind=8)
           EROOTS(3)  = cmplx(cos(5d0/40d0*pi),0d0,kind=8)
           EROOTS(4)  = cmplx(cos(7d0/40d0*pi),0d0,kind=8)
           EROOTS(5)  = cmplx(cos(9d0/40d0*pi),0d0,kind=8)
           EROOTS(6)  = cmplx(cos(11d0/40d0*pi),0d0,kind=8)
           EROOTS(7)  = cmplx(cos(13d0/40d0*pi),0d0,kind=8)
           EROOTS(8)  = cmplx(cos(15d0/40d0*pi),0d0,kind=8)
           EROOTS(9)  = cmplx(cos(17d0/40d0*pi),0d0,kind=8)
           EROOTS(10)  = cmplx(cos(19d0/40d0*pi),0d0,kind=8)
           EROOTS(11)  = cmplx(cos(21d0/40d0*pi),0d0,kind=8)
           EROOTS(12)  = cmplx(cos(23d0/40d0*pi),0d0,kind=8)
           EROOTS(13)  = cmplx(cos(25d0/40d0*pi),0d0,kind=8)
           EROOTS(14)  = cmplx(cos(27d0/40d0*pi),0d0,kind=8)
           EROOTS(15)  = cmplx(cos(29d0/40d0*pi),0d0,kind=8)
           EROOTS(16)  = cmplx(cos(31d0/40d0*pi),0d0,kind=8)
           EROOTS(17)  = cmplx(cos(33d0/40d0*pi),0d0,kind=8)
           EROOTS(18)  = cmplx(cos(35d0/40d0*pi),0d0,kind=8)
           EROOTS(19)  = cmplx(cos(37d0/40d0*pi),0d0,kind=8)
           EROOTS(20)  = cmplx(cos(39d0/40d0*pi),0d0,kind=8)         
           
           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx( 0.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx(-5.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx( 0.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  5) = cmplx( 1.062500000000000000D+01, 0.000000000000000000D+00,kind=8)
           COEFFS(  6) = cmplx( 0.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  7) = cmplx(-1.250000000000000000D+01, 0.000000000000000000D+00,kind=8)
           COEFFS(  8) = cmplx( 0.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  9) = cmplx( 8.886718750000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 10) = cmplx( 0.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 11) = cmplx(-3.910156250000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 12) = cmplx( 0.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 13) = cmplx( 1.047363281250000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 14) = cmplx( 0.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 15) = cmplx(-1.611328125000000000D-01, 0.000000000000000000D+00,kind=8)
           COEFFS( 16) = cmplx( 0.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 17) = cmplx( 1.258850097656250000D-02, 0.000000000000000000D+00,kind=8)
           COEFFS( 18) = cmplx( 0.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 19) = cmplx(-3.814697265625000000D-04, 0.000000000000000000D+00,kind=8)
           COEFFS( 20) = cmplx( 0.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 21) = cmplx( 1.907348632812500000D-06, 0.000000000000000000D+00,kind=8)


        case (11)   ! z^20 + z^19 + ... + z + 1
           print*, "z^20 + z^19 + ... + z + 1"
           eigsknown = .TRUE.
           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  5) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  6) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  7) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  8) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  9) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 10) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 11) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 12) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 13) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 14) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 15) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 16) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 17) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 18) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 19) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 20) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 21) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)

           ! roots of unity 21 without 1
           EROOTS(1)  = cmplx(cos(2d0/21d0*pi),sin(2d0/21d0*pi),kind=8)
           EROOTS(2)  = cmplx(cos(4d0/21d0*pi),sin(4d0/21d0*pi),kind=8)
           EROOTS(3)  = cmplx(cos(6d0/21d0*pi),sin(6d0/21d0*pi),kind=8)
           EROOTS(4)  = cmplx(cos(8d0/21d0*pi),sin(8d0/21d0*pi),kind=8)
           EROOTS(5)  = cmplx(cos(10d0/21d0*pi),sin(10d0/21d0*pi),kind=8)
           EROOTS(6)  = cmplx(cos(12d0/21d0*pi),sin(12d0/21d0*pi),kind=8)
           EROOTS(7)  = cmplx(cos(14d0/21d0*pi),sin(14d0/21d0*pi),kind=8)
           EROOTS(8)  = cmplx(cos(16d0/21d0*pi),sin(16d0/21d0*pi),kind=8)
           EROOTS(9)  = cmplx(cos(18d0/21d0*pi),sin(18d0/21d0*pi),kind=8)
           EROOTS(10)  = cmplx(cos(20d0/21d0*pi),sin(20d0/21d0*pi),kind=8)
           EROOTS(11)  = cmplx(cos(22d0/21d0*pi),sin(22d0/21d0*pi),kind=8)
           EROOTS(12)  = cmplx(cos(24d0/21d0*pi),sin(24d0/21d0*pi),kind=8)
           EROOTS(13)  = cmplx(cos(26d0/21d0*pi),sin(26d0/21d0*pi),kind=8)
           EROOTS(14)  = cmplx(cos(28d0/21d0*pi),sin(28d0/21d0*pi),kind=8)
           EROOTS(15)  = cmplx(cos(30d0/21d0*pi),sin(30d0/21d0*pi),kind=8)
           EROOTS(16)  = cmplx(cos(32d0/21d0*pi),sin(32d0/21d0*pi),kind=8)
           EROOTS(17)  = cmplx(cos(34d0/21d0*pi),sin(34d0/21d0*pi),kind=8)
           EROOTS(18)  = cmplx(cos(36d0/21d0*pi),sin(36d0/21d0*pi),kind=8)
           EROOTS(19)  = cmplx(cos(38d0/21d0*pi),sin(38d0/21d0*pi),kind=8)
           EROOTS(20)  = cmplx(cos(40d0/21d0*pi),sin(40d0/21d0*pi),kind=8)

        case (12)
           print*, "C. Traverso 24 MPSolve"
           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx(-1.728000000000000000D+03, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx( 3.284992000000000000D+06, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx(-3.949133824000000000D+09, 0.000000000000000000D+00,kind=8)
           COEFFS(  5) = cmplx( 1.015406084096000000D+12, 0.000000000000000000D+00,kind=8)
           COEFFS(  6) = cmplx(-9.688873555722240000D+14, 0.000000000000000000D+00,kind=8)
           COEFFS(  7) = cmplx( 1.265052493274939392D+18, 0.000000000000000000D+00,kind=8)
           COEFFS(  8) = cmplx( 4.537862510900726989D+20, 0.000000000000000000D+00,kind=8)
           COEFFS(  9) = cmplx( 2.021488014436448023D+22, 0.000000000000000000D+00,kind=8)
           COEFFS( 10) = cmplx(-1.431502639275795929D+26, 0.000000000000000000D+00,kind=8)
           COEFFS( 11) = cmplx(-2.047976871739763688D+29, 0.000000000000000000D+00,kind=8)
           COEFFS( 12) = cmplx(-7.561785527781800281D+31, 0.000000000000000000D+00,kind=8)
           COEFFS( 13) = cmplx( 2.206941937668751764D+33, 0.000000000000000000D+00,kind=8)
           COEFFS( 14) = cmplx( 2.337502650696833264D+37, 0.000000000000000000D+00,kind=8)
           COEFFS( 15) = cmplx( 1.657650689150882436D+40, 0.000000000000000000D+00,kind=8)
           COEFFS( 16) = cmplx( 5.968852409133617603D+42, 0.000000000000000000D+00,kind=8)
           COEFFS( 17) = cmplx( 1.971787195201967199D+44, 0.000000000000000000D+00,kind=8)
           COEFFS( 18) = cmplx(-7.269074194037149894D+47, 0.000000000000000000D+00,kind=8)
           COEFFS( 19) = cmplx(-1.955142887477579909D+50, 0.000000000000000000D+00,kind=8)
           COEFFS( 20) = cmplx( 6.611223380089858717D+51, 0.000000000000000000D+00,kind=8)
           COEFFS( 21) = cmplx( 7.337981286595499208D+54, 0.000000000000000000D+00,kind=8)
           COEFFS( 22) = cmplx( 5.750602254715702767D+56, 0.000000000000000000D+00,kind=8)
           COEFFS( 23) = cmplx(-3.196998408115594373D+58, 0.000000000000000000D+00,kind=8)
           COEFFS( 24) = cmplx(-4.052135566767965996D+60, 0.000000000000000000D+00,kind=8)
           COEFFS( 25) = cmplx(-5.476529142819802481D+61, 0.000000000000000000D+00,kind=8)
           eigsknown = .TRUE.         
           EROOTS( 1 )=cmplx(-3.52d2, 0d0,kind=8)
           EROOTS( 2 )=cmplx(-3.52d2, 0d0,kind=8)
           EROOTS( 3 )=cmplx(-2.8371450777d2, -2.9920517772d2,kind=8)
           EROOTS( 4 )=cmplx(-2.8371450777d2,  2.9920517772d2,kind=8)
           EROOTS( 5 )=cmplx(-2.7867414048d2,  6.1005469197d2,kind=8)
           EROOTS( 6 )=cmplx(-2.7867414048d2, -6.1005469197d2,kind=8)
           EROOTS( 7 )=cmplx(-2.74892372d2, 0d0,kind=8)
           EROOTS( 8 )=cmplx(-2.014171531d2, 0d0,kind=8)
           EROOTS( 9 )=cmplx(-1.255366582d2, 0d0,kind=8)
           EROOTS( 10 )=cmplx(-9.599999999d1, 0d0,kind=8)
           EROOTS( 11 )=cmplx(-8.8692435121d1,  5.5009607430d2,kind=8)
           EROOTS( 12 )=cmplx(-8.869243512d1, -5.5009607430d2,kind=8)
           EROOTS( 13 )=cmplx(-1.6000000000d1, 0d0,kind=8)
           EROOTS( 14 )=cmplx( 8.23178509855d1, 0d0,kind=8)
           EROOTS( 15 )=cmplx( 8.8692435121d1, -5.50096074303d2,kind=8)
           EROOTS( 16 )=cmplx( 8.8692435121d1,  5.5009607430d2,kind=8)
           EROOTS( 17 )=cmplx( 1.9293739373d2,  1.60865921259d3,kind=8)
           EROOTS( 18 )=cmplx( 1.929373937d2, -1.6086592125d3,kind=8)
           EROOTS( 19 )=cmplx( 2.0141715312d2, 0d0,kind=8)
           EROOTS( 20 )=cmplx( 2.7489237213d2, 0d0,kind=8)
           EROOTS( 21 )=cmplx( 7.52d2, 0d0,kind=8)
           EROOTS( 22 )=cmplx( 7.52d2, 0d0,kind=8)
           EROOTS( 23 )=cmplx( 9.1106065d2,  1.5722d0,kind=8)
           EROOTS( 24 )=cmplx( 9.1106065d2, -1.5722d0,kind=8)

        case (13) ! Mandelbort 31 MPSolve
           print*, "Mandelbrot 31 MPSolve"
           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx( 1.600000000000000000D+01, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx( 1.200000000000000000D+02, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx( 5.680000000000000000D+02, 0.000000000000000000D+00,kind=8)
           COEFFS(  5) = cmplx( 1.932000000000000000D+03, 0.000000000000000000D+00,kind=8)
           COEFFS(  6) = cmplx( 5.096000000000000000D+03, 0.000000000000000000D+00,kind=8)
           COEFFS(  7) = cmplx( 1.094800000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS(  8) = cmplx( 1.978800000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS(  9) = cmplx( 3.078200000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 10) = cmplx( 4.194400000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 11) = cmplx( 5.078800000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 12) = cmplx( 5.530800000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 13) = cmplx( 5.474600000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 14) = cmplx( 4.970000000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 15) = cmplx( 4.165800000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 16) = cmplx( 3.239800000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 17) = cmplx( 2.346100000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 18) = cmplx( 1.586400000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 19) = cmplx( 1.006800000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 20) = cmplx( 6.036000000000000000D+03, 0.000000000000000000D+00,kind=8)
           COEFFS( 21) = cmplx( 3.434000000000000000D+03, 0.000000000000000000D+00,kind=8)
           COEFFS( 22) = cmplx( 1.860000000000000000D+03, 0.000000000000000000D+00,kind=8)
           COEFFS( 23) = cmplx( 9.580000000000000000D+02, 0.000000000000000000D+00,kind=8)
           COEFFS( 24) = cmplx( 4.700000000000000000D+02, 0.000000000000000000D+00,kind=8)
           COEFFS( 25) = cmplx( 2.210000000000000000D+02, 0.000000000000000000D+00,kind=8)
           COEFFS( 26) = cmplx( 1.000000000000000000D+02, 0.000000000000000000D+00,kind=8)
           COEFFS( 27) = cmplx( 4.200000000000000000D+01, 0.000000000000000000D+00,kind=8)
           COEFFS( 28) = cmplx( 1.400000000000000000D+01, 0.000000000000000000D+00,kind=8)
           COEFFS( 29) = cmplx( 5.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 30) = cmplx( 2.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 31) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 32) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)

           eigsknown = .TRUE.          
           EROOTS( 1 ) = cmplx(-1.996376137,0d0,kind=8)
           EROOTS( 2 ) = cmplx(-1.966773216,0d0,kind=8)
           EROOTS( 3 ) = cmplx(-1.907280091,0d0,kind=8)
           EROOTS( 4 ) = cmplx(-1.772892903,0d0,kind=8)
           EROOTS( 5 ) = cmplx(-1.754877666,0d0,kind=8)
           EROOTS( 6 ) = cmplx(-1.47601464272,0d0,kind=8)
           EROOTS( 7 ) = cmplx(-1.284084925525, 4.272688960406d-1,kind=8)
           EROOTS( 8 ) = cmplx(-1.284084925525,-4.272688960406d-1,kind=8)
           EROOTS( 9 ) = cmplx(-1.138000666650,-2.403324012620d-1,kind=8)
           EROOTS( 10 )= cmplx(-1.138000666650, 2.403324012620d-1,kind=8)
           EROOTS( 11 )= cmplx(-1d0,0d0,kind=8)
           EROOTS( 12 )= cmplx(-5.968916446451269d-1, 6.629807445770295d-1,kind=8)
           EROOTS( 13 )= cmplx(-5.968916446451269d-1,-6.629807445770295d-1,kind=8)
           EROOTS( 14 )= cmplx(-2.17526747030511d-1,-1.11445426587329,kind=8)
           EROOTS( 15 )= cmplx(-2.17526747030511d-1, 1.11445426587329,kind=8)
           EROOTS( 16 )= cmplx(-1.6359826155202d-1, 1.09778064288827,kind=8)
           EROOTS( 17 )= cmplx(-1.6359826155202d-1,-1.09778064288827,kind=8)
           EROOTS( 18 )= cmplx(-1.225611668766536d-1,-7.4486176661974423d-1,kind=8)
           EROOTS( 19 )= cmplx(-1.225611668766536d-1, 7.4486176661974423d-1,kind=8)
           EROOTS( 20 )= cmplx(-1.13418655949436d-1,-8.605694725015730d-1,kind=8)
           EROOTS( 21 )= cmplx(-1.13418655949436d-1,8.605694725015730d-1,kind=8)
           EROOTS( 22 )= cmplx(-1.5570386020902d-2, 1.020497366498289d0,kind=8)
           EROOTS( 23 )= cmplx(-1.5570386020902d-2,-1.020497366498289d0,kind=8)
           EROOTS( 24 )= cmplx(3.59892739012579001d-1, 6.84762020211812856d-1,kind=8)
           EROOTS( 25 )= cmplx(3.59892739012579001d-1,-6.84762020211812856d-1,kind=8)
           EROOTS( 26 )= cmplx(3.8900684056977123543d-1,-2.1585065087081910777d-1,kind=8)
           EROOTS( 27 )= cmplx(3.8900684056977123543d-1, 2.1585065087081910777d-1,kind=8)
           EROOTS( 28 )= cmplx(3.96534570032415023d-1, 6.04181810488988837d-1,kind=8)
           EROOTS( 29 )= cmplx(3.96534570032415023d-1,-6.04181810488988837d-1,kind=8)
           EROOTS( 30 )= cmplx(4.433256333996235387d-1, 3.729624166628465083d-1,kind=8)
           EROOTS( 31 )= cmplx(4.433256333996235387d-1,-3.729624166628465083d-1,kind=8)

      case (14) ! Mandelbort 63 MPSolve
           print*, "Mandelbrot 63 MPSolve"
           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx( 3.200000000000000000D+01, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx( 4.960000000000000000D+02, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx( 4.976000000000000000D+03, 0.000000000000000000D+00,kind=8)
           COEFFS(  5) = cmplx( 3.644000000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS(  6) = cmplx( 2.083360000000000000D+05, 0.000000000000000000D+00,kind=8)
           COEFFS(  7) = cmplx( 9.712720000000000000D+05, 0.000000000000000000D+00,kind=8)
           COEFFS(  8) = cmplx( 3.807704000000000000D+06, 0.000000000000000000D+00,kind=8)
           COEFFS(  9) = cmplx( 1.284398000000000000D+07, 0.000000000000000000D+00,kind=8)
           COEFFS( 10) = cmplx( 3.794590400000000000D+07, 0.000000000000000000D+00,kind=8)
           COEFFS( 11) = cmplx( 9.958292000000000000D+07, 0.000000000000000000D+00,kind=8)
           COEFFS( 12) = cmplx( 2.348135920000000000D+08, 0.000000000000000000D+00,kind=8)
           COEFFS( 13) = cmplx( 5.021965000000000000D+08, 0.000000000000000000D+00,kind=8)
           COEFFS( 14) = cmplx( 9.819001680000000000D+08, 0.000000000000000000D+00,kind=8)
           COEFFS( 15) = cmplx( 1.766948340000000000D+09, 0.000000000000000000D+00,kind=8)
           COEFFS( 16) = cmplx( 2.943492972000000000D+09, 0.000000000000000000D+00,kind=8)
           COEFFS( 17) = cmplx( 4.562339774000000000D+09, 0.000000000000000000D+00,kind=8)
           COEFFS( 18) = cmplx( 6.609143792000000000D+09, 0.000000000000000000D+00,kind=8)
           COEFFS( 19) = cmplx( 8.984070856000000000D+09, 0.000000000000000000D+00,kind=8)
           COEFFS( 20) = cmplx( 1.150090186400000000D+10, 0.000000000000000000D+00,kind=8)
           COEFFS( 21) = cmplx( 1.391004352400000000D+10, 0.000000000000000000D+00,kind=8)
           COEFFS( 22) = cmplx( 1.594168477600000000D+10, 0.000000000000000000D+00,kind=8)
           COEFFS( 23) = cmplx( 1.735793770800000000D+10, 0.000000000000000000D+00,kind=8)
           COEFFS( 24) = cmplx( 1.799943337200000000D+10, 0.000000000000000000D+00,kind=8)
           COEFFS( 25) = cmplx( 1.781377799400000000D+10, 0.000000000000000000D+00,kind=8)
           COEFFS( 26) = cmplx( 1.685941079200000000D+10, 0.000000000000000000D+00,kind=8)
           COEFFS( 27) = cmplx( 1.528606570000000000D+10, 0.000000000000000000D+00,kind=8)
           COEFFS( 28) = cmplx( 1.329936233200000000D+10, 0.000000000000000000D+00,kind=8)
           COEFFS( 29) = cmplx( 1.112013616200000000D+10, 0.000000000000000000D+00,kind=8)
           COEFFS( 30) = cmplx( 8.948546308000000000D+09, 0.000000000000000000D+00,kind=8)
           COEFFS( 31) = cmplx( 6.939692682000000000D+09, 0.000000000000000000D+00,kind=8)
           COEFFS( 32) = cmplx( 5.193067630000000000D+09, 0.000000000000000000D+00,kind=8)
           COEFFS( 33) = cmplx( 3.754272037000000000D+09, 0.000000000000000000D+00,kind=8)
           COEFFS( 34) = cmplx( 2.625062128000000000D+09, 0.000000000000000000D+00,kind=8)
           COEFFS( 35) = cmplx( 1.777171560000000000D+09, 0.000000000000000000D+00,kind=8)
           COEFFS( 36) = cmplx( 1.166067016000000000D+09, 0.000000000000000000D+00,kind=8)
           COEFFS( 37) = cmplx( 7.421792840000000000D+08, 0.000000000000000000D+00,kind=8)
           COEFFS( 38) = cmplx( 4.585914320000000000D+08, 0.000000000000000000D+00,kind=8)
           COEFFS( 39) = cmplx( 2.752767160000000000D+08, 0.000000000000000000D+00,kind=8)
           COEFFS( 40) = cmplx( 1.606178600000000000D+08, 0.000000000000000000D+00,kind=8)
           COEFFS( 41) = cmplx( 9.114311400000000000D+07, 0.000000000000000000D+00,kind=8)
           COEFFS( 42) = cmplx( 5.032349600000000000D+07, 0.000000000000000000D+00,kind=8)
           COEFFS( 43) = cmplx( 2.704919600000000000D+07, 0.000000000000000000D+00,kind=8)
           COEFFS( 44) = cmplx( 1.416222000000000000D+07, 0.000000000000000000D+00,kind=8)
           COEFFS( 45) = cmplx( 7.228014000000000000D+06, 0.000000000000000000D+00,kind=8)
           COEFFS( 46) = cmplx( 3.598964000000000000D+06, 0.000000000000000000D+00,kind=8)
           COEFFS( 47) = cmplx( 1.749654000000000000D+06, 0.000000000000000000D+00,kind=8)
           COEFFS( 48) = cmplx( 8.310140000000000000D+05, 0.000000000000000000D+00,kind=8)
           COEFFS( 49) = cmplx( 3.857410000000000000D+05, 0.000000000000000000D+00,kind=8)
           COEFFS( 50) = cmplx( 1.750480000000000000D+05, 0.000000000000000000D+00,kind=8)
           COEFFS( 51) = cmplx( 7.768400000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 52) = cmplx( 3.370800000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 53) = cmplx( 1.429000000000000000D+04, 0.000000000000000000D+00,kind=8)
           COEFFS( 54) = cmplx( 5.916000000000000000D+03, 0.000000000000000000D+00,kind=8)
           COEFFS( 55) = cmplx( 2.398000000000000000D+03, 0.000000000000000000D+00,kind=8)
           COEFFS( 56) = cmplx( 9.500000000000000000D+02, 0.000000000000000000D+00,kind=8)
           COEFFS( 57) = cmplx( 3.650000000000000000D+02, 0.000000000000000000D+00,kind=8)
           COEFFS( 58) = cmplx( 1.320000000000000000D+02, 0.000000000000000000D+00,kind=8)
           COEFFS( 59) = cmplx( 4.200000000000000000D+01, 0.000000000000000000D+00,kind=8)
           COEFFS( 60) = cmplx( 1.400000000000000000D+01, 0.000000000000000000D+00,kind=8)
           COEFFS( 61) = cmplx( 5.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 62) = cmplx( 2.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 63) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS( 64) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           
           eigsknown = .TRUE.
           EROOTS( 1 )=cmplx(-1.999095682327018473210d0,0d0,kind=8)
           EROOTS( 2 )=cmplx(-1.9918141725491222157325609498622881d0,0d0,kind=8)
           EROOTS( 3 )=cmplx(-1.977179587006257387346088520662828616836d0,0d0,kind=8)
           EROOTS( 4 )=cmplx(-1.953705894284396245427622199013653238901d0,0d0,kind=8)
           EROOTS( 5 )=cmplx(-1.927147709363950262460068188946594278007d0,0d0,kind=8)
           EROOTS( 6 )=cmplx(-1.8848035715866817923294780929158396496359d0,0d0,kind=8)
           EROOTS( 7 )=cmplx(-1.8323152027512291920848975260425181432293d0,0d0,kind=8)
           EROOTS( 8 )=cmplx(-1.76926167027683114607548022863625740038777d0, 5.6919500395600315304900187298015859319654d-2,kind=8)
           EROOTS( 9 )=cmplx(-1.76926167027683114607548022863625740038777d0,-5.6919500395600315304900187298015859319654d-2,kind=8)
           EROOTS( 10 )=cmplx(-1.674066091474787971565296029172325596206403d0,0d0,kind=8)
           EROOTS( 11 )=cmplx(-1.5748891397523009698199655524959742837719482d0,0d0,kind=8)
           EROOTS( 12 )=cmplx(-1.408446485740072654917577008805998851928020904d0, &
                -1.36171997304659915684707793608163610038822995d-1,kind=8)
           EROOTS( 13 )=cmplx(-1.408446485740072654917577008805998851928020904d0, &
                1.36171997304659915684707793608163610038822995d-1,kind=8)
           EROOTS( 14 )=cmplx(-1.29255806103352208716418470636149411998013630326d0, &
                4.3819881608663183712973712432734844004535476504d-1,kind=8)
           EROOTS( 15 )=cmplx(-1.29255806103352208716418470636149411998013630326d0, &
                -4.3819881608663183712973712432734844004535476504d-1,kind=8)
           EROOTS( 16 )=cmplx(-1.26228728143847254301011194120806575232050489502d0, &
                4.0810432411269038329016065742601506306041169168d-1,kind=8)
           EROOTS( 17 )=cmplx(-1.26228728143847254301011194120806575232050489502d0, &
                -4.0810432411269038329016065742601506306041169168d-1,kind=8)
           EROOTS( 18 )=cmplx(-1.25273588401203794629581100256433997387062287256d0, &
                -3.4247064788975089386187578687092843396383393805d-1,kind=8)
           EROOTS( 19 )=cmplx(-1.25273588401203794629581100256433997387062287256d0, &
                3.4247064788975089386187578687092843396383393805d-1,kind=8)
           EROOTS( 20 )=cmplx(-1.02819385245481759930249745596731843328070508279421d0, &
                -3.61376517118561592479460832997830315786692639704085d-1,kind=8)
           EROOTS( 21 )=cmplx(-1.02819385245481759930249745596731843328070508279421d0, &
                3.61376517118561592479460832997830315786692639704085d-1,kind=8)
           EROOTS( 22 )=cmplx(-6.23532485956252757990016587001026776428072703359878868d-1, &
                6.81064414225239608090835812686561539088332735217609127d-1,kind=8)
           EROOTS( 23 )=cmplx(-6.23532485956252757990016587001026776428072703359878868d-1, &
                -6.81064414225239608090835812686561539088332735217609127d-1,kind=8)
           EROOTS( 24 )=cmplx(-6.2243629504129358796016350694723840189750985673649588591d-1, &
                4.2487843647562918431157443880525338683545992964599689876d-1,kind=8)
           EROOTS( 25 )=cmplx(-6.2243629504129358796016350694723840189750985673649588591d-1, &
                -4.2487843647562918431157443880525338683545992964599689876d-1,kind=8)
           EROOTS( 26 )=cmplx(-5.308278048599427289214772971196026578135170949646890946d-1, &
                6.682887255592057714440924655647011851367651843270734380d-1,kind=8)
           EROOTS( 27 )=cmplx(-5.308278048599427289214772971196026578135170949646890946d-1, &
                -6.682887255592057714440924655647011851367651843270734380d-1,kind=8)
           EROOTS( 28 )=cmplx(-2.72102461488938894219383324518026874585894699621947085d-1, &
                -8.42364690294128145503155708242929569550778268698265965d-1,kind=8)
           EROOTS( 29 )=cmplx(-2.72102461488938894219383324518026874585894699621947085d-1, &
                8.42364690294128145503155708242929569550778268698265965d-1,kind=8)
           EROOTS( 30 )=cmplx(-2.24915951286740054685326255204118310792682454680693d-1, &
                1.11626015745499183500126825424467009109873946082435d0,kind=8)
           EROOTS( 31 )=cmplx(-2.24915951286740054685326255204118310792682454680693d-1, &
                -1.11626015745499183500126825424467009109873946082435d0,kind=8)
           EROOTS( 32 )=cmplx(-2.0728383545566641282413385018667121332401155604017d-1, &
                1.11748077249496291137377567312207879579746389236127d0,kind=8)
           EROOTS( 33 )=cmplx(-2.0728383545566641282413385018667121332401155604017d-1, &
                -1.11748077249496291137377567312207879579746389236127d0,kind=8)
           EROOTS( 34 )=cmplx(-1.7457822113571696945156643266162905020167505710204d-1, &
                1.07142767145403118922964631021955987671322451961088d0,kind=8)
           EROOTS( 35 )=cmplx(-1.7457822113571696945156643266162905020167505710204d-1, &
                -1.07142767145403118922964631021955987671322451961088d0,kind=8)
           EROOTS( 36 )=cmplx(-1.57516053475965356164335109644674141293297577896685d-1, &
                -1.10900651411360717797175198615475582901468585712356d0,kind=8) 
           EROOTS( 37 )=cmplx(-1.57516053475965356164335109644674141293297577896685d-1, &
                1.10900651411360717797175198615475582901468585712356d0,kind=8)
           EROOTS( 38 )=cmplx(-1.274999735463630001995395653459879637298616757217284d-1, &
                9.874609094894567922074076807929788675642068522522938d-1,kind=8)
           EROOTS( 39 )=cmplx(-1.274999735463630001995395653459879637298616757217284d-1, &
                -9.874609094894567922074076807929788675642068522522938d-1,kind=8)
           EROOTS( 40 )=cmplx(-1.42334819203540667677618453136202688358025283954839d-2, &
                -1.0329147752136441093950134026551104360994260360822540d0,kind=8)
           EROOTS( 41 )=cmplx(-1.42334819203540667677618453136202688358025283954839d-2, &
                1.0329147752136441093950134026551104360994260360822540d0,kind=8)
           EROOTS( 42 )=cmplx(-6.98356849626139181796649107548406610452886379651341d-3, &
                -1.0036038622882895485307049669513531297649273745391915d0,kind=8)
           EROOTS( 43 )=cmplx(-6.98356849626139181796649107548406610452886379651341d-3, &
                1.0036038622882895485307049669513531297649273745391915d0,kind=8)
           EROOTS( 44 )=cmplx( 1.4895466603687646529815779208794106185666477731693128d-2, &
                -8.481487619084165277193311117832376290806619901265058603d-1,kind=8)
           EROOTS( 45 )=cmplx( 1.4895466603687646529815779208794106185666477731693128d-2, &
                8.481487619084165277193311117832376290806619901265058603d-1,kind=8)
           EROOTS( 46 )=cmplx( 1.211927861059064863147044434105037593859287800520963579338d-1, &
                6.1061169221075421167538724415035774824319702690063863369691d-1,kind=8)
           EROOTS( 47 )=cmplx( 1.211927861059064863147044434105037593859287800520963579338d-1, &
                -6.1061169221075421167538724415035774824319702690063863369691d-1,kind=8)
           EROOTS( 48 )=cmplx( 3.52482539722363278193253964052161589243593334212239870706d-1, &
                -6.98337239583330331258141954760484537633150485928512286760d-1,kind=8)
           EROOTS( 49 )=cmplx( 3.52482539722363278193253964052161589243593334212239870706d-1, &
                6.98337239583330331258141954760484537633150485928512286760d-1,kind=8)
           EROOTS( 50 )=cmplx( 3.7600868184676755970480431772902286888800357334637481632029182d-1, &
                -1.4474937132163286474711018201298830556966056842762643026975894d-1,kind=8)
           EROOTS( 51 )=cmplx( 3.7600868184676755970480431772902286888800357334637481632029182d-1, &
                1.4474937132163286474711018201298830556966056842762643026975894d-1,kind=8)
           EROOTS( 52 )=cmplx( 3.76893240379311323690004017968512473363482317941533875341d-1, &
                6.78568693190448141957540792996773280196881194582788907016d-1,kind=8)
           EROOTS( 53 )=cmplx( 3.76893240379311323690004017968512473363482317941533875341d-1, &
                -6.78568693190448141957540792996773280196881194582788907016d-1,kind=8)
           EROOTS( 54 )=cmplx( 3.865391765961580265082930869043677799284877313516569138807d-1, &
                5.693247113031029032137923571351905081619323911951388853856d-1,kind=8)
           EROOTS( 55 )=cmplx( 3.865391765961580265082930869043677799284877313516569138807d-1, &
                -5.693247113031029032137923571351905081619323911951388853856d-1,kind=8)
           EROOTS( 56 )=cmplx( 4.12916024722700479197334566382612257174765142865547121703d-1, &
                6.148067601433856949545497204007997358291659758563137777616d-1,kind=8)
           EROOTS( 57 )=cmplx( 4.12916024722700479197334566382612257174765142865547121703d-1, &
                -6.148067601433856949545497204007997358291659758563137777616d-1,kind=8)
           EROOTS( 58 )=cmplx( 4.3237619264199450782466964808692137388785063987699403620424125d-1, &
                2.267599044353486186978765599716989721202321914603899690444951d-1,kind=8)
           EROOTS( 59 )=cmplx( 4.3237619264199450782466964808692137388785063987699403620424125d-1, &
                -2.267599044353486186978765599716989721202321914603899690444951d-1,kind=8)
           EROOTS( 60 )=cmplx( 4.52774498724915493508803077732546131473562000961307327749350d-1, &
                -3.96170128033165002412596877271155937712569079351815707744770d-1,kind=8)
           EROOTS( 61 )=cmplx( 4.52774498724915493508803077732546131473562000961307327749350d-1, &
                3.96170128033165002412596877271155937712569079351815707744770d-1,kind=8)
           EROOTS( 62 )=cmplx( 4.56823285823316651283953236253270107801699459631329688710054d-1, &
                3.47758700883481983632188723200264206004781117755664551397643d-1,kind=8)
           EROOTS( 63 )=cmplx( 4.56823285823316651283953236253270107801699459631329688710054d-1, &
                -3.47758700883481983632188723200264206004781117755664551397643d-1,kind=8)

        case (15)  ! Vanni Noferini's example
           print*, "Vanni Noferini's example degree 12"
           eigsknown = .TRUE.
           EROOTS(  1) = cmplx(-7.049873051880334616D+08, 0.000000000000000000D+00,kind=8)
           EROOTS(  2) = cmplx(-2.180034227531623779D+12, 0.000000000000000000D+00,kind=8)
           EROOTS(  3) = cmplx( 1.682782742108504814D-01, 0.000000000000000000D+00,kind=8)
           EROOTS(  4) = cmplx( 1.596055627509852748D+00, 0.000000000000000000D+00,kind=8)
           EROOTS(  5) = cmplx(-6.644021946056043459D-01, 0.000000000000000000D+00,kind=8)
           EROOTS(  6) = cmplx(-1.319743558965743135D+00, 0.000000000000000000D+00,kind=8)
           EROOTS(  7) = cmplx(-1.778074633616713474D+00, 0.000000000000000000D+00,kind=8)
           EROOTS(  8) = cmplx(-8.664335223749246584D-01, 0.000000000000000000D+00,kind=8)
           EROOTS(  9) = cmplx( 2.235053144268618619D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 10) = cmplx( 5.640539003364333226D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 11) = cmplx( 9.982648688423080563D-02, 0.000000000000000000D+00,kind=8)
           EROOTS( 12) = cmplx(-1.159101677081605164D+00, 0.000000000000000000D+00,kind=8)

           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx( 2.180739214839947754D+12, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx( 1.536896462124072305D+21, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx( 4.819762585735845118D+21, 0.000000000000000000D+00,kind=8)
           COEFFS(  5) = cmplx(-4.287885194098189926D+20, 0.000000000000000000D+00,kind=8)
           COEFFS(  6) = cmplx(-1.398982408375965948D+22, 0.000000000000000000D+00,kind=8)
           COEFFS(  7) = cmplx(-1.095743846272546007D+22, 0.000000000000000000D+00,kind=8)
           COEFFS(  8) = cmplx( 4.616345464143469347D+21, 0.000000000000000000D+00,kind=8)
           COEFFS(  9) = cmplx( 5.478694646332004499D+21, 0.000000000000000000D+00,kind=8)
           COEFFS( 10) = cmplx(-6.016228610670696858D+20, 0.000000000000000000D+00,kind=8)
           COEFFS( 11) = cmplx(-6.476339169892343480D+20, 0.000000000000000000D+00,kind=8)
           COEFFS( 12) = cmplx( 1.463412202439570883D+20, 0.000000000000000000D+00,kind=8)
           COEFFS( 13) = cmplx(-8.133965147212722176D+18, 0.000000000000000000D+00,kind=8)

        case (16)  ! Vanni Noferini's example
           print*, "Vanni Noferini's example degree 35"
           EROOTS(  1) = cmplx(-1.291692558483668089D+09, 0.000000000000000000D+00,kind=8)
           EROOTS(  2) = cmplx( 1.055688140578689087D+12, 0.000000000000000000D+00,kind=8)
           EROOTS(  3) = cmplx( 1.675039009280836266D+00, 0.000000000000000000D+00,kind=8)
           EROOTS(  4) = cmplx(-1.735618896149361268D+00, 0.000000000000000000D+00,kind=8)
           EROOTS(  5) = cmplx(-8.736248137695643567D-01, 0.000000000000000000D+00,kind=8)
           EROOTS(  6) = cmplx( 1.567247101836536094D-01, 0.000000000000000000D+00,kind=8)
           EROOTS(  7) = cmplx(-8.123807295252576111D-02, 0.000000000000000000D+00,kind=8)
           EROOTS(  8) = cmplx(-1.700107112233965490D+00, 0.000000000000000000D+00,kind=8)
           EROOTS(  9) = cmplx( 1.629659881926279530D+00, 0.000000000000000000D+00,kind=8)
           EROOTS( 10) = cmplx( 4.754351337267898048D-02, 0.000000000000000000D+00,kind=8)
           EROOTS( 11) = cmplx(-6.385833288358082616D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 12) = cmplx( 1.580580687553951202D+00, 0.000000000000000000D+00,kind=8)
           EROOTS( 13) = cmplx(-1.179899725058404680D+00, 0.000000000000000000D+00,kind=8)
           EROOTS( 14) = cmplx(-1.296753657684092564D+00, 0.000000000000000000D+00,kind=8)
           EROOTS( 15) = cmplx(-3.873421536428965362D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 16) = cmplx( 5.507548547376082126D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 17) = cmplx(-6.857497741235893951D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 18) = cmplx( 4.313897594242247641D-02, 0.000000000000000000D+00,kind=8)
           EROOTS( 19) = cmplx(-1.165802346145184964D+00, 0.000000000000000000D+00,kind=8)
           EROOTS( 20) = cmplx(-2.201353048373792642D+00, 0.000000000000000000D+00,kind=8)
           EROOTS( 21) = cmplx( 1.549696913065910842D+00, 0.000000000000000000D+00,kind=8)
           EROOTS( 22) = cmplx(-1.197850949097038153D+00, 0.000000000000000000D+00,kind=8)
           EROOTS( 23) = cmplx( 1.948799513303653308D+00, 0.000000000000000000D+00,kind=8)
           EROOTS( 24) = cmplx(-1.591087298555633689D+00, 0.000000000000000000D+00,kind=8)
           EROOTS( 25) = cmplx(-1.681364135528388148D+00, 0.000000000000000000D+00,kind=8)
           EROOTS( 26) = cmplx(-9.057681844549185790D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 27) = cmplx(-6.522467879996737272D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 28) = cmplx( 1.536844931675364623D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 29) = cmplx(-9.619673309645394577D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 30) = cmplx( 7.991892105919597977D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 31) = cmplx(-1.469823840934013415D+00, 0.000000000000000000D+00,kind=8)
           EROOTS( 32) = cmplx(-7.001492492072077800D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 33) = cmplx(-4.611880142571331276D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 34) = cmplx(-6.726581586759780018D-01, 0.000000000000000000D+00,kind=8)
           EROOTS( 35) = cmplx( 9.285560636708366511D-01, 0.000000000000000000D+00,kind=8)
           eigsknown = .TRUE.
           
           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx(-1.054396448009028564D+12, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx(-1.363624527049741042D+21, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx(-1.524097086752196185D+22, 0.000000000000000000D+00,kind=8)
           COEFFS(  5) = cmplx(-5.430331895826900792D+22, 0.000000000000000000D+00,kind=8)
           COEFFS(  6) = cmplx( 1.939601603024919514D+22, 0.000000000000000000D+00,kind=8)
           COEFFS(  7) = cmplx( 6.382421578453396740D+23, 0.000000000000000000D+00,kind=8)
           COEFFS(  8) = cmplx( 1.387099440741049132D+24, 0.000000000000000000D+00,kind=8)
           COEFFS(  9) = cmplx(-1.190491576707463537D+24, 0.000000000000000000D+00,kind=8)
           COEFFS( 10) = cmplx(-9.200565166506923742D+24, 0.000000000000000000D+00,kind=8)
           COEFFS( 11) = cmplx(-9.735270106399667711D+24, 0.000000000000000000D+00,kind=8)
           COEFFS( 12) = cmplx( 1.777363646001393349D+25, 0.000000000000000000D+00,kind=8)
           COEFFS( 13) = cmplx( 5.168204627563281514D+25, 0.000000000000000000D+00,kind=8)
           COEFFS( 14) = cmplx( 1.967380995894478195D+25, 0.000000000000000000D+00,kind=8)
           COEFFS( 15) = cmplx(-8.420111546546592614D+25, 0.000000000000000000D+00,kind=8)
           COEFFS( 16) = cmplx(-1.232468300699993727D+26, 0.000000000000000000D+00,kind=8)
           COEFFS( 17) = cmplx( 2.433691458775549011D+24, 0.000000000000000000D+00,kind=8)
           COEFFS( 18) = cmplx( 1.559414230336801483D+26, 0.000000000000000000D+00,kind=8)
           COEFFS( 19) = cmplx( 1.328754831802398326D+26, 0.000000000000000000D+00,kind=8)
           COEFFS( 20) = cmplx(-2.923746138127169960D+25, 0.000000000000000000D+00,kind=8)
           COEFFS( 21) = cmplx(-1.195013474572605217D+26, 0.000000000000000000D+00,kind=8)
           COEFFS( 22) = cmplx(-6.625916327707706196D+25, 0.000000000000000000D+00,kind=8)
           COEFFS( 23) = cmplx( 1.655173557000522824D+25, 0.000000000000000000D+00,kind=8)
           COEFFS( 24) = cmplx( 3.693093358048776752D+25, 0.000000000000000000D+00,kind=8)
           COEFFS( 25) = cmplx( 1.502676482564606996D+25, 0.000000000000000000D+00,kind=8)
           COEFFS( 26) = cmplx(-1.890306739562183708D+24, 0.000000000000000000D+00,kind=8)
           COEFFS( 27) = cmplx(-3.782276693804713849D+24, 0.000000000000000000D+00,kind=8)
           COEFFS( 28) = cmplx(-1.226808372999457977D+24, 0.000000000000000000D+00,kind=8)
           COEFFS( 29) = cmplx(-4.247883694827938356D+22, 0.000000000000000000D+00,kind=8)
           COEFFS( 30) = cmplx( 5.982240413017354273D+22, 0.000000000000000000D+00,kind=8)
           COEFFS( 31) = cmplx( 9.855623614277378638D+21, 0.000000000000000000D+00,kind=8)
           COEFFS( 32) = cmplx(-7.935576108985024512D+20, 0.000000000000000000D+00,kind=8)
           COEFFS( 33) = cmplx(-2.160544389022461460D+20, 0.000000000000000000D+00,kind=8)
           COEFFS( 34) = cmplx( 5.489320495834040320D+18, 0.000000000000000000D+00,kind=8)
           COEFFS( 35) = cmplx( 8.609523091780359680D+17, 0.000000000000000000D+00,kind=8)
           COEFFS( 36) = cmplx(-2.909234290103115200D+16, 0.000000000000000000D+00,kind=8)         

        case (17)  
           print*, "cubic rcoeffsnomial small roots"
           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx(-1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx(-1.000000000000000102D-16, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx( 1.000000000000000102D-16, 0.000000000000000000D+00,kind=8)
           EROOTS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           EROOTS(  2) = cmplx( 1.000000000000000021D-08, 0.000000000000000000D+00,kind=8)
           EROOTS(  3) = cmplx(-1.000000000000000021D-08, 0.000000000000000000D+00,kind=8)
           eigsknown = .TRUE.


        case (18)  
           print*, "cubic rcoeffsnomial very small roots"
           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx(-1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx(-1.000000000000000083D-30, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx( 1.000000000000000083D-30, 0.000000000000000000D+00,kind=8)
           EROOTS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           EROOTS(  2) = cmplx( 1.000000000000000078D-15, 0.000000000000000000D+00,kind=8)
           EROOTS(  3) = cmplx(-1.000000000000000078D-15, 0.000000000000000000D+00,kind=8)
           eigsknown = .TRUE.

        case (19)  
           print*, "cubic rcoeffsnomial large roots"
           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx(-1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx(-1.000000000000000000D+16, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx( 1.000000000000000000D+16, 0.000000000000000000D+00,kind=8)
           EROOTS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           EROOTS(  2) = cmplx( 1.000000000000000000D+08, 0.000000000000000000D+00,kind=8)
           EROOTS(  3) = cmplx(-1.000000000000000000D+08, 0.000000000000000000D+00,kind=8)
           eigsknown = .TRUE.

        case (20)  
           print*, "cubic rcoeffsnomial very large roots"
           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx(-1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx(-1.000000000000000020D+30, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx( 1.000000000000000020D+30, 0.000000000000000000D+00,kind=8)
           EROOTS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           EROOTS(  2) = cmplx( 1.000000000000000000D+15, 0.000000000000000000D+00,kind=8)
           EROOTS(  3) = cmplx(-1.000000000000000000D+15, 0.000000000000000000D+00,kind=8)
           eigsknown = .TRUE.

        case (21)  
           print*, "underflow test"
           COEFFS(  1) = cmplx( 1.000000000000000000D+00, 0.000000000000000000D+00,kind=8)
           COEFFS(  2) = cmplx(-1.111111111000000068D-01, 0.000000000000000000D+00,kind=8)
           COEFFS(  3) = cmplx( 1.122334454433221093D-03, 0.000000000000000000D+00,kind=8)
           COEFFS(  4) = cmplx(-1.123457901110987560D-06, 0.000000000000000000D+00,kind=8)
           COEFFS(  5) = cmplx( 1.123570145779775593D-10, 0.000000000000000000D+00,kind=8)
           COEFFS(  6) = cmplx(-1.123580258012209910D-15, 0.000000000000000000D+00,kind=8)
           COEFFS(  7) = cmplx( 1.123570145779775695D-21, 0.000000000000000000D+00,kind=8)
           COEFFS(  8) = cmplx(-1.123457901110987674D-28, 0.000000000000000000D+00,kind=8)
           COEFFS(  9) = cmplx( 1.122334454433221294D-36, 0.000000000000000000D+00,kind=8)
           COEFFS( 10) = cmplx(-1.111111111000000245D-45, 0.000000000000000000D+00,kind=8)
           COEFFS( 11) = cmplx( 1.000000000000000176D-55, 0.000000000000000000D+00,kind=8)
           EROOTS(  1) = cmplx( 1.000000000000000056D-01, 0.000000000000000000D+00,kind=8)
           EROOTS(  2) = cmplx( 1.000000000000000021D-02, 0.000000000000000000D+00,kind=8)
           EROOTS(  3) = cmplx( 1.000000000000000021D-03, 0.000000000000000000D+00,kind=8)
           EROOTS(  4) = cmplx( 1.000000000000000048D-04, 0.000000000000000000D+00,kind=8)
           EROOTS(  5) = cmplx( 1.000000000000000082D-05, 0.000000000000000000D+00,kind=8)
           EROOTS(  6) = cmplx( 9.999999999999999547D-07, 0.000000000000000000D+00,kind=8)
           EROOTS(  7) = cmplx( 9.999999999999999547D-08, 0.000000000000000000D+00,kind=8)
           EROOTS(  8) = cmplx( 1.000000000000000021D-08, 0.000000000000000000D+00,kind=8)
           EROOTS(  9) = cmplx( 1.000000000000000062D-09, 0.000000000000000000D+00,kind=8)
           EROOTS( 10) = cmplx( 1.000000000000000036D-10, 0.000000000000000000D+00,kind=8)
           eigsknown = .TRUE.
           
        case (22)
           print*, "underflow test"
           COEFFS(  1) = cmplx(     1.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  2) = cmplx(    -1.111111111111111188D-001,     0.000000000000000000D+000,kind=8)
           COEFFS(  3) = cmplx(     1.122334455667789097D-003,     0.000000000000000000D+000,kind=8)
           COEFFS(  4) = cmplx(    -1.123457913581370410D-006,     0.000000000000000000D+000,kind=8)
           COEFFS(  5) = cmplx(     1.123570270608431373D-010,     0.000000000000000000D+000,kind=8)
           COEFFS(  6) = cmplx(    -1.123581506423495493D-015,     0.000000000000000000D+000,kind=8)
           COEFFS(  7) = cmplx(     1.123582630006124464D-021,     0.000000000000000000D+000,kind=8)
           COEFFS(  8) = cmplx(    -1.123582742364387289D-028,     0.000000000000000000D+000,kind=8)
           COEFFS(  9) = cmplx(     1.123582753600102633D-036,     0.000000000000000000D+000,kind=8)
           COEFFS( 10) = cmplx(    -1.123582754722561813D-045,     0.000000000000000000D+000,kind=8)
           COEFFS( 11) = cmplx(     1.123582754823684361D-055,     0.000000000000000000D+000,kind=8)
           COEFFS( 12) = cmplx(    -1.123582754722561718D-066,     0.000000000000000000D+000,kind=8)
           COEFFS( 13) = cmplx(     1.123582753600102447D-078,     0.000000000000000000D+000,kind=8)
           COEFFS( 14) = cmplx(    -1.123582742364387387D-091,     0.000000000000000000D+000,kind=8)
           COEFFS( 15) = cmplx(     1.123582630006124457D-105,     0.000000000000000000D+000,kind=8)
           COEFFS( 16) = cmplx(    -1.123581506423495557D-120,     0.000000000000000000D+000,kind=8)
           COEFFS( 17) = cmplx(     1.123570270608431512D-136,     0.000000000000000000D+000,kind=8)
           COEFFS( 18) = cmplx(    -1.123457913581370721D-153,     0.000000000000000000D+000,kind=8)
           COEFFS( 19) = cmplx(     1.122334455667789376D-171,     0.000000000000000000D+000,kind=8)
           COEFFS( 20) = cmplx(    -1.111111111111111547D-190,     0.000000000000000000D+000,kind=8)
           COEFFS( 21) = cmplx(     1.000000000000000382D-210,     0.000000000000000000D+000,kind=8)
           EROOTS(  1) = cmplx(     1.000000000000000056D-001,     0.000000000000000000D+000,kind=8)
           EROOTS(  2) = cmplx(     1.000000000000000021D-002,     0.000000000000000000D+000,kind=8)
           EROOTS(  3) = cmplx(     1.000000000000000021D-003,     0.000000000000000000D+000,kind=8)
           EROOTS(  4) = cmplx(     1.000000000000000048D-004,     0.000000000000000000D+000,kind=8)
           EROOTS(  5) = cmplx(     1.000000000000000082D-005,     0.000000000000000000D+000,kind=8)
           EROOTS(  6) = cmplx(     9.999999999999999547D-007,     0.000000000000000000D+000,kind=8)
           EROOTS(  7) = cmplx(     9.999999999999999547D-008,     0.000000000000000000D+000,kind=8)
           EROOTS(  8) = cmplx(     1.000000000000000021D-008,     0.000000000000000000D+000,kind=8)
           EROOTS(  9) = cmplx(     1.000000000000000062D-009,     0.000000000000000000D+000,kind=8)
           EROOTS( 10) = cmplx(     1.000000000000000036D-010,     0.000000000000000000D+000,kind=8)
           EROOTS( 11) = cmplx(     9.999999999999999395D-012,     0.000000000000000000D+000,kind=8)
           EROOTS( 12) = cmplx(     9.999999999999999799D-013,     0.000000000000000000D+000,kind=8)
           EROOTS( 13) = cmplx(     1.000000000000000030D-013,     0.000000000000000000D+000,kind=8)
           EROOTS( 14) = cmplx(     9.999999999999999988D-015,     0.000000000000000000D+000,kind=8)
           EROOTS( 15) = cmplx(     1.000000000000000078D-015,     0.000000000000000000D+000,kind=8)
           EROOTS( 16) = cmplx(     9.999999999999999791D-017,     0.000000000000000000D+000,kind=8)
           EROOTS( 17) = cmplx(     1.000000000000000072D-017,     0.000000000000000000D+000,kind=8)
           EROOTS( 18) = cmplx(     1.000000000000000072D-018,     0.000000000000000000D+000,kind=8)
           EROOTS( 19) = cmplx(     9.999999999999999752D-020,     0.000000000000000000D+000,kind=8)
           EROOTS( 20) = cmplx(     9.999999999999999452D-021,     0.000000000000000000D+000,kind=8)
           eigsknown = .TRUE.


        case (23)  
           print*, "deflation stability test"
           COEFFS(  1) = cmplx(     1.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  2) = cmplx(    -1.000100009999999929D+004,     0.000000000000000000D+000,kind=8)
           COEFFS(  3) = cmplx(     1.000100009999999929D+004,     0.000000000000000000D+000,kind=8)
           COEFFS(  4) = cmplx(    -1.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           EROOTS(  1) = cmplx(     1.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           EROOTS(  2) = cmplx(     1.000000000000000000D+004,     0.000000000000000000D+000,kind=8)
           EROOTS(  3) = cmplx(     1.000000000000000048D-004,     0.000000000000000000D+000,kind=8)
           eigsknown = .TRUE.

        case (24)  
           print*, "deflation stability test"
           COEFFS(  1) = cmplx(     1.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  2) = cmplx(    -1.000000100000010058D+007,     0.000000000000000000D+000,kind=8)
           COEFFS(  3) = cmplx(     1.000000100000010058D+007,     0.000000000000000000D+000,kind=8)
           COEFFS(  4) = cmplx(    -1.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           EROOTS(  1) = cmplx(     1.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           EROOTS(  2) = cmplx(     1.000000000000000000D+007,     0.000000000000000000D+000,kind=8)
           EROOTS(  3) = cmplx(     9.999999999999999547D-008,     0.000000000000000000D+000,kind=8)
           eigsknown = .TRUE.

        case (25)  
           print*, "deflation stability test"
           COEFFS(  1) = cmplx(     1.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  2) = cmplx(    -1.000000000100000000D+010,     0.000000000000000000D+000,kind=8)
           COEFFS(  3) = cmplx(     1.000000000100000000D+010,     0.000000000000000000D+000,kind=8)
           COEFFS(  4) = cmplx(    -1.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           EROOTS(  1) = cmplx(     1.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           EROOTS(  2) = cmplx(     1.000000000000000000D+010,     0.000000000000000000D+000,kind=8)
           EROOTS(  3) = cmplx(     1.000000000000000036D-010,     0.000000000000000000D+000,kind=8)
           eigsknown = .TRUE.

        case (26)  
           print*, "deflation stability test M = 15"
           COEFFS(  1) = cmplx(     1.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  2) = cmplx(     2.008113668772817828D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  3) = cmplx(     2.062602533561145313D-001,     0.000000000000000000D+000,kind=8)
           COEFFS(  4) = cmplx(    -2.765068035781857336D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  5) = cmplx(    -2.297742563797626314D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  6) = cmplx(     1.407243980445040066D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  7) = cmplx(     2.658998457037264451D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  8) = cmplx(    -4.334487113945562631D-001,     0.000000000000000000D+000,kind=8)
           COEFFS(  9) = cmplx(    -2.524246561104071862D+000,     0.000000000000000000D+000,kind=8)
           COEFFS( 10) = cmplx(    -1.855934805163346468D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 11) = cmplx(     2.241640989150738505D+000,     0.000000000000000000D+000,kind=8)
           COEFFS( 12) = cmplx(     5.601915780500317243D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 13) = cmplx(    -1.933624888209345638D+000,     0.000000000000000000D+000,kind=8)
           COEFFS( 14) = cmplx(    -7.731250263122397826D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 15) = cmplx(     1.645784111358635116D+000,     0.000000000000000000D+000,kind=8)
           COEFFS( 16) = cmplx(     8.802665861453735641D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 17) = cmplx(    -1.393601953893695677D+000,     0.000000000000000000D+000,kind=8)
           COEFFS( 18) = cmplx(    -9.187425196735220512D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 19) = cmplx(     1.180022353380246969D+000,     0.000000000000000000D+000,kind=8)
           COEFFS( 20) = cmplx(     9.132397847613086261D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 21) = cmplx(    -1.002785735053218108D+000,     0.000000000000000000D+000,kind=8)
           COEFFS( 22) = cmplx(    -8.802134896365174654D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 23) = cmplx(     8.576837238433498989D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 24) = cmplx(     8.306303014588767297D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 25) = cmplx(    -7.400370187436519087D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 26) = cmplx(    -7.717651854435200187D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 27) = cmplx(     6.453538130763551983D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 28) = cmplx(     7.083910365342065418D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 29) = cmplx(    -5.695972366284772770D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 30) = cmplx(    -6.435751942125034208D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 31) = cmplx(     5.092640848134358800D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 32) = cmplx(     5.792176747912534118D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 33) = cmplx(    -4.613737616690670174D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 34) = cmplx(    -5.164170656334372556D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 35) = cmplx(     4.234166367593968361D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 36) = cmplx(     4.557196243525453716D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 37) = cmplx(    -3.932860132781448304D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 38) = cmplx(    -3.972878982338478027D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 39) = cmplx(     3.692047196652574192D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 40) = cmplx(     3.410127405793773048D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 41) = cmplx(    -3.496497658528884278D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 42) = cmplx(    -2.865843212290505315D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 43) = cmplx(     3.332731663024025570D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 44) = cmplx(     2.335320215791015508D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 45) = cmplx(    -3.188114266405422947D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 46) = cmplx(    -1.812390839665631581D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 47) = cmplx(     3.049671184839135019D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 48) = cmplx(     1.289354654036876935D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 49) = cmplx(    -2.902267223969822729D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 50) = cmplx(    -7.567357554348332072D-002,     0.000000000000000000D+000,kind=8)
           COEFFS( 51) = cmplx(     2.725312122540412552D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 52) = cmplx(     2.030745102749671738D-002,     0.000000000000000000D+000,kind=8)
           COEFFS( 53) = cmplx(    -2.485804563312331850D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 54) = cmplx(     3.841629107132993326D-002,     0.000000000000000000D+000,kind=8)
           COEFFS( 55) = cmplx(     2.120988490647549574D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 56) = cmplx(    -1.010257621752959356D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 57) = cmplx(    -1.484590286505819978D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 58) = cmplx(     1.607879790762276884D-001,     0.000000000000000000D+000,kind=8)
           COEFFS( 59) = cmplx(     1.079458153815517421D-002,     0.000000000000000000D+000,kind=8)
           COEFFS( 60) = cmplx(    -9.458473818619324291D-002,     0.000000000000000000D+000,kind=8)
           COEFFS( 61) = cmplx(     4.239115827521624386D-002,     0.000000000000000000D+000,kind=8)
           EROOTS(  1) = cmplx(     9.407561694088811821D-002,    -8.950697058314459609D-001,kind=8)
           EROOTS(  2) = cmplx(     1.871205217359833139D-001,    -8.803328406604251644D-001,kind=8)
           EROOTS(  3) = cmplx(     2.781152949374527394D-001,    -8.559508646656381892D-001,kind=8)
           EROOTS(  4) = cmplx(     3.660629787682203595D-001,    -8.221909118783408132D-001,kind=8)
           EROOTS(  5) = cmplx(     4.500000000000001221D-001,    -7.794228634059947591D-001,kind=8)
           EROOTS(  6) = cmplx(     5.290067270632258234D-001,    -7.281152949374527505D-001,kind=8)
           EROOTS(  7) = cmplx(     6.022175457229724804D-001,    -6.688303429296547087D-001,kind=8)
           EROOTS(  8) = cmplx(     6.688303429296548197D-001,    -6.022175457229724804D-001,kind=8)
           EROOTS(  9) = cmplx(     7.281152949374527505D-001,    -5.290067270632258234D-001,kind=8)
           EROOTS( 10) = cmplx(     7.794228634059948702D-001,    -4.499999999999999556D-001,kind=8)
           EROOTS( 11) = cmplx(     8.221909118783408132D-001,    -3.660629787682201375D-001,kind=8)
           EROOTS( 12) = cmplx(     8.559508646656381892D-001,    -2.781152949374526839D-001,kind=8)
           EROOTS( 13) = cmplx(     8.803328406604251644D-001,    -1.871205217359833972D-001,kind=8)
           EROOTS( 14) = cmplx(     8.950697058314459609D-001,    -9.407561694088811821D-002,kind=8)
           EROOTS( 15) = cmplx(     9.000000000000000222D-001,     0.000000000000000000D+000,kind=8)
           EROOTS( 16) = cmplx(     9.000000000000000222D-001,    -0.000000000000000000D+000,kind=8)
           EROOTS( 17) = cmplx(     8.950697058314459609D-001,     9.407561694088811821D-002,kind=8)
           EROOTS( 18) = cmplx(     8.803328406604251644D-001,     1.871205217359833972D-001,kind=8)
           EROOTS( 19) = cmplx(     8.559508646656381892D-001,     2.781152949374526839D-001,kind=8)
           EROOTS( 20) = cmplx(     8.221909118783408132D-001,     3.660629787682201375D-001,kind=8)
           EROOTS( 21) = cmplx(     7.794228634059948702D-001,     4.499999999999999556D-001,kind=8)
           EROOTS( 22) = cmplx(     7.281152949374527505D-001,     5.290067270632258234D-001,kind=8)
           EROOTS( 23) = cmplx(     6.688303429296548197D-001,     6.022175457229724804D-001,kind=8)
           EROOTS( 24) = cmplx(     6.022175457229724804D-001,     6.688303429296547087D-001,kind=8)
           EROOTS( 25) = cmplx(     5.290067270632258234D-001,     7.281152949374527505D-001,kind=8)
           EROOTS( 26) = cmplx(     4.500000000000001221D-001,     7.794228634059947591D-001,kind=8)
           EROOTS( 27) = cmplx(     3.660629787682203595D-001,     8.221909118783408132D-001,kind=8)
           EROOTS( 28) = cmplx(     2.781152949374527394D-001,     8.559508646656381892D-001,kind=8)
           EROOTS( 29) = cmplx(     1.871205217359833139D-001,     8.803328406604251644D-001,kind=8)
           EROOTS( 30) = cmplx(     9.407561694088811821D-002,     8.950697058314459609D-001,kind=8)
           EROOTS( 31) = cmplx(    -1.045284632676533321D-001,     9.945218953682734009D-001,kind=8)
           EROOTS( 32) = cmplx(    -2.079116908177593426D-001,     9.781476007338056888D-001,kind=8)
           EROOTS( 33) = cmplx(    -3.090169943749473402D-001,     9.510565162951536422D-001,kind=8)
           EROOTS( 34) = cmplx(    -4.067366430758000417D-001,     9.135454576426009776D-001,kind=8)
           EROOTS( 35) = cmplx(    -4.999999999999997780D-001,     8.660254037844387076D-001,kind=8)
           EROOTS( 36) = cmplx(    -5.877852522924730261D-001,     8.090169943749474513D-001,kind=8)
           EROOTS( 37) = cmplx(    -6.691306063588579045D-001,     7.431448254773944662D-001,kind=8)
           EROOTS( 38) = cmplx(    -7.431448254773940221D-001,     6.691306063588583486D-001,kind=8)
           EROOTS( 39) = cmplx(    -8.090169943749473402D-001,     5.877852522924732481D-001,kind=8)
           EROOTS( 40) = cmplx(    -8.660254037844387076D-001,     4.999999999999999445D-001,kind=8)
           EROOTS( 41) = cmplx(    -9.135454576426009776D-001,     4.067366430758000417D-001,kind=8)
           EROOTS( 42) = cmplx(    -9.510565162951535312D-001,     3.090169943749475068D-001,kind=8)
           EROOTS( 43) = cmplx(    -9.781476007338056888D-001,     2.079116908177593148D-001,kind=8)
           EROOTS( 44) = cmplx(    -9.945218953682734009D-001,     1.045284632676532904D-001,kind=8)
           EROOTS( 45) = cmplx(    -1.000000000000000000D+000,     5.665538897647979615D-016,kind=8)
           EROOTS( 46) = cmplx(    -1.000000000000000000D+000,    -5.665538897647979615D-016,kind=8)
           EROOTS( 47) = cmplx(    -9.945218953682734009D-001,    -1.045284632676532904D-001,kind=8)
           EROOTS( 48) = cmplx(    -9.781476007338056888D-001,    -2.079116908177593148D-001,kind=8)
           EROOTS( 49) = cmplx(    -9.510565162951535312D-001,    -3.090169943749475068D-001,kind=8)
           EROOTS( 50) = cmplx(    -9.135454576426009776D-001,    -4.067366430758000417D-001,kind=8)
           EROOTS( 51) = cmplx(    -8.660254037844387076D-001,    -4.999999999999999445D-001,kind=8)
           EROOTS( 52) = cmplx(    -8.090169943749473402D-001,    -5.877852522924732481D-001,kind=8)
           EROOTS( 53) = cmplx(    -7.431448254773940221D-001,    -6.691306063588583486D-001,kind=8)
           EROOTS( 54) = cmplx(    -6.691306063588579045D-001,    -7.431448254773944662D-001,kind=8)
           EROOTS( 55) = cmplx(    -5.877852522924730261D-001,    -8.090169943749474513D-001,kind=8)
           EROOTS( 56) = cmplx(    -4.999999999999997780D-001,    -8.660254037844387076D-001,kind=8)
           EROOTS( 57) = cmplx(    -4.067366430758000417D-001,    -9.135454576426009776D-001,kind=8)
           EROOTS( 58) = cmplx(    -3.090169943749473402D-001,    -9.510565162951536422D-001,kind=8)
           EROOTS( 59) = cmplx(    -2.079116908177593426D-001,    -9.781476007338056888D-001,kind=8)
           EROOTS( 60) = cmplx(    -1.045284632676533321D-001,    -9.945218953682734009D-001,kind=8)
           eigsknown = .TRUE.

        case (27)  ! Bernoulli rcoeffsnomial of degree 20  
           print*, "Bernoulli rcoeffsnomial of degree 20"
           COEFFS(  1) = cmplx(     1.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  2) = cmplx(    -1.000000000000000000D+001,     0.000000000000000000D+000,kind=8)
           COEFFS(  3) = cmplx(     3.166666666666666785D+001,     0.000000000000000000D+000,kind=8)
           COEFFS(  4) = cmplx(     0.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  5) = cmplx(    -1.615000000000000000D+002,     0.000000000000000000D+000,kind=8)
           COEFFS(  6) = cmplx(     0.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  7) = cmplx(     9.228571428571428896D+002,     0.000000000000000000D+000,kind=8)
           COEFFS(  8) = cmplx(     0.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  9) = cmplx(    -4.199000000000000000D+003,     0.000000000000000000D+000,kind=8)
           COEFFS( 10) = cmplx(     0.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS( 11) = cmplx(     1.399666666666666606D+004,     0.000000000000000000D+000,kind=8)
           COEFFS( 12) = cmplx(     0.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS( 13) = cmplx(    -3.188471428571428623D+004,     0.000000000000000000D+000,kind=8)
           COEFFS( 14) = cmplx(     0.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS( 15) = cmplx(     4.522000000000000000D+004,     0.000000000000000000D+000,kind=8)
           COEFFS( 16) = cmplx(     0.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS( 17) = cmplx(    -3.436150000000000000D+004,     0.000000000000000000D+000,kind=8)
           COEFFS( 18) = cmplx(     0.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS( 19) = cmplx(     1.044452380952380918D+004,     0.000000000000000000D+000,kind=8)
           COEFFS( 20) = cmplx(     0.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS( 21) = cmplx(    -5.291242424242424249D+002,     0.000000000000000000D+000,kind=8)
           eigsknown = .FALSE.

        case (28)   ! p(z) = (20!) sum_{k=0}^{20} z^k/k!
           print*, "p(z) = (20!) sum_{k=0}^{20} z^k/k!"
           COEFFS(  1) = cmplx(     1.000000000000000000D+000,     0.000000000000000000D+000,kind=8)
           COEFFS(  2) = cmplx(     2.000000000000000000D+001,     0.000000000000000000D+000,kind=8)
           COEFFS(  3) = cmplx(     3.800000000000000000D+002,     0.000000000000000000D+000,kind=8)
           COEFFS(  4) = cmplx(     6.840000000000000000D+003,     0.000000000000000000D+000,kind=8)
           COEFFS(  5) = cmplx(     1.162800000000000000D+005,     0.000000000000000000D+000,kind=8)
           COEFFS(  6) = cmplx(     1.860480000000000000D+006,     0.000000000000000000D+000,kind=8)
           COEFFS(  7) = cmplx(     2.790720000000000000D+007,     0.000000000000000000D+000,kind=8)
           COEFFS(  8) = cmplx(     3.907008000000000000D+008,     0.000000000000000000D+000,kind=8)
           COEFFS(  9) = cmplx(     5.079110400000000000D+009,     0.000000000000000000D+000,kind=8)
           COEFFS( 10) = cmplx(     6.094932480000000000D+010,     0.000000000000000000D+000,kind=8)
           COEFFS( 11) = cmplx(     6.704425728000000000D+011,     0.000000000000000000D+000,kind=8)
           COEFFS( 12) = cmplx(     6.704425728000000000D+012,     0.000000000000000000D+000,kind=8)
           COEFFS( 13) = cmplx(     6.033983155200000000D+013,     0.000000000000000000D+000,kind=8)
           COEFFS( 14) = cmplx(     4.827186524160000000D+014,     0.000000000000000000D+000,kind=8)
           COEFFS( 15) = cmplx(     3.379030623518720000D+015,     0.000000000000000000D+000,kind=8)
           COEFFS( 16) = cmplx(     2.027418266737049600D+016,     0.000000000000000000D+000,kind=8)
           COEFFS( 17) = cmplx(     1.013709176318197760D+017,     0.000000000000000000D+000,kind=8)
           COEFFS( 18) = cmplx(     4.054836705272791040D+017,     0.000000000000000000D+000,kind=8)
           COEFFS( 19) = cmplx(     1.216451011581837312D+018,     0.000000000000000000D+000,kind=8)
           COEFFS( 20) = cmplx(     2.432902023163674624D+018,     0.000000000000000000D+000,kind=8)
           COEFFS( 21) = cmplx(     2.432902023163674624D+018,     0.000000000000000000D+000,kind=8)
           eigsknown = .FALSE.

        case (29,30,31,32,33)  
           print*, "Bevilacqua, Del Corso, Gemignani P1"
           eigsknown = .FALSE.
           COEFFS = cmplx(0d0,0d0,kind=8)
           COEFFS(1) = cmplx(1d0,0d0,kind=8)
           COEFFS(N+1) = cmplx(1d0,0d0,kind=8)
           COEFFS(N/2+1) = cmplx((1d0*N)/(1d0*N+1d0) + (1d0*N+1d0)/(1d0*N),0d0,kind=8)

        case (34,35,36,37,38)  
           print*, "Bevilacqua, Del Corso, Gemignani P2"
           eigsknown = .FALSE.
           COEFFS(1) = cmplx(1d0,0d0,kind=8)
           do ii = 1,n/2-1
              COEFFS(ii+1) = cmplx(1d0/n*(n+ii),0d0,kind=8)
              COEFFS(n-ii+1) = cmplx(1d0/n*(n+ii),0d0,kind=8)
           end do
           COEFFS(n+1)= cmplx(1d0,0d0,kind=8)
           COEFFS(n/2+1) = cmplx((n+1d0)/(1d0*n),0d0,kind=8)

        case (39,40,41,42,43,44,45,46,47,48)
           print*, "Bevilacqua, Del Corso, Gemignani P3"
           eigsknown = .FALSE.
           if (ll<=43) then
              lambda = 0.9d0
           else
              lambda = 0.999d0
           end if
           COEFFS = cmplx(0d0,0d0,kind=8)
           COEFFS(1) = cmplx(1d0,0d0,kind=8)
           COEFFS(N+1) = cmplx(-1d0,0d0,kind=8)
           COEFFS(N-1+1) = cmplx((lambda+1d0)/(1d0-lambda),0d0,kind=8)
           COEFFS(1+1) = cmplx(-(lambda+1d0)/(1d0-lambda),0d0,kind=8)
     end select
     
     
     
     
     
     
         
     ! call roots
     call z_poly_roots(N,COEFFS,ROOTS,RESIDUALS,INFO)
     
     print*, ROOTS

     ! check INFO
     if (INFO.NE.0) then
        !call u_test_failed(__LINE__)
        print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print*, "!   INFO not 0                                   !"
        print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print*, INFO
     end if
     
     ! compute normofp
     normofp = abs(COEFFS(1))**2
     do ii=1,N
        normofp = normofp + abs(COEFFS(ii+1))**2
     end do
     normofp = dsqrt(normofp)

     print*, normofp

     ! maximum residuals
     res = 0d0
     do ii=1,N
        if (RESIDUALS(ii) >= res) then
           res = RESIDUALS(ii)
        end if
     end do
     
     TABLE_B(kk) = res

     ! maximum forward error
     if (eigsknown) then
        forw = 0d0
        ARGS = 0
        do ii = 1,N
           
           err = EISCOR_DBL_INF
           ll = 0
           do jj = 1,N
              !if ((ARGS(ll).EQ.0).AND.(err .GT. abs(ROOTS(ii)-EROOTS(jj)))) then
              if (err .GT. abs(ROOTS(ii)-EROOTS(jj))) then
                 err = abs(ROOTS(ii)-EROOTS(jj))
                 ll = jj
              end if
           end do
           print*, ii, ll, ROOTS(ii), EROOTS(ll), abs(ROOTS(ii)-EROOTS(ll)), err
           
           ARGS(ll) = ARGS(ll)+1
           
           if (err/abs(EROOTS(ll)) .GT. forw) then
              forw = err/abs(EROOTS(ll))
           end if
        end do
        
        TABLE_F(kk) = forw
     else
        TABLE_F(kk) = -1d0!EISCOR_DBL_INF+1
     end if
     EK(kk) = eigsknown
     
     deallocate(RESIDUALS,FORWARD,COEFFS,EROOTS,ROOTS,ARGS)
  end do
  ! end check 1)


  ! print table
  print*, "Relative Forward Error"
  do kk=1,48
     if (EK(kk)) then
        write (*,"(I3,1x,A,1x,ES10.4E2,1x,A)") kk, "&", TABLE_F(kk), "\\%"
     end if
  end do

  print*, "Relative Backward Error (residual based, not comparable with old paper)"
  do kk=1,48
     write (*,"(I3,1x,A,1x,ES10.4E2,1x,A)") kk, "&", TABLE_B(kk), "\\%"
  end do

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  print*, "Runtime of this example: ", dble(c_stop-c_start)/dble(c_rate)
  print*,""

     
end program example_z_poly_roots
