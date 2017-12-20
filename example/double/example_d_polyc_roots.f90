!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_polyc_roots
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the roots of a polynomial in Chebyshev basis.
! The polynomial is of dimension N.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_polyc_roots

  implicit none
  
  ! compute variables
  integer, parameter :: N = 10
  integer :: ii, jj, ij
  real(8) :: COEFFS(N), RESIDUALS(N), a, b, RES(N,3)
  complex(8) :: ROOTS(N), E(4096), C(N+3), ac
  complex(8) :: CCOEFFS(N), RECUR(N,3), ALLROOTS (N,1)


  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  print*,""
  print*,"example_d_polyc_roots:"
  print*,""


  select case (N)
  case (3)
     E(1) = cmplx( 0.557453770738378d0,0d0,kind=8)
     E(2) = cmplx(-0.627050844182526d0,0d0,kind=8)
     E(3) = cmplx(-1.430402926555852d0,0d0,kind=8)
  case (4) 
     E(1) = cmplx( 0.762168899481069d0,0d0,kind=8)
     E(2) = cmplx(-0.103856504303369d0,0d0,kind=8)
     E(3) = cmplx(-0.896143495696630d0,0d0,kind=8)
     E(4) = cmplx(-1.762168899481096d0,0d0,kind=8)
  case (6)     
     E(1) = cmplx( 0.902287022664004d0,0d0,kind=8)
     E(2) = cmplx( 0.449789049847388d0,0d0,kind=8)
     E(3) = cmplx(-0.097141724794192d0,0d0,kind=8)
     E(4) = cmplx(-0.961042100503975d0,0d0,kind=8)
     E(5) = cmplx(-0.616008095011604d0,0d0,kind=8)
     E(6) = cmplx(-2.677884152201622d0,0d0,kind=8)  
  case (16)
     E(1)  = cmplx( 0.988204083122132d0,0d0,kind=8)
     E(2)  = cmplx( 0.920275424337706d0,0d0,kind=8)
     E(3)  = cmplx( 0.829340725875603d0,0d0,kind=8)
     E(4)  = cmplx( 0.692506316000414d0,0d0,kind=8)
     E(5)  = cmplx( 0.536555483924795d0,0d0,kind=8)
     E(6)  = cmplx( 0.351115892910315d0,0d0,kind=8)
     E(7)  = cmplx( 0.155268963574267d0,0d0,kind=8)
     E(8)  = cmplx(-0.049030568926040d0,0d0,kind=8)
     E(9)  = cmplx(-0.251674067719672d0,0d0,kind=8)
     E(10) = cmplx(-0.995211608857433d0,0d0,kind=8)
     E(11) = cmplx(-0.954157319555161d0,0d0,kind=8)
     E(12) = cmplx(-0.877096681676023d0,0d0,kind=8)
     E(13) = cmplx(-0.759966324925241d0,0d0,kind=8)
     E(14) = cmplx(-0.441233277211359d0,0d0,kind=8)
     E(15) = cmplx(-0.615854972922036d0,0d0,kind=8)
     E(16) = cmplx(-7.529042067952273d0,0d0,kind=8)
  
  end select


  print*, "Eigenvalues computed with MATLAB"
  do ii = 1,N
     if (N.LE.16) then
        print*, ii, E(ii)
     end if
     COEFFS(ii) = 1d0*(N+1-ii)
     CCOEFFS(ii) = cmplx(COEFFS(ii),0d0,kind=8)
     !COEFFS(ii) = 1d0*(ii)
  end do
  
  call d_polyc_roots(N,COEFFS,ROOTS,RESIDUALS)

  if (N.LE.16) then
     a = 0d0
     do ii=1,N
        b = abs(E(1)-ROOTS(ii))
        ij = 1
        do jj = 2,N
           if (abs(E(jj)-ROOTS(ii)).LT.b) then
              b = abs(E(jj)-ROOTS(ii))
              ij = jj
           end if
        end do
        print*, ROOTS(ii), E(ij), abs(E(ij)-ROOTS(ii))
        a = a + abs(E(ij)-ROOTS(ii))
        !ROOTS(ii) = E(ii)
     end do
     
     print*, "sum of errors", a

  ! check polynomial value in roots
  ! using Clenshaws algorithm
     a = 0d0
     do ii=1,N
        C = cmplx(0d0,0d0,kind=8)
        C(N+1) = cmplx(1d0,0d0,kind=8)
        do jj=N-1,1,-1
           C(jj+1) = cmplx(COEFFS(N+1-(jj+1)),0d0,kind=8) + cmplx(2d0,0d0,kind=8)*ROOTS(ii)*C(jj+2) - C(jj+3)
        end do
        C(1) = cmplx(2d0*COEFFS(N+1-1),0d0,kind=8) + cmplx(2d0,0d0,kind=8)*ROOTS(ii)*C(2) - C(3)
        
        ac = (C(1) - C(3))/cmplx(2d0,0d0,kind=8)
        print*, ROOTS(ii),ac, abs(ac)
        if (abs(ac).GT.a) then
           a = abs(ac)
        end if
     end do
     print*, "maximal polynomial value in roots", a
     
     print*, "Real Roots"
     ! using Clenshaws algorithm
     do ii=1,N
        C = cmplx(0d0,0d0,kind=8)
        C(N+1) = cmplx(1d0,0d0,kind=8)
        do jj=N-1,1,-1
           C(jj+1) = cmplx(COEFFS(N+1-(jj+1)),0d0,kind=8) + cmplx(2d0,0d0,kind=8)*E(ii)*C(jj+2) - C(jj+3)
        end do
        C(1) = cmplx(2d0*COEFFS(N+1-1),0d0,kind=8) + cmplx(2d0,0d0,kind=8)*E(ii)*C(2) - C(3)
        
        ac = (C(1) - C(3))/cmplx(2d0,0d0,kind=8)
        print*, E(ii),ac, abs(ac)
     end do
     
     RECUR = cmplx(0d0,0d0,kind=8)
     RECUR(:,1) = cmplx(.5d0,0d0,kind=8)
     RECUR(:,3) = cmplx(.5d0,0d0,kind=8)
     RECUR(N,1) = cmplx(1d0,0d0,kind=8)
     RECUR(N,3) = cmplx(0d0,0d0,kind=8)
     
     
     call z_polyc_residuals(N,3,0,CCOEFFS,RECUR,ROOTS,ALLROOTS,RES)

     print*, RES(:,1)

     b = 0d0
     do ii=1,N
        b = b + abs(RES(ii,1))
     end do
     
     print*, "sum of residuals", b

  end if

  
  ! stop timer
  call system_clock(count=c_stop)
  print*, "Test took ", dble(c_stop-c_start)/dble(c_rate), "s"

end program example_d_polyc_roots
