#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_symtrid_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the eigenvalues of the tridiagonal matrix
! [-0.5 0 -0.5] testing the forward error and of a random normally 
! distributed matrix testing the backward error.  
!
! check 1) [-0.5 0 -0.5] plus random spike
! check 2) random tridiagonal matrix + random spike
! check 3) random tridiagonal matrix with first row/column zero
! check 4) random tridiagonal matrix with last row/column zero
! check 5) random tridiagonal matrix with row 5 and column 5 zero
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_symtrid_qr  

  implicit none
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! parameters
  integer, parameter :: N = 4096
  logical, parameter :: sca = .FALSE.
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute variables
  integer :: M
  integer :: ii, jj, kk, INFO
  real(8) :: WORK(5*N), D(N), E(N-1), eig(N), t, t1, t2, t3, nrm
  real(8) :: Ds(N), Es(N-1), pi = EISCOR_DBL_PI
  complex(8) :: Z(N,N), v(N), U(N), Us(N)
  integer :: ITS(N-1)
  logical :: backward
  !
  logical :: P(N-2)
  real(8) :: Q(3*(N-1)), C1(3*N), B1(3*N), D1(2*(N+1))
  real(8) :: D2(2*(N+1)), W, B2(3*N), C2(3*N), LWORK(2*N)
  complex(8) :: ROOTS(N)
  ! RES
  complex(8) :: CCOEFFS(N), RECUR(N,3), ALLROOTS (N,1)
  real (8) :: RES(N,3)
  ! 
  complex(8) :: H(N,N), HT(N,N), S(N,N)
  ! timing variables
  integer:: c_start, c_stop, c_rate


  interface
    function l_upr1fact_hess(m,flags)
      logical :: l_upr1fact_hess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_hess
  end interface


  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)
  
  
  do kk=1,5
     call u_fixedseed_initialize(INFO)
     do ii=1,50*kk*N
        call random_number(t)
     end do
     CCOEFFS = cmplx(0d0,0d0,kind=8)
     select case (kk)
     case (1)
        ! check 1) [-0.5 0 -0.5] plus spike
        ! initialize T to be a tridiagonal matrix of the form
        !1   0 -1
        !-  -1  0 -1
        !2     -1  0 ...
        Ds = 0d0
        Es = 5d-1
        Es(1) = sqrt(2d0)/2d0
        !call z_1Darray_random_normal(N,Us)
        Us = 0d0
        Us(1) = cmplx(sqrt(2d0)/2d0,0d0,kind=8)
        CCOEFFS(N) = cmplx(-1d0,0d0,kind=8)
        
     case (2)        
        CCOEFFS = cmplx(0d0,0d0,kind=8)
        ! check 2) random tridiagonal matrix + random spike
        do ii=1,N-1
           call random_number(t)
           Ds(ii) = t
           call random_number(t)
           Es(ii) = t
           call random_number(t)
           Us(ii) = cmplx(t,0d0,kind=8)
           CCOEFFS(N+1-ii) = -Us(ii)
        end do
        Us(1) = Us(1)/sqrt(2d0)
        call random_number(t)
        Ds(N) = t
        call random_number(t)
        Us(N) = cmplx(t,0d0,kind=8)
        CCOEFFS(1) = -Us(N)

     case (3)        
        CCOEFFS = cmplx(0d0,0d0,kind=8)
        ! check 3) random tridiagonal matrix with first row/column zero
        do ii=1,N-1
           call random_number(t)
           Ds(ii) = t
           call random_number(t)
           Es(ii) = t
           call random_number(t)
           Us(ii) = cmplx(t,0d0,kind=8)
           CCOEFFS(N+1-ii) = -Us(ii)
        end do
        Us(1) = Us(1)/sqrt(2d0)
        call random_number(t)
        Ds(N) = t
        call random_number(t)
        Us(N) = cmplx(t,0d0,kind=8)
        CCOEFFS(1) = -Us(N)
        !Ds(1) = 0d0
        !Es(1) = 0d0
        Us(1) = 0d0

     case (4)        
        CCOEFFS = cmplx(0d0,0d0,kind=8)
        ! check 4) random tridiagonal matrix with last row/column zero
        do ii=1,N-1
           call random_number(t)
           Ds(ii) = t
           call random_number(t)
           Es(ii) = t
           call random_number(t)
           Us(ii) = cmplx(t,0d0,kind=8)
           CCOEFFS(N+1-ii) = -Us(ii)
        end do
        Us(1) = Us(1)/sqrt(2d0)
        call random_number(t)
        Ds(N) = t
        call random_number(t)
        Us(N) = cmplx(t,0d0,kind=8)
        CCOEFFS(1) = -Us(N)
        !Ds(N) = 0d0
        !Es(N-1) = 0d0
        Us(N) = 0d0


     case (5)        
        CCOEFFS = cmplx(0d0,0d0,kind=8)
        ! check 5) random tridiagonal matrix with row 5 and column 5 zero
        do ii=1,N-1
           call random_number(t)
           Ds(ii) = t
           call random_number(t)
           Es(ii) = t
           call random_number(t)
           Us(ii) = cmplx(t,0d0,kind=8)
           CCOEFFS(N+1-ii) = -Us(ii)
        end do
        Us(1) = Us(1)/sqrt(2d0)
        call random_number(t)
        Ds(N) = t
        call random_number(t)
        Us(N) = cmplx(t,0d0,kind=8)
        CCOEFFS(1) = -Us(N)
        !Ds(5) = 0d0
        !Es(5) = 0d0
        Us(5) = 0d0
     end select

  U = Us
  D = Ds
  E = Es

  ! fill P
  P = .FALSE.

  call d_spr1_factor(.TRUE.,.TRUE.,.FALSE.,N,D,E,U,Q,D1,C1,B1,t2,N,Z,LWORK,INFO)

  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if

  call z_upr1fact_twistedqz(.FALSE.,.TRUE.,.FALSE.,l_upr1fact_hess,N,P,Q,D1,C1,B1,D2,C2,B2,Z,W,ITS,INFO)

  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if


  ! extract roots
  call z_upr1fact_extracttri(.TRUE.,N,D1,C1,B1,ROOTS)
  ! back transformation
  do ii=1,N
     ROOTS(ii) = cmplx(0d0,1d0,kind=8)*(cmplx(1d0,0d0,kind=8)-ROOTS(ii))/(cmplx(1d0,0d0,kind=8)+ROOTS(ii))
     !print*, ii, ROOTS(ii)
  end do
  

  if (kk.EQ.1) then
     RECUR = cmplx(0d0,0d0,kind=8)
     RECUR(:,1) = cmplx(.5d0,0d0,kind=8)
     RECUR(:,3) = cmplx(.5d0,0d0,kind=8)
     RECUR(N,1) = cmplx(1d0,0d0,kind=8)
     RECUR(N,3) = cmplx(0d0,0d0,kind=8)
     
     RES = -10d0
     
     call z_polyc_residuals(N,3,0,CCOEFFS,RECUR,ROOTS,ALLROOTS,RES)
     
     !print*, "RESIDUAL", RES(:,3)
     t2 = 0d0
     do ii=1,N
        t2 = t2 + RES(ii,3)
     end do
     
        print*, "case", kk, "norm", t2
     if ((t2>1d-14*N*N).OR.(t2.NE.t2)) then
        ! backward error test failed
        print*, "case", kk, "norm", t2
        call u_test_failed(__LINE__)
     end if     
  end if

  call z_upr1fact_extracttri(.FALSE.,N,D1,C1,B1,HT)

  ! HT is part of the Schur form of the upr1 problem, 
  ! HT is upper triangular, now we back transform HT
  
  H = HT
  S = -HT
  do ii=1,N
     H(ii,ii) = H(ii,ii)+cmplx(1d0,0d0,kind=8)
     S(ii,ii) = S(ii,ii)+cmplx(1d0,0d0,kind=8)
  end do

  ! S = S(H)^-1 = (1-HT)(1+HT)^(-1)
  call ztrsm('R', 'U', 'N', 'N', N, N, cmplx(1d0,0d0,kind=8), H, N, S, N)
  
  ! S = i S
  S = cmplx(0d0,1d0,kind=8)*S

  ! H is original matrix
  H = 0d0
  do ii=1,N-1
     H(ii,ii) = H(ii,ii) + Ds(ii)
     H(ii+1,ii) = H(ii+1,ii) + Es(ii)
     H(ii,ii+1) = H(ii,ii+1) + Es(ii)
     H(ii,N) = H(ii,N) + Us(ii)
  end do
  H(N,N) = H(N,N) + Ds(ii) + Us(ii)

  ! HT = Z^H
  HT = 0d0
  do ii=1,N
     do jj=1,N
        HT(ii,jj) = conjg(Z(jj,ii))
     end do
  end do
  
  H = H - matmul(Z,matmul(S,HT))

   
  ! compute the Frobenius norm of H
  t = 0d0
  do ii=1,N
     do jj=1,N
        t = t + abs(H(ii,jj))**2
     end do
  end do
  t = sqrt(t)

  if ((t>1d-14*N*N).OR.(t.NE.t)) then
     ! backward error test (involving the Schur vectors) failed
     print*, "case", kk, "norm", t
     call u_test_failed(__LINE__)
  end if

  end do

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))  

end program test_d_symtrid_qr
