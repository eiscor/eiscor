#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_symtrid_qr3
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the eigenvalues of the tridiagonal matrix
! [-0.5 0 -0.5] testing the forward error and of a random normally 
! distributed matrix testing the backward error.  
!
! check 1) [-0.5 0 -0.5]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_symtrid_qr3  

  implicit none
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! parameters
  integer, parameter :: N = 2**09
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute variables
  integer :: M
  integer :: ii, jj, id, INFO
  real(8) :: D(N), E(N), VV(N)
  real(8) :: Q(3*(N-1)), T(2*N)
  real(8) :: small, pi = EISCOR_DBL_PI
  complex(8) :: U(N), swap
  integer :: ITS(N-1)
  real(8) :: xx
  complex(8) :: rho, block(2,2), t1(2,2), t2(2,2)
  complex(8) :: H1(N), H2(N), U2(N)
  real(kind=8) :: tV(N)
  complex(kind=8) :: tU(N)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)

  ! check 1) [-0.5 0 -0.5]
  ! initialize T to be a tridiagonal matrix of the form
  !  0 -1
  ! -1  0 -1
  !     -1 0 ...
  D = 0d0
  E = -5d-1

  call d_symtrid_factor4(.FALSE.,N,D,E,U,VV,1d0,INFO)

  print*,""
  do ii=1,N
    tU(ii) = U(ii)
    tV(ii) = sqrt(VV(ii))
!    print*,tU(ii),U(ii)
!    print*,tV(ii),VV(ii)
  end do

  Q = 0d0
  T = 0d0
  do ii=1,N-1
    Q(3*ii-2) = dble(U(ii))
    Q(3*ii-1) = aimag(U(ii))
    Q(3*ii) = sqrt(VV(ii))
    T(2*ii-1) = 1d0
  end do
  T(2*N-1) = dble(U(N))
  T(2*N) = aimag(U(N))
    
  ! H1
  H1(1) = U(1)
  do ii = 2,N
    H1(ii) = conjg(U(ii-1))*U(ii)
  end do

  ! H2
  H2(1) = cmplx(Q(1),Q(2),kind=8)*cmplx(T(1),T(2),kind=8)
  do ii = 2,N-1
    H2(ii) = cmplx(Q(3*(ii-1)-2),-Q(3*(ii-1)-1),kind=8)*cmplx(Q(3*ii-2),Q(3*ii-1),kind=8)*cmplx(T(2*ii-1),T(2*ii),kind=8)
  end do
  H2(N) = cmplx(Q(3*(N-1)-2),-Q(3*(N-1)-1),kind=8)*cmplx(T(2*N-1),T(2*N),kind=8)

!print*,""
!do ii=1,N
!  print*,H1(ii),H2(ii),abs(H1(ii)-H2(ii))
!end do

  ! get 2x2 block
  block(1,1) =  U(N-1)
  block(2,2) =  conjg(U(N-1))
  block(1,2) = -sqrt(VV(N-1))
  block(2,1) =  sqrt(VV(N-1))
  block(:,2) =  block(:,2)*U(N)
  if (N > 2) then
    xx = abs(U(N-2))
    if (xx.GT.0) then
      block(1,:) = conjg(U(N-2))*block(1,:)/xx
    end if
  end if
    
  ! compute eigenvalues and eigenvectors
  t1 = block
  call z_2x2array_eig(.FALSE.,t1,t1,t2,t2)
    
  ! choose wikinson shift
  ! complex abs does not matter here
  if(abs(block(2,2)-t1(1,1)) < abs(block(2,2)-t1(2,2)))then
    rho = t1(1,1)
  else
    rho = t1(2,2)
  end if

  ! compute a nonzero shift
  ! random shift
  xx = abs(rho)
  if (xx == 0) then
    call random_number(xx)
    rho = cmplx(cos(xx),sin(xx),kind=8)
  ! wilkinson shift
  else
    rho = rho/xx
  end if

!  call z_unifact_singlestep_shift(.FALSE.,N,Q,T,N,D,ITS,rho)
  call z_urffact_singlestep_shift(N,U,VV,ITS,rho)
  call z_urffact_singlestep_shift2(N,tU,tV,ITS,rho)

  print*,""
  do ii=1,N
    print*,U(ii),tU(ii),abs(tU(ii)-U(ii))
!    print*,tV(ii),VV(ii)
  end do

!  ! H1
  H1(1) = U(1)
  do ii = 2,N
    H1(ii) = conjg(U(ii-1))*U(ii)
  end do

  ! H2
  H2(1) = cmplx(Q(1),Q(2),kind=8)*cmplx(T(1),T(2),kind=8)
  do ii = 2,N-1
    H2(ii) = cmplx(Q(3*(ii-1)-2),-Q(3*(ii-1)-1),kind=8)*cmplx(Q(3*ii-2),Q(3*ii-1),kind=8)*cmplx(T(2*ii-1),T(2*ii),kind=8)
  end do
  H2(N) = cmplx(Q(3*(N-1)-2),-Q(3*(N-1)-1),kind=8)*cmplx(T(2*N-1),T(2*N),kind=8)

  ! U2
  U2(1) = H2(1)
  do ii = 2,N
    U2(ii) = H2(ii)/conjg(U2(ii-1))
  end do

!print*,""
!print*,rho
!print*,""
!do ii=1,N-1
!  !print*,H1(ii),H2(ii),abs(H1(ii)-H2(ii))
!  print*,abs(U(ii)-tU(ii)),U(ii),U2(ii),abs(tU(ii)-U2(ii))
!!  print*,cmplx(VV(ii),0d0,kind=8),cmplx(Q(3*ii)**2,0d0,kind=8),abs(VV(ii)-Q(3*ii)**2)
!!print*,""
!end do
!print*,U(N),U2(N),abs(U(N)-U2(N))
!print*,""
!do ii=1,N-1
!  print*,VV(ii),Q(3*ii)**2,abs(VV(ii)-Q(3*ii)**2),abs(U(ii))**2+VV(ii),abs(U2(ii))**2+Q(3*ii)**2
!end do


  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))  

end program test_d_symtrid_qr3
