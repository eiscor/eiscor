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
! check 1) [-0.5 0 -0.5]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_symtrid_qr  

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
  logical :: backward
  
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
    
  call z_unifact_qr(.FALSE.,.FALSE.,N,Q,T,N,D,ITS,INFO)

!  print*,""
!  print*,"Unitary"
!  do ii=1,N
!    print*,U(ii),VV(ii)
!  end do

  call z_urffact_rfqr(N,U,VV,ITS,INFO)

  ! transform eigenvalues
  do ii=1,N
     D(ii) = aimag(U(ii))/(1d0+dble(U(ii)))
     E(ii) = T(2*ii)/(1d0+T(2*ii-1))
  end do

!  ! sort by real eigenvalue
!  do ii = 1,N
!    small = 1d0
!    id = ii
!    do jj = ii,N
!      if ( D(jj) < small ) then
!        id = jj
!        small = D(id)
!      end if
!    end do
!    D(id) = D(ii)
!    D(ii) = small
!    swap = U(id)
!    U(id) = U(ii)
!    U(ii) = swap
!  end do
!    
!  ! sort by real eigenvalue
!  do ii = 1,N
!    small = 1d0
!    id = ii
!    do jj = ii,N
!      if ( E(jj) < small ) then
!        id = jj
!        small = E(id)
!      end if
!    end do
!    E(id) = E(ii)
!    E(ii) = small
!    swap = cmplx(T(2*id-1),T(2*id),kind=8)
!    T(2*id-1) = T(2*ii-1)
!    T(2*id)   = T(2*ii)
!    T(2*ii-1) = dble(swap)
!    T(2*ii)   = aimag(swap)
!  end do
    
  ! exact eigenvalues
  do ii=1,N
     VV(ii) = cos(dble(N+1-ii)*pi/(dble(N)+1d0))
  end do
  
  print*,""
  print*,"Eigs"
  do ii=1,N
    !swap = cmplx(-VV(ii),1d0,kind=8)/cmplx(VV(ii),1d0,kind=8)
    !print*,abs(swap-U(ii)),U(ii),cmplx(T(2*ii-1),T(2*ii),kind=8),abs(swap-cmplx(T(2*ii-1),T(2*ii),kind=8))
    print*,U(ii),cmplx(T(2*ii-1),T(2*ii),kind=8),abs(U(ii)-cmplx(T(2*ii-1),T(2*ii),kind=8))
    !print*,abs(D(ii)-VV(ii)),D(ii),E(ii),abs(VV(ii)-E(ii))
  end do

!  print*,""
!  print*,"Eigs"
!  do ii=1,N
!    print*,U(ii),(cmplx(0d0,1d0,kind=8)-VV(ii))/(cmplx(0d0,1d0,kind=8)+VV(ii)),D(ii),VV(ii),abs(D(ii)-VV(ii))
!  end do

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))  

end program test_d_symtrid_qr
