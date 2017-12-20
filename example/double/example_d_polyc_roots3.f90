#include "eiscor.h"
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
  integer, parameter :: NEWTONSTEPS = 1
  !integer, parameter :: N = 5
  !integer, parameter :: N = 15
  !integer, parameter :: N = 20
  !integer, parameter :: N = 35
  !integer, parameter :: N = 50
  !integer, parameter :: N = 100
  !integer, parameter :: N = 1600 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !integer, parameter :: N = 4
  !integer, parameter :: N = 8
  !integer, parameter :: N = 16
  !integer, parameter :: N = 32
  !integer, parameter :: N = 64
  !integer, parameter :: N = 128
  !integer, parameter :: N = 256
  integer, parameter :: N = 512
  !integer, parameter :: N = 1024
  !integer, parameter :: N = 2048
  !integer, parameter :: N = 4096
  !integer, parameter :: N = 8192
  !integer, parameter :: N = 16384
  real(8), parameter :: scale = 1d0
  integer :: ii, jj!, ij
  real(8), allocatable :: COEFFS(:), RESIDUALS(:), RES(:,:)
  real(8) :: a, b
  complex(8), allocatable :: ROOTS(:)!, C(:)
  !complex(8) :: ac
  complex(8), allocatable :: CCOEFFS(:), RECUR(:,:), ALLROOTS (:,:)

  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! allocate memory
  allocate(COEFFS(N),RESIDUALS(N),RES(N,(3*(NEWTONSTEPS+1))),ROOTS(N))
  allocate(CCOEFFS(N))
  allocate(RECUR(N,3),ALLROOTS(N,(NEWTONSTEPS+1)))
  !,C(N+3))

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  print*,""
  print*,"example_d_polyc_roots:"
  print*,""

  COEFFS = 0d0
  COEFFS(N) = 1d0
  RES = -1d0

  do ii = 1,N
     COEFFS(ii) = 1d0*(N+1-ii)
  end do
  !call d_1Darray_random_normal(N,COEFFS)

  !COEFFS(1) = 0d0

  do ii=1,N
     CCOEFFS(ii) = cmplx(COEFFS(ii),0d0,kind=8)
  end do


!!$  if (N.LT.40) then
!!$     print*, COEFFS
!!$     print*, CCOEFFS
!!$  end if
  
  print*, "N = ", N
  
  call d_polyc_roots(N,COEFFS,ROOTS,RESIDUALS)

  RECUR = cmplx(0d0,0d0,kind=8)
  RECUR(:,1) = cmplx(.5d0,0d0,kind=8)
  RECUR(:,3) = cmplx(.5d0,0d0,kind=8)
  RECUR(N,1) = cmplx(1d0,0d0,kind=8)
  RECUR(N,2) = cmplx(0d0,0d0,kind=8)    
  RECUR(N,3) = cmplx(0d0,0d0,kind=8)    

  print*, "residuals"

!!$  ROOTS(15)=cmplx( -4.7807767719522110d-002, 0.99641554167562285d0,kind=8     ) 
!!$  ROOTS(16)=cmplx( -4.4193964153617302d-002, -1.0014695287012816d0,kind=8     )
!!$  ROOTS(17)=cmplx( 0.89240144335611238d0     ,-0.44930917523233938d0,kind=8     ) 
!!$  ROOTS(18)=cmplx(-0.83364271299236237d0     , 0.55388066597595420d0,kind=8     )
!!$  ROOTS(19)=cmplx( -4.8778646827590436E-002, 0.99792329719991701d0,kind=8     ) 
!!$  ROOTS(20)=cmplx( -4.6209996764836697E-002,-0.99981873869128890d0,kind=8     )
!!$
!!$  ROOTS(15) = cmplx(0d0,1d0,kind=8)*(cmplx(1d0,0d0,kind=8)-ROOTS(15))/(cmplx(1d0,0d0,kind=8)+ROOTS(15))
!!$  ROOTS(16) = cmplx(0d0,1d0,kind=8)*(cmplx(1d0,0d0,kind=8)-ROOTS(16))/(cmplx(1d0,0d0,kind=8)+ROOTS(16))
!!$  ROOTS(17) = cmplx(0d0,1d0,kind=8)*(cmplx(1d0,0d0,kind=8)-ROOTS(17))/(cmplx(1d0,0d0,kind=8)+ROOTS(17))
!!$  ROOTS(18) = cmplx(0d0,1d0,kind=8)*(cmplx(1d0,0d0,kind=8)-ROOTS(18))/(cmplx(1d0,0d0,kind=8)+ROOTS(18))
!!$  ROOTS(19) = cmplx(0d0,1d0,kind=8)*(cmplx(1d0,0d0,kind=8)-ROOTS(19))/(cmplx(1d0,0d0,kind=8)+ROOTS(19))
!!$  ROOTS(20) = cmplx(0d0,1d0,kind=8)*(cmplx(1d0,0d0,kind=8)-ROOTS(20))/(cmplx(1d0,0d0,kind=8)+ROOTS(20))

!!$  a = 0d0
!!$  do ii=1,N
!!$     C = cmplx(0d0,0d0,kind=8)
!!$     C(N+1) = cmplx(1d0,0d0,kind=8)
!!$     do jj=N-1,1,-1
!!$        C(jj+1) = cmplx(COEFFS(N+1-(jj+1)),0d0,kind=8) + cmplx(2d0,0d0,kind=8)*ROOTS(ii)*C(jj+2) - C(jj+3)
!!$     end do
!!$     C(1) = cmplx(2d0*COEFFS(N+1-1),0d0,kind=8) + cmplx(2d0,0d0,kind=8)*ROOTS(ii)*C(2) - C(3)
!!$     
!!$     ac = (C(1) - C(3))/cmplx(2d0,0d0,kind=8)
!!$     print*, ii, ROOTS(ii), ac, abs(ac)
!!$     if (abs(ac).GT.a) then
!!$        a = abs(ac)
!!$     end if
!!$  end do
!!$  print*, "maximal polynomial value in roots", a
  


  !do ii=1,N
  !   ROOTS(ii) = cmplx(dble(ROOTS(ii)),0d0,kind=8)
  !end do
  call z_polyc_residuals(N,3,NEWTONSTEPS,CCOEFFS,RECUR,ROOTS,ALLROOTS,RES)
  
!!$  print*, "roots"
!!$  do ii=1,N
!!$     !print*, ii, ROOTS(ii), (ROOTS(ii)-cmplx(0.5d0,0d0,kind=8))*2d0, &
!!$     !     &abs(ROOTS(ii)), abs((ROOTS(ii)-cmplx(0.5d0,0d0,kind=8))*2d0),RES(ii,1)
!!$   
!!$     if (abs(RES(ii,3)).GT.1d-12) then
!!$        print*, ii, ROOTS(ii), RES(ii,3)
!!$     end if
!!$     !ROOTS(ii) = cmplx(dble(ROOTS(ii)),0d0,kind=8)
!!$  end do

  if (N.LT.10) then
     do jj=1,(NEWTONSTEPS+1)
        print*, RES(:,3*jj-2) 
        print*, RES(:,3*jj-1)
        print*, RES(:,3*jj)
     end do
  end if
  
  a = 0d0
  b = 0d0
  do ii=1,N
     if (abs(RES(ii,3)).GT.a) then
        a = abs(RES(ii,3))
        print*, ii, ROOTS(ii), RES(ii,3)
        print*, RES(ii,:)
     end if
     b = b + abs(RES(ii,3))
  end do
  
  print*, "max residual", a, "sum of residuals", b
  
  
  do ii=1,N
!!$     if (RES(ii,1).NE.RES(ii,1)) then
!!$        print*, ii, ROOTS(ii), RES(ii,1)
!!$        print*, RES(ii,:)
!!$     end if
     print*, ii, ROOTS(ii), RES(ii,1)
!!$     !print*, RES(ii,:)
  end do
  
  ! stop timer
  call system_clock(count=c_stop)
  print*, "Test took ", dble(c_stop-c_start)/dble(c_rate), "s"

  !print*, "MATLAB"
  !print*, "ROOTS=["
  !do ii=1,N
  !   print*, dble(ROOTS(ii)), aimag(ROOTS(ii))
  !end do
  !print*, "];"

  
end program example_d_polyc_roots
