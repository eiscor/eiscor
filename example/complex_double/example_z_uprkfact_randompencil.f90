#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_z_uprkfact_randommatr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine chooses d k x k random matrices A_0, ..., A_k-1 and computes 
! the eigenvalues of I \lambda^d + A_d-1 \lambda^d-1 + ... + A_1 \lambda + A_0
! using z_uprkfact_twistedqz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_z_uprkfact_randommatr


  implicit none
  
  ! compute variables
  integer, parameter :: dd = 50
  integer, parameter :: k = 5
  logical, parameter :: output=.FALSE.
  !logical, parameter :: output=.TRUE.
  integer :: N = dd*k
  integer :: ii, jj, ll, INFO, lwork, it
  complex(8), allocatable :: MA(:,:),MB(:,:), EIGS(:), REIGS(:), EIGSA(:), EIGSB(:)
  complex(8), allocatable :: REIGS2(:), V(:,:),W(:,:), WORK(:)
  complex(8), allocatable :: REIGSA(:), REIGSB(:), VL(:,:),VR(:,:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:), D1(:), C1(:), B1(:)
  real(8), allocatable :: D2(:), C2(:), B2(:), RWORK(:)
  logical, allocatable :: P(:)
  real(8) :: h, maxerr
  complex(8) :: Gf(2,2)

  ! real and imag part of eigenvalues
  double precision, allocatable :: rev(:), iev(:)
  ! interface
  interface
    function l_upr1fact_upperhess(m,flags)
      logical :: l_upr1fact_upperhess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_upperhess
  end interface

  ! timing variables
  integer:: c_start2, c_stop2, c_start, c_stop, c_rate
  ! start timer
  call system_clock(count_rate=c_rate)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print banner
  print*,""
  print*,"example_z_uprkfact_randompencil:"
  print*,""
  
  allocate(MA(N,k),MB(N,k),rev(N),iev(N),ITS(N))

  allocate(Q(3*k*(N-1)),D1(2*k*(N+1)),C1(3*k*N),B1(3*k*N),P(N-2))
  allocate(D2(2*k*(N+1)),C2(3*k*N),B2(3*k*N),V(N,N),W(N,N),EIGS(N),REIGS(N))
  allocate(WORK(25*N+25),RWORK(8*N),REIGS2(N))
  allocate(REIGSA(N),REIGSB(N),VL(N,N),VR(N,N),EIGSA(N*k),EIGSB(N*k))

  call u_fixedseed_initialize(INFO)  
  call z_2Darray_random_normal(N,k,MA)
  call z_2Darray_random_normal(N,k,MB)


  ! make P0 lower triangular
  do ii=1,k-1
     do jj=ii+1,k
        MA(ii,jj)=cmplx(0d0,0d0,kind=8)
     end do
  end do

  ! make Pd upper triangular
  do ii=2,k
     do jj=1,ii-1
        MB(N-k+ii,jj)=cmplx(0d0,0d0,kind=8)
     end do
  end do

  !print*, MB

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute roots with LAPACK
  call system_clock(count=c_start2)
  V = cmplx(0d0,0d0,kind=8)
  W = cmplx(0d0,0d0,kind=8)

  do ii=1,N-k
     V(ii+k,ii)=cmplx(1d0,0d0,kind=8)
     W(ii,ii)=cmplx(1d0,0d0,kind=8)
  end do
  V(:,N-k+1:N) = MA
  W(:,N-k+1:N) = MB


  if (DEBUG) then
     ! Check V and W
     print*, "V"
     do ii=1,N
        write (*,"(ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3)") V(ii,:)
     end do
     
     print*, "W"
     do ii=1,N
        write (*,"(ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3)") W(ii,:)
     end do
  end if

  lwork = 25*N+25
  call zggev('N','N', N, V, N, W, N, REIGSA, REIGSB, &
       &VL, N, VR, N, WORK, lwork, RWORK, INFO)
  
  do jj= 1,N
     REIGS(jj) = REIGSA(jj)/REIGSB(jj)
  end do

  if (output) then
     do jj= 1,N
        print*, jj, REIGS(jj),REIGSA(jj),REIGSB(jj)
     end do
  end if
  call system_clock(count=c_stop2)
  print*, "Runtime unstructured QR solver (LAPACK) ",(dble(c_stop2-c_start2)/dble(c_rate))

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Start Times
  call system_clock(count=c_start)

  call z_uprkdense_qz(.TRUE.,.FALSE.,N,k,MA,MB,EIGSA,EIGSB,V,W,INFO)

!!$  call z_uprkdense_factor(.TRUE.,.FALSE.,.FALSE.,N,k,MA,MB,P,Q,&
!!$       &D1,C1,B1,D2,C2,B2,V,W,INFO)
!!$  if (INFO.NE.0) then
!!$     print*, "Info code from z_uprkdense_factor: ", INFO
!!$  end if
!!$
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! check factorization
!!$  do ll = 1,k
!!$     call z_upr1fact_extracttri(.FALSE.,N,&
!!$          &D1(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C1(((ll-1)*3*N+1):(ll*3*N)),&
!!$          &B1(((ll-1)*3*N+1):(ll*3*N)),VL)
!!$     if (ll.EQ.1) then
!!$        V = VL
!!$     else
!!$        V = matmul(V,VL)
!!$     end if     
!!$  end do
!!$  ! Apply Q
!!$  do jj = N-1,1,-1
!!$     ! Apply Q(jj)
!!$     Gf(1,1) = cmplx(Q(3*jj-2),Q(3*jj-1),kind=8)
!!$     Gf(2,1) = cmplx(Q(3*jj),0d0,kind=8)
!!$     Gf(1,2) = -Gf(2,1)
!!$     Gf(2,2) = conjg(Gf(1,1))
!!$     
!!$     V((jj):(jj+1),:) = matmul(Gf,V((jj):(jj+1),:))
!!$  end do
!!$
!!$  do ll = 1,k
!!$     call z_upr1fact_extracttri(.FALSE.,N,&
!!$          &D2(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C2(((ll-1)*3*N+1):(ll*3*N)),&
!!$          &B2(((ll-1)*3*N+1):(ll*3*N)),VR)
!!$     if (ll.EQ.1) then
!!$        W = VR
!!$     else
!!$        W = matmul(W,VR)
!!$     end if     
!!$  end do
!!$
!!$  if (DEBUG) then
!!$     ! Check V and W
!!$     print*, "V"
!!$     do ii=1,N
!!$        write (*,"(ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3)") V(ii,:)
!!$     end do
!!$     
!!$     print*, "W"
!!$     do ii=1,N
!!$        write (*,"(ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3,ES13.4E3)") W(ii,:)
!!$     end do
!!$
!!$     call zggev('N','N', N, V, N, W, N, REIGSA, REIGSB, &
!!$          &VL, N, VR, N, WORK, lwork, RWORK, INFO)     
!!$     
!!$     do jj= 1,N
!!$        REIGS(jj) = REIGSA(jj)/REIGSB(jj)
!!$        print*, jj, REIGS(jj),REIGSA(jj),REIGSB(jj)
!!$     end do
!!$  end if
!!$
!!$  
!!$
!!$
!!$  call z_uprkfact_twistedqz(.TRUE.,.FALSE.,.FALSE.,l_upr1fact_upperhess,N,k,&
!!$       &P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
!!$  if (INFO.NE.0) then
!!$     print*, "Info code from z_uprkfact_twistedqz: ", INFO
!!$  end if
!!$
!!$  it = 0
!!$  do ll = 1,N
!!$     it = it + ITS(ll)
!!$     !print*, ITS(ll)
!!$  end do
!!$  print*, "Iterations per eigenvalue: ", (1d0*it)/N
!!$
!!$  do ll = 1,k
!!$     call z_upr1fact_extracttri(.TRUE.,N,&
!!$          &D1(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C1(((ll-1)*3*N+1):(ll*3*N)),&
!!$          &B1(((ll-1)*3*N+1):(ll*3*N)),V)
!!$     EIGSA(((ll-1)*N+1):ll*N) = V(1:N,1)
!!$  end do
!!$  do ll = 2,k
!!$     do ii=1,N
!!$        EIGSA(ii) = EIGSA(ii)*EIGSA((ll-1)*N+ii)     
!!$     end do
!!$  end do
!!$
!!$  do ll = 1,k
!!$     call z_upr1fact_extracttri(.TRUE.,N,&
!!$          &D2(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C2(((ll-1)*3*N+1):(ll*3*N)),&
!!$          &B2(((ll-1)*3*N+1):(ll*3*N)),V)
!!$     EIGSB(((ll-1)*N+1):ll*N) = V(1:N,1)
!!$  end do
!!$  do ll = 2,k
!!$     do ii=1,N
!!$        EIGSB(ii) = EIGSB(ii)*EIGSB((ll-1)*N+ii)     
!!$     end do
!!$  end do

  do ii=1,N
     EIGS(ii) = EIGSA(ii)/EIGSB(ii)     
  end do

  call system_clock(count=c_stop)


  maxerr = 0d0
  do ii = 1,N
    jj = 1
    h = abs(EIGS(ii)-REIGS(jj))
    do ll = 2,N
       if (h>abs(EIGS(ii)-REIGS(ll))) then
          jj = ll
          h = abs(EIGS(ii)-REIGS(jj))
       end if
    end do
    if (N.LT.100) then
       print*, ii, EIGS(ii), REIGS(jj), abs(EIGS(ii)-REIGS(jj))
    end if
    if (abs(EIGS(ii)-REIGS(jj)).GT.maxerr) then
       maxerr = abs(EIGS(ii)-REIGS(jj))
    end if
  end do

  print*, "Maxmimal error vs. LAPACK", maxerr
  print*, "Runtime   structured QR solver (eiscor) ",(dble(c_stop-c_start)/dble(c_rate))
  if (N.LT.100) then
     print*, "Runtime unstructured QR solver (LAPACK) ",(dble(c_stop2-c_start2)/dble(c_rate))
  end if
  

  if (output) then
     print*, "The Eigenvalues"
     do ii = 1,N
        print*, ii,  REIGS(ii)
     end do
  end if


  print*,""
  deallocate(MA,MB,rev,iev,ITS)
  deallocate(Q,D1,C1,B1,P)
  deallocate(D2,C2,B2,V,W,EIGS,REIGS,REIGS2)
  deallocate(WORK,RWORK)
  deallocate(REIGSA,REIGSB,VL,VR,EIGSA,EIGSB)

end program example_z_uprkfact_randommatr
