#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_z_uprkfpen_scaled
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine chooses d k x k random matrices A_0, ..., A_d and computes 
! the eigenvalues of A_d \lambda^d + A_d-1 \lambda^d-1 + ... + A_1 \lambda + A_0
! using the QZ part of z_uprkfact_twistedqz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_z_uprkfpen_scaled


  implicit none
  
  ! compute variables
  integer, parameter :: dd = 10
  integer, parameter :: k = 10
  real(8), parameter :: norm = 1
  logical, parameter :: output=.FALSE.
  !logical, parameter :: output=.TRUE.
  integer :: N = dd*k
  integer :: ii, jj, ll, INFO, lwork, it
  complex(8), allocatable :: MA(:,:),MB(:,:), EIGS(:), REIGS(:), EIGSA(:), EIGSB(:)
  complex(8), allocatable :: REIGS2(:), V(:,:),W(:,:), WORK(:)
  complex(8), allocatable :: CDA(:,:), CDB(:,:), MC(:,:), MD(:,:)
  complex(8), allocatable :: REIGSA(:), REIGSB(:), VL(:,:),VR(:,:)
  complex(8), allocatable :: TA(:,:), TB(:,:), TL(:,:), Vt(:,:), Wt(:,:)

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
    function l_upr1fact_hess(m,flags)
      logical :: l_upr1fact_hess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_hess
  end interface

  ! timing variables
  integer:: c_start3, c_stop3, c_start2, c_stop2, c_start, c_stop, c_rate
  ! start timer
  call system_clock(count_rate=c_rate)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print banner
  print*,""
  print*,"example_z_uprkfact_scaledrandompencil:"
  print*, "K", k, "degree", dd 
  print*, "norm", norm
  print*,""
  
  allocate(MA(N,k),MB(N,k),rev(N),iev(N),ITS(N))

  allocate(Q(3*k*(N-1)),D1(2*k*(N+1)),C1(3*k*N),B1(3*k*N),P(N-2))
  allocate(D2(2*k*(N+1)),C2(3*k*N),B2(3*k*N),V(N,N),W(N,N),EIGS(N),REIGS(N))
  allocate(WORK(25*N+25),RWORK(8*N),REIGS2(N))
  allocate(REIGSA(N),REIGSB(N),VL(N,N),VR(N,N),EIGSA(N*k),EIGSB(N*k))
  allocate(CDA(N,N), CDB(N,N), MC(N,k), MD(N,k))
  allocate(TA(N,N), TB(N,N), TL(N,N), Vt(N,N), Wt(N,N))

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

  !MB(N,k)=cmplx(0d0,0d0,kind=8)
  

  ! scale columns to have norm norm
  do ii=1,k
     h = 0d0
     do ll=1,k
        h = h + abs(MA(ll,ii))**2
     end do

     do ll=k+1,dd*k
        h = h + abs(MA(ll,ii)+MB(ll-k,ii))**2
     end do

     do ll=dd*k+1,(dd+1)*k
        h = h + abs(MB(ll-k,ii))**2
     end do
     do ll=1,N
        MA(ll,ii) = MA(ll,ii)/h*norm
        MB(ll,ii) = MB(ll,ii)/h*norm
     end do
  end do
           

  MC = MA
  MD = MB
  
  
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

  CDA = V
  CDB = W

!!$  lwork = 25*N+25
!!$  call zggev('N','N', N, V, N, W, N, REIGSA, REIGSB, &
!!$       &VL, N, VR, N, WORK, lwork, RWORK, INFO) ! LAPACK
!!$  
!!$  do jj= 1,N
!!$     REIGS(jj) = REIGSA(jj)/REIGSB(jj)
!!$  end do
!!$
!!$  if (output) then
!!$     do jj= 1,N
!!$        print*, jj, REIGS(jj),REIGSA(jj),REIGSB(jj)
!!$     end do
!!$  end if
!!$  call system_clock(count=c_stop2)
!!$  print*, "Runtime unstructured QR solver (LAPACK) ",(dble(c_stop2-c_start2)/dble(c_rate))
!!$
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! Start Times
!!$  call system_clock(count=c_start)
!!$
!!$  call z_uprkdense_qz(.TRUE.,.FALSE.,N,k,MA,MB,EIGSA,EIGSB,V,W,INFO)
!!$
!!$  do ii=1,N
!!$     EIGS(ii) = EIGSA(ii)/EIGSB(ii)     
!!$  end do
!!$
!!$  call system_clock(count=c_stop)
!!$
!!$  maxerr = 0d0
!!$  do ii = 1,N
!!$    jj = 1
!!$    h = abs(EIGS(ii)-REIGS(jj))
!!$    do ll = 2,N
!!$       if (h>abs(EIGS(ii)-REIGS(ll))) then
!!$          jj = ll
!!$          h = abs(EIGS(ii)-REIGS(jj))
!!$       end if
!!$    end do
!!$    if (N.LT.100) then
!!$       print*, ii, EIGS(ii), REIGS(jj), abs(EIGS(ii)-REIGS(jj))
!!$    end if
!!$    if (abs(EIGS(ii)-REIGS(jj)).GT.maxerr) then
!!$       maxerr = abs(EIGS(ii)-REIGS(jj))
!!$    end if
!!$  end do
!!$
!!$  print*, "Maxmimal error vs. LAPACK", maxerr
!!$  print*, "Runtime   structured QR solver (eiscor) ",(dble(c_stop-c_start)/dble(c_rate))
!!$  if (N.LT.100) then
!!$     print*, "Runtime unstructured QR solver (LAPACK) ",(dble(c_stop2-c_start2)/dble(c_rate))
!!$  end if
!!$
!!$
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! Start Times
!!$  ! slow variant
!!$  MA = MC
!!$  MB = MD
!!$  call system_clock(count=c_start)
!!$
!!$  call z_uprkdense_qz_slow(.TRUE.,.FALSE.,N,k,MA,MB,EIGSA,EIGSB,V,W,INFO)
!!$
!!$  do ii=1,N
!!$     EIGS(ii) = EIGSA(ii)/EIGSB(ii)     
!!$  end do
!!$
!!$  call system_clock(count=c_stop)
!!$
!!$  maxerr = 0d0
!!$  do ii = 1,N
!!$    jj = 1
!!$    h = abs(EIGS(ii)-REIGS(jj))
!!$    do ll = 2,N
!!$       if (h>abs(EIGS(ii)-REIGS(ll))) then
!!$          jj = ll
!!$          h = abs(EIGS(ii)-REIGS(jj))
!!$       end if
!!$    end do
!!$    if (N.LT.100) then
!!$       print*, ii, EIGS(ii), REIGS(jj), abs(EIGS(ii)-REIGS(jj))
!!$    end if
!!$    if (abs(EIGS(ii)-REIGS(jj)).GT.maxerr) then
!!$       maxerr = abs(EIGS(ii)-REIGS(jj))
!!$    end if
!!$  end do
!!$
!!$  print*, "Maxmimal error vs. LAPACK", maxerr
!!$  print*, "Runtime   structured QR solver (eiscor_slow) ",(dble(c_stop-c_start)/dble(c_rate))
!!$  if (N.LT.100) then
!!$     print*, "Runtime unstructured QR solver (LAPACK) ",(dble(c_stop2-c_start2)/dble(c_rate))
!!$  end if
!!$
!!$
!!$  MA = MC
!!$  MB = MD
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Schur decomposition
  call system_clock(count=c_start3)

  call z_uprk_compress(.TRUE.,.TRUE.,.TRUE.,N,k,MA,MB,P,Q,&
       &D1,C1,B1,D2,C2,B2,V,W,INFO)
  if (INFO.NE.0) then
     print*, "Info code from z_uprkdense_factor: ", INFO
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! check factorization
  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D1,C1,B1,TA)
  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D2,C2,B2,TB)

!!$  call z_upr1fact_extracttri(.FALSE.,.TRUE.,N,&
!!$       &D1(1:(2*(N+1))),C1(1:(3*N)),&
!!$       &B1(1:(3*N)),TA)
!!$  do ll = 2,k
!!$     call z_upr1fact_extracttri(.FALSE.,.FALSE.,N,&
!!$          &D1(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C1(((ll-1)*3*N+1):(ll*3*N)),&
!!$          &B1(((ll-1)*3*N+1):(ll*3*N)),TA)
!!$  end do
!!$  do ll = 1,k
!!$     call z_upr1fact_extracttri(.FALSE.,(ll.EQ.1),N,&
!!$          &D2(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C2(((ll-1)*3*N+1):(ll*3*N)),&
!!$          &B2(((ll-1)*3*N+1):(ll*3*N)),TB)
!!$  end do

  ! Apply Q
  do jj = N-1,1,-1
     ! Apply Q(jj)
     Gf(1,1) = cmplx(Q(3*jj-2),Q(3*jj-1),kind=8)
     Gf(2,1) = cmplx(Q(3*jj),0d0,kind=8)
     Gf(1,2) = -Gf(2,1)
     Gf(2,2) = conjg(Gf(1,1))
     
     TA((jj):(jj+1),:) = matmul(Gf,TA((jj):(jj+1),:))
  end do

  do ii = 1,N
     do jj = 1,N
        Vt(ii,jj) = conjg(V(jj,ii))
     end do
  end do
  do ii = 1,N
     do jj = 1,N
        Wt(ii,jj) = conjg(W(jj,ii))
     end do
  end do

  TA = matmul(W,matmul(TA,Vt))

  !print*, TA
  !print*, CDA
  
  h = 0d0
  do jj=1,N
     do ii=1,N
        if (output.AND.abs(TA(ii,jj)-CDA(ii,jj)).GT.1d-13) then
           print*, ii, jj, TA(ii,jj), CDA(ii,jj), abs(TA(ii,jj)-CDA(ii,jj))
        end if
        h = h + abs(TA(ii,jj)-CDA(ii,jj))**2
     end do
  end do
  print*, "Backward Error (factorization, TA): ", sqrt(h)

  TB = matmul(W,matmul(TB,Vt))


  h = 0d0
  do jj=1,N
     do ii=1,N
        if (output.AND.abs(TB(ii,jj)-CDB(ii,jj)).GT.1d-13) then
           print*, ii, jj, TB(ii,jj), CDB(ii,jj), abs(TB(ii,jj)-CDB(ii,jj))
        end if
        h = h + abs(TB(ii,jj)-CDB(ii,jj))**2
     end do
  end do
  print*, "Backward Error (factorization, TB): ", sqrt(h)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! (twisted) Hessenberg QZ
  call z_uprkfpen_qz(.TRUE.,.FALSE.,l_upr1fact_hess,N,k,&
       &P,Q,D1,C1,B1,D2,C2,B2,N,V,W,ITS,INFO)
  if (INFO.NE.0) then
     print*, "Info code from z_uprkfact_twistedqz: ", INFO
  end if

  it = 0
  do ll = 1,N
     it = it + ITS(ll)
     !print*, ITS(ll)
  end do
  print*, "Iterations per eigenvalue: ", (1d0*it)/N

  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D1,C1,B1,TA)
  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D2,C2,B2,TB)
!!$  do ll = 1,k
!!$     call z_upr1fact_extracttri(.FALSE.,(ll.EQ.1),N,&
!!$          &D1(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C1(((ll-1)*3*N+1):(ll*3*N)),&
!!$          &B1(((ll-1)*3*N+1):(ll*3*N)),TA)
!!$  end do
!!$  do ll = 1,k
!!$     call z_upr1fact_extracttri(.FALSE.,(ll.EQ.1),N,&
!!$          &D2(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C2(((ll-1)*3*N+1):(ll*3*N)),&
!!$          &B2(((ll-1)*3*N+1):(ll*3*N)),TB)
!!$  end do
 
  do ii = 1,N
     do jj = 1,N
        Vt(ii,jj) = conjg(V(jj,ii))
     end do
  end do
  do ii = 1,N
     do jj = 1,N
        Wt(ii,jj) = conjg(W(jj,ii))
     end do
  end do
  call system_clock(count=c_stop3)

  TA = matmul(W,matmul(TA,Vt))
  h = 0d0
  do jj=1,N
     do ii=1,N
        if (output.AND.abs(TA(ii,jj)-CDA(ii,jj)).GT.1d-13) then
           print*, ii, jj, TA(ii,jj), CDA(ii,jj), abs(TA(ii,jj)-CDA(ii,jj))
        end if
        h = h + abs(TA(ii,jj)-CDA(ii,jj))**2
     end do
  end do
  print*, "Backward Error (Schur decomposition, TA): ", sqrt(h)

  TB = matmul(W,matmul(TB,Vt))
  h = 0d0
  do jj=1,N
     do ii=1,N
        if (output.AND.abs(TB(ii,jj)-CDB(ii,jj)).GT.1d-13) then
           print*, ii, jj, TB(ii,jj), CDB(ii,jj), abs(TB(ii,jj)-CDB(ii,jj))
        end if
        h = h + abs(TB(ii,jj)-CDB(ii,jj))**2
     end do
  end do
  print*, "Backward Error (Schur decomposition, TB): ", sqrt(h)
  print*, "Runtime structured QR solver (eiscor) with eigenvectors ",(dble(c_stop3-c_start3)/dble(c_rate))

  call zggev('N','N', N, TA, N, TB, N, REIGSA, REIGSB, &
       &VL, N, VR, N, WORK, lwork, RWORK, INFO)  

  do jj= 1,N
     REIGS(jj) = REIGSA(jj)/REIGSB(jj)
  end do


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
  deallocate(CDA,CDB,MC,MD)
  deallocate(TA,TB,TL,Vt,Wt)

end program example_z_uprkfpen_scaled
