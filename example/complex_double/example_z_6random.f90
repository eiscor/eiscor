#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_z_uprkfpen_2random
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine chooses d k x k random matrices A_0, ..., A_d and computes 
! the eigenvalues of A_d \lambda^d + A_d-1 \lambda^d-1 + ... + A_1 \lambda + A_0
! using z_uprkdense_qz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_z_uprkfpen_6random


  implicit none
  
  ! compute variables
  integer, parameter :: dd = 20
  integer, parameter :: k = 8
  logical, parameter :: output=.FALSE.
  !logical, parameter :: output=.TRUE.
  integer :: N = dd*k
  integer :: ii, jj, ll, INFO, lwork, it, qq
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
  print*,"example_z_uprkfact_2randompencil:"
  print*, "K", k, "degree", dd, "N", N
  print*,""
  
  allocate(MA(N,k),MB(N,k),rev(N),iev(N),ITS(N))

  allocate(Q(3*k*(N-1)),D1(2*k*N),C1(3*k*N),B1(3*k*N),P(N-2))
  allocate(D2(2*k*N),C2(3*k*N),B2(3*k*N),V(N,N),W(N,N),EIGS(N),REIGS(N))
  allocate(WORK(25*N+25),RWORK(8*N),REIGS2(N))
  allocate(REIGSA(N),REIGSB(N),VL(N,N),VR(N,N),EIGSA(N),EIGSB(N))
  allocate(CDA(N,N), CDB(N,N), MC(N,k), MD(N,k))
  allocate(TA(N,N), TB(N,N), TL(N,N), Vt(N,N), Wt(N,N))
  call u_fixedseed_initialize(INFO)  
  call z_2Darray_random_normal(N,k,MA)
  call z_2Darray_random_normal(N,k,MB)


  ! make P0 upper triangular
  do ii=2,k
     do jj=1,ii-1
        MA(ii,jj)=cmplx(0d0,0d0,kind=8)
     end do
  end do

  ! make Pd upper triangular
  do ii=2,k
     do jj=1,ii-1
        MB(N-k+ii,jj)=cmplx(0d0,0d0,kind=8)
     end do
  end do
  
  MC = MA
  MD = MB

  do qq=1,2
  
     if (qq.EQ.2) then
        
        MA = MC
        MB = MD
        
        !MB(N,k)=cmplx(0d0,0d0,kind=8)
        MA(1,1)=cmplx(0d0,0d0,kind=8)
        
     end if
  
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

  lwork = 25*N+25
  call zggev('N','N', N, V, N, W, N, REIGSA, REIGSB, &
       &VL, N, VR, N, WORK, lwork, RWORK, INFO) ! LAPACK
  
  do jj= 1,N
     REIGS(jj) = REIGSA(jj)/REIGSB(jj)
     !print*, "LAPACK", jj, REIGS(jj), REIGSA(jj), REIGSB(jj)
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

  call z_uprkdense_qz2(.FALSE.,.FALSE.,N,K,MA,MB,N,EIGSA,EIGSB,V,W,TA,TB,INFO)
  !print*, "         No EIGSA                                                 EIGSB"

  do ii=1,N
     !print*, ii, EIGSA(ii), EIGSB(ii)
     EIGS(ii) = EIGSA(ii)/EIGSB(ii)
     !print*, ii, EIGS(ii), EIGSA(ii), EIGSB(ii)
  end do

  call system_clock(count=c_stop)

  maxerr = 0d0
  
  do ii = 1,N
     if (abs(EIGSB(ii)).LE.EISCOR_DBL_EPS*abs(EIGSA(ii))) then
        jj = 1
        h = abs(REIGSB(jj))
        do ll = 2,N
           if (h>abs(REIGSB(ll))) then
              jj = ll
              h = abs(REIGSB(jj))
           end if
        end do
        if (abs(REIGSB(ii)).LE.EISCOR_DBL_EPS*abs(REIGSA(ii))) then
           maxerr = EISCOR_DBL_INF
        end if
          
     else  
        jj = 1
        do while ((jj.LE.N).AND.(REIGS(jj).NE.REIGS(jj)))
           jj = jj+1
        end do
        h = abs(EIGS(ii)-REIGS(jj))
        ll = 2
        do while (ll.LE.N)
           do while ((ll.LE.N).AND.(REIGS(ll).NE.REIGS(ll)))
              ll = ll+1
           end do
           if (h>abs(EIGS(ii)-REIGS(ll))) then
              jj = ll
              h = abs(EIGS(ii)-REIGS(jj))
           end if
           ll = ll+1
        end do
        if (N.LT.100) then
           print*, ii, EIGS(ii), REIGS(jj), abs(EIGS(ii)-REIGS(jj))/abs(EIGS(jj))
        end if
        !if (abs(EIGS(ii)).EQ.abs(EIGS(ii))) then
        !   if (abs(REIGS(jj)).EQ.abs(REIGS(jj))) then
        if (abs(EIGS(ii)-REIGS(jj)).GT.maxerr) then
           maxerr = abs(EIGS(ii)-REIGS(jj))
        end if
        !   end if
        !end if
     end if
     
  end do
  
 
  
!!$  maxerr = 0d0
!!$  !print*, "         No EIGS                                                  ",&
!!$  !       &"REIGS                                                  DIFF"
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
!!$       print*, ii, EIGS(ii), REIGS(jj), abs(EIGS(ii)-REIGS(jj))/abs(EIGS(jj))
!!$    end if
!!$    if (abs(EIGS(ii)).EQ.abs(EIGS(ii))) then
!!$       if (abs(REIGS(jj)).EQ.abs(REIGS(jj))) then
!!$          if (abs(EIGS(ii)-REIGS(jj)).GT.maxerr) then
!!$             maxerr = abs(EIGS(ii)-REIGS(jj))
!!$          end if
!!$       end if
!!$    end if
!!$  end do
!!$
!!$  print*, "Maximal error vs. LAPACK", maxerr
!!$
!!$  maxerr = 0d0
!!$  do ii = 1,N
!!$    jj = 1
!!$    h = abs(REIGS(ii)-EIGS(jj))
!!$    do ll = 2,N
!!$       if (h>abs(REIGS(ii)-EIGS(ll))) then
!!$          jj = ll
!!$          h = abs(REIGS(ii)-EIGS(jj))
!!$       end if
!!$    end do
!!$    if (N.LT.100) then
!!$       print*, ii, EIGS(ii), REIGS(jj), abs(EIGS(ii)-REIGS(jj))/abs(EIGS(jj))
!!$    end if
!!$    if (abs(EIGS(jj)).EQ.abs(EIGS(jj))) then
!!$       if (abs(REIGS(ii)).EQ.abs(REIGS(ii))) then
!!$          if (abs(REIGS(ii)-EIGS(jj)).GT.maxerr) then
!!$             maxerr = abs(REIGS(ii)-EIGS(jj))
!!$          end if
!!$       end if
!!$    end if
!!$  end do

  print*, "Maximal error vs. LAPACK", maxerr
  print*, "Runtime   structured QR solver (eiscor) ",(dble(c_stop-c_start)/dble(c_rate))
  if (N.LT.100) then
     print*, "Runtime unstructured QR solver (LAPACK) ",(dble(c_stop2-c_start2)/dble(c_rate))
  end if


!  do ii = 1,N
!     print*, "LAPACK",ii, REIGS(ii), REIGSA(ii), REIGSB(ii)
!  end do
!
!    do ii = 1,N
!     print*, "EISCOR",ii, EIGS(ii), EIGSA(ii), EIGSB(ii)
!  end do

!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! Start Times
!!$  ! slow variant
!!$  MA = MC
!!$  MB = MD
!!$  call system_clock(count=c_start)
!!$
!!$  call z_uprkdense_qz2_slow(.FALSE.,N,k,MA,MB,N,EIGSA,EIGSB,V,W,TA,TB,INFO)
!!$
!!$  do ii=1,N
!!$     EIGS(ii) = EIGSA(ii)/EIGSB(ii)     
!!$     if (abs(EIGSB(ii))<1e-20) then
!!$        EIGS(ii) = cmplx(EISCOR_DBL_INF,EISCOR_DBL_INF,kind=8)
!!$     end if
!!$  end do
!!$
!!$  call system_clock(count=c_stop)
!!$
!!$  maxerr = 0d0
!!$  
!!$  do ii = 1,N
!!$     if (abs(EIGSB(ii)).LE.EISCOR_DBL_EPS*abs(EIGSA(ii))) then
!!$        jj = 1
!!$        h = abs(REIGSB(jj))
!!$        do ll = 2,N
!!$           if (h>abs(REIGSB(ll))) then
!!$              jj = ll
!!$              h = abs(REIGSB(jj))
!!$           end if
!!$        end do
!!$        if (abs(REIGSB(ii)).LE.EISCOR_DBL_EPS*abs(REIGSA(ii))) then
!!$           maxerr = EISCOR_DBL_INF
!!$        end if
!!$          
!!$     else  
!!$        jj = 1
!!$        do while ((jj.LE.N).AND.(REIGS(jj).NE.REIGS(jj)))
!!$           jj = jj+1
!!$        end do
!!$        h = abs(EIGS(ii)-REIGS(jj))
!!$        ll = 2
!!$        do while (ll.LE.N)
!!$           do while ((ll.LE.N).AND.(REIGS(ll).NE.REIGS(ll)))
!!$              ll = ll+1
!!$           end do
!!$           if (h>abs(EIGS(ii)-REIGS(ll))) then
!!$              jj = ll
!!$              h = abs(EIGS(ii)-REIGS(jj))
!!$           end if
!!$           ll = ll+1
!!$        end do
!!$        if (N.LT.100) then
!!$           print*, ii, EIGS(ii), REIGS(jj), abs(EIGS(ii)-REIGS(jj))/abs(EIGS(jj))
!!$        end if
!!$        !if (abs(EIGS(ii)).EQ.abs(EIGS(ii))) then
!!$        !   if (abs(REIGS(jj)).EQ.abs(REIGS(jj))) then
!!$        if (abs(EIGS(ii)-REIGS(jj)).GT.maxerr) then
!!$           maxerr = abs(EIGS(ii)-REIGS(jj))
!!$        end if
!!$        !   end if
!!$        !end if
!!$     end if
!!$     
!!$  end do
!!$  print*, "Maximal error vs. LAPACK", maxerr
!!$
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

  print*, "Maximal error vs. LAPACK", maxerr

  maxerr = 0d0
  do ii = 1,N
    jj = 1
    h = abs(REIGS(ii)-EIGS(jj))
    do ll = 2,N
       if (h>abs(REIGS(ii)-EIGS(ll))) then
          jj = ll
          h = abs(REIGS(ii)-EIGS(jj))
       end if
    end do
    if (N.LT.100) then
       print*, ii, EIGS(ii), REIGS(jj), abs(EIGS(ii)-REIGS(jj))/abs(EIGS(jj))
    end if
    if (abs(EIGS(jj)).EQ.abs(EIGS(jj))) then
       if (abs(REIGS(ii)).EQ.abs(REIGS(ii))) then
          if (abs(REIGS(ii)-EIGS(jj)).GT.maxerr) then
             maxerr = abs(REIGS(ii)-EIGS(jj))
          end if
       end if
    end if
  end do

  print*, "Maximal error vs. LAPACK", maxerr
!!$
!!$  print*, "Runtime   structured QR solver (eiscor_slow) ",(dble(c_stop-c_start)/dble(c_rate))
!!$  if (N.LT.100) then
!!$     print*, "Runtime unstructured QR solver (LAPACK) ",(dble(c_stop2-c_start2)/dble(c_rate))
!!$  end if
!!$
!!$
!!$  MA = MC
!!$  MB = MD
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! Schur decomposition
!!$  call system_clock(count=c_start3)
!!$
!!$  call z_uprk_compress2(.TRUE.,.TRUE.,.TRUE.,N,K,MA,MB,P,Q,&
!!$       &D1,C1,B1,D2,C2,B2,V,W,INFO)
!!$  if (INFO.NE.0) then
!!$     print*, "Info code from z_uprkdense_factor: ", INFO
!!$  end if
!!$
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! check factorization
!!$  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D1,C1,B1,TA)
!!$  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D2,C2,B2,TB)
  do ll = 1,k
     ! z_uprkutri_decompress(DIAG,N,K,STR,STP,D,C,B,T)
     call z_upr1utri_decompress(.TRUE.,N,&
          &D1(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C1(((ll-1)*3*N+1):(ll*3*N)),&
          &B1(((ll-1)*3*N+1):(ll*3*N)),TL)
     if (ll.EQ.1) then
        TA = TL
     else
        TA = matmul(TA,TL)
     end if     
  end do
  do ll = 1,k
     call z_upr1utri_decompress(.TRUE.,N,&
          &D2(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C2(((ll-1)*3*N+1):(ll*3*N)),&
          &B2(((ll-1)*3*N+1):(ll*3*N)),TL)
     if (ll.EQ.1) then
        TB = TL
     else
        TB = matmul(TB,TL)
     end if     
  end do
!!$
!!$  ! Apply Q
!!$  do jj = N-1,1,-1
!!$     ! Apply Q(jj)
!!$     Gf(1,1) = cmplx(Q(3*jj-2),Q(3*jj-1),kind=8)
!!$     Gf(2,1) = cmplx(Q(3*jj),0d0,kind=8)
!!$     Gf(1,2) = -Gf(2,1)
!!$     Gf(2,2) = conjg(Gf(1,1))
!!$     
!!$     TA((jj):(jj+1),:) = matmul(Gf,TA((jj):(jj+1),:))
!!$  end do
!!$
!!$  do ii = 1,N
!!$     do jj = 1,N
!!$        Vt(ii,jj) = conjg(V(jj,ii))
!!$     end do
!!$  end do
!!$  do ii = 1,N
!!$     do jj = 1,N
!!$        Wt(ii,jj) = conjg(W(jj,ii))
!!$     end do
!!$  end do
!!$
!!$  TA = matmul(W,matmul(TA,Vt))
!!$
!!$  h = 0d0
!!$  do jj=1,N
!!$     do ii=1,N
!!$        if (output.AND.abs(TA(ii,jj)-CDA(ii,jj)).GT.1d-13) then
!!$           print*, ii, jj, TA(ii,jj), CDA(ii,jj), abs(TA(ii,jj)-CDA(ii,jj))
!!$        end if
!!$        h = h + abs(TA(ii,jj)-CDA(ii,jj))**2
!!$     end do
!!$  end do
!!$  print*, "Backward Error (factorization, TA): ", sqrt(h)
!!$
!!$  TB = matmul(W,matmul(TB,Vt))
!!$
!!$
!!$  h = 0d0
!!$  do jj=1,N
!!$     do ii=1,N
!!$        if (output.AND.abs(TB(ii,jj)-CDB(ii,jj)).GT.1d-13) then
!!$           print*, ii, jj, TB(ii,jj), CDB(ii,jj), abs(TB(ii,jj)-CDB(ii,jj))
!!$        end if
!!$        h = h + abs(TB(ii,jj)-CDB(ii,jj))**2
!!$     end do
!!$  end do
!!$  print*, "Backward Error (factorization, TB): ", sqrt(h)
!!$
!!$
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! (twisted) Hessenberg QZ
!!$  call z_uprkfpen_qz(.TRUE.,.FALSE.,l_upr1fact_hess,N,k,&
!!$       &P,Q,D1,C1,B1,D2,C2,B2,N,V,W,ITS,INFO)
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
!!$  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D1,C1,B1,TA)
!!$  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D2,C2,B2,TB)
  do ll = 1,k
     call z_upr1utri_decompress(.TRUE.,N,&
          &D1(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C1(((ll-1)*3*N+1):(ll*3*N)),&
          &B1(((ll-1)*3*N+1):(ll*3*N)),TL)
     if (ll.EQ.1) then
        TA = TL
     else
        TA = matmul(TA,TL)
     end if     
  end do
  do ll = 1,k
     call z_upr1utri_decompress(.TRUE.,N,&
          &D2(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C2(((ll-1)*3*N+1):(ll*3*N)),&
          &B2(((ll-1)*3*N+1):(ll*3*N)),TL)
     if (ll.EQ.1) then
        TB = TL
     else
        TB = matmul(TB,TL)
     end if     
  end do
!!$ 
!!$  do ii = 1,N
!!$     do jj = 1,N
!!$        Vt(ii,jj) = conjg(V(jj,ii))
!!$     end do
!!$  end do
!!$  do ii = 1,N
!!$     do jj = 1,N
!!$        Wt(ii,jj) = conjg(W(jj,ii))
!!$     end do
!!$  end do
!!$  call system_clock(count=c_stop3)
!!$
!!$  TA = matmul(W,matmul(TA,Vt))
!!$  h = 0d0
!!$  do jj=1,N
!!$     do ii=1,N
!!$        if (output.AND.abs(TA(ii,jj)-CDA(ii,jj)).GT.1d-13) then
!!$           print*, ii, jj, TA(ii,jj), CDA(ii,jj), abs(TA(ii,jj)-CDA(ii,jj))
!!$        end if
!!$        h = h + abs(TA(ii,jj)-CDA(ii,jj))**2
!!$     end do
!!$  end do
!!$  print*, "Backward Error (Schur decomposition, TA): ", sqrt(h)
!!$
!!$  TB = matmul(W,matmul(TB,Vt))
!!$  h = 0d0
!!$  do jj=1,N
!!$     do ii=1,N
!!$        if (output.AND.abs(TB(ii,jj)-CDB(ii,jj)).GT.1d-13) then
!!$           print*, ii, jj, TB(ii,jj), CDB(ii,jj), abs(TB(ii,jj)-CDB(ii,jj))
!!$        end if
!!$        h = h + abs(TB(ii,jj)-CDB(ii,jj))**2
!!$     end do
!!$  end do
!!$  print*, "Backward Error (Schur decomposition, TB): ", sqrt(h)
!!$  print*, "Runtime structured QR solver (eiscor) with eigenvectors ",(dble(c_stop3-c_start3)/dble(c_rate))

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

end do
  
  print*,""
  deallocate(MA,MB,rev,iev,ITS)
  deallocate(Q,D1,C1,B1,P)
  deallocate(D2,C2,B2,V,W,EIGS,REIGS,REIGS2)
  deallocate(WORK,RWORK)
  deallocate(REIGSA,REIGSB,VL,VR,EIGSA,EIGSB)
  deallocate(CDA,CDB,MC,MD)
  deallocate(TA,TB,TL,Vt,Wt)

end program example_z_uprkfpen_6random
