#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_z_uprkfpen_2manyscaled
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine chooses d k x k random matrices A_0, ..., A_d and computes 
! the eigenvalues of A_d \lambda^d + A_d-1 \lambda^d-1 + ... + A_1 \lambda + A_0
! using the QZ part of z_uprkfact_twistedqz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_z_uprkfpen_2manyscaled


  implicit none
  
  ! compute variables
  integer, parameter :: dd = 10
  integer, parameter :: kpara = 100
  real(8) :: norm, norma, expo, normb
  logical, parameter :: output=.FALSE., david = .TRUE.
  !logical, parameter :: output=.TRUE.
  integer :: N = dd*kpara, maxnumber, k, kin, Nin
  integer :: ii, jj, ll, co,co2, INFO, lwork, it, sumit
  complex(8), allocatable :: MA(:,:),MB(:,:), EIGS(:), REIGS(:), EIGSA(:), EIGSB(:)
  complex(8), allocatable :: REIGS2(:), V(:,:),W(:,:), WORK(:)
  complex(8), allocatable :: CDA(:,:), CDB(:,:), MC(:,:), MD(:,:)
  complex(8), allocatable :: REIGSA(:), REIGSB(:), VL(:,:),VR(:,:)
  complex(8), allocatable :: TA(:,:), TB(:,:), TL(:,:), Vt(:,:), Wt(:,:)

  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:), D1(:), C1(:), B1(:)
  real(8), allocatable :: D2(:), C2(:), B2(:), RWORK(:)
  logical, allocatable :: P(:)
  real(8) :: h, maxerr, ru,rv,rw,s,afactbe,bfactbe
  complex(8) :: Gf(2,2)

  ! real and imag part of eigenvalues
  double precision, allocatable :: rev(:), iev(:)
  character(len=32) :: arg
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

!  if (iargc()>0) then
!     call getarg(1, arg)
!     read (arg,'(I10)') maxnumber
!     if (iargc()>1) then
!        call getarg(2, arg)
!        read (arg,'(I10)') k
!        print*, k
!        if (k.GT.kpara) then
!           k = kpara
!        end if
!     else
!        k = kpara
!     end if
!  else
!     maxnumber = -1
!     k = kpara
!  end if

  maxnumber = 10000
  k = 4
  
  N = dd*k  
  allocate(MA(N,k),MB(N,k),rev(N),iev(N),ITS(N))

  allocate(Q(3*k*(N-1)),D1(2*k*(N+1)),C1(3*k*N),B1(3*k*N),P(N-2))
  allocate(D2(2*k*(N+1)),C2(3*k*N),B2(3*k*N),V(N,N),W(N,N),EIGS(N),REIGS(N))
  allocate(WORK(25*N+25),RWORK(8*N),REIGS2(N))
  allocate(REIGSA(N),REIGSB(N),VL(N,N),VR(N,N),EIGSA(N*k),EIGSB(N*k))
  allocate(CDA(N,N), CDB(N,N), MC(N,k), MD(N,k))
  allocate(TA(N,N), TB(N,N), TL(N,N), Vt(N,N), Wt(N,N))


  open (unit=7, file="err.txt", status='unknown')
  open (unit=8, file="errA.txt", status='unknown')
  open (unit=9, file="errB.txt", status='unknown')
  
  open (unit=21, file="MA.txt", status='unknown')
  open (unit=22, file="MB.txt", status='unknown')
  open (unit=23, file="nk.txt", status='unknown')
  
  read (23,"(I16)") kin
  read (23,"(I16)") Nin
  if (kin.NE.k) then
     k = kin
     if (Nin.NE.N) then
        N = k*dd
     end if
  end if
  
  
  
  if (david) then
     open (unit=10, file="errFactAdavid.dat", status='unknown')
     open (unit=11, file='errFactBdavid.dat', status='unknown')
     open (unit=12, file='errAdavid.dat', status = 'unknown')
     open (unit=13, file='errBdavid.dat', status = 'unknown')
  else
     open (unit=10, file="errFactA.dat", status='unknown')
     open (unit=11, file='errFactB.dat', status='unknown')
     open (unit=12, file='errA.dat', status = 'unknown')
     open (unit=13, file='errB.dat', status = 'unknown')
  end if
  
  !call u_fixedseed_initialize(INFO)  
  call u_randomseed_initialize(INFO)
  
  write (7,*) ""
  write (8,*) ""
  write (9,*) ""

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print banner
  print*,""
  print*,"example_z_uprkfact_randompencil:"
  print*, "K", k, "degree", dd 
  print*, "norm", norm
  print*,""
  write (7,*) "% K", k, "degree", dd, "norm", norm
  write (8,*) "% K", k, "degree", dd, "norm", norm
  write (9,*) "% K", k, "degree", dd, "norm", norm
  
  read (21, *) MA
  read (22, *) MB
  
  
  h = 0d0
  do ii=1,k
     do ll = 1, N
        h = h + abs(MA(ll,ii))**2
     end do
  end do
  
  norma = sqrt(h)
  
  h = 0d0
  do ii=1,k
     do ll = 1, N
        h = h + abs(MB(ll,ii))**2
     end do
  end do
  
  normb = sqrt(h) 
  
  
  MC = MA
  MD = MB

  print*, MA
  print*, MB
  
  
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

  MA = MC
  MB = MD
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Schur decomposition
  call system_clock(count=c_start3)

  if (david) then
     call z_uprk_compress2(.TRUE.,.TRUE.,.TRUE.,N,K,MA,MB,P,Q,&
          &D1,C1,B1,D2,C2,B2,V,W,INFO)
  else
     call z_uprk_compress(.TRUE.,.TRUE.,.TRUE.,N,K,MA,MB,P,Q,&
          &D1,C1,B1,D2,C2,B2,V,W,INFO)
  end if
  if (INFO.NE.0) then
     print*, "Info code from z_uprkdense_factor: ", INFO
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! check factorization
  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D1,C1,B1,TA)
  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D2,C2,B2,TB)

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
  
  h = 0d0
  do jj=1,N
     do ii=jj+1,N
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
     print*, ITS(ll)
  end do
  print*, "Iterations per eigenvalue: ", (1d0*it)/N

  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D1,C1,B1,TA)
  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D2,C2,B2,TB)
 
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


  print*,""
  deallocate(MA,MB,rev,iev,ITS)
  deallocate(Q,D1,C1,B1,P)
  deallocate(D2,C2,B2,V,W,EIGS,REIGS,REIGS2)
  deallocate(WORK,RWORK)
  deallocate(REIGSA,REIGSB,VL,VR,EIGSA,EIGSB)
  deallocate(CDA,CDB,MC,MD)
  deallocate(TA,TB,TL,Vt,Wt)

end program example_z_uprkfpen_2manyscaled
