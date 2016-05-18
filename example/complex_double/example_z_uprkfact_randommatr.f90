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
  complex(8), allocatable :: MA(:,:),MB(:,:), EIGS(:), REIGS(:)
  complex(8), allocatable :: REIGS2(:), V(:,:),W(:,:), WORK(:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:), D(:), C(:), B(:)
  real(8), allocatable :: D2(:),C2(:),B2(:), RWORK(:)
  logical, allocatable :: P(:)
  real(8) :: h, maxerr

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
  print*,"example_z_uprkfact_randommatr:"
  print*,""
  
  allocate(MA(N,k),MB(N,k),rev(N),iev(N),ITS(N))

  allocate(Q(3*k*(N-1)),D(2*k*(N+1)),C(3*k*N),B(3*k*N),P(N-2))
  allocate(D2(2*k*(N+1)),C2(3*k*N),B2(3*k*N),V(N,N),W(N,N),EIGS(k*N),REIGS(N))
  allocate(WORK(25*N+25),RWORK(2*N),REIGS2(N))

  call u_fixedseed_initialize(INFO)  
  call z_2Darray_random_normal(N,k,MA)


  ! make P0 lower triangular
  do ii=1,k-1
     do jj=ii+1,k
        MA(ii,jj)=cmplx(0d0,0d0,kind=8)
     end do
  end do
  MB = MA
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute roots with LAPACK
  call system_clock(count=c_start2)
  V = cmplx(0d0,0d0,kind=8)

  do ii=1,N-k
     V(ii+k,ii)=cmplx(1d0,0d0,kind=8)
  end do
  V(:,N-k+1:N) = MA

  lwork = 25*N+25
  call zgeev('N','N', N, V, N, REIGS, W, N, W, N, WORK, lwork, RWORK, INFO)
  
  if (output) then
     do jj= 1,N
        print*, jj, REIGS(jj)
     end do
  end if
  call system_clock(count=c_stop2)
  print*, "Runtime unstructured QR solver (LAPACK) ",(dble(c_stop2-c_start2)/dble(c_rate))

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Start Times
  call system_clock(count=c_start)

  call z_uprkdense_qz(.FALSE.,.FALSE.,N,k,MA,MB,EIGS,EIGS,V,W,INFO)

!!$  call z_uprkdense_factor(.FALSE.,.FALSE.,.FALSE.,N,k,MA,MB,P,Q,D,C,B,D2,C2,B2,V,W,INFO)
!!$  if (INFO.NE.0) then
!!$     print*, "Info code from z_uprkdense_factor: ", INFO
!!$  end if
!!$
!!$  call z_uprkfact_twistedqz(.FALSE.,.FALSE.,.FALSE.,l_upr1fact_upperhess,N,k,&
!!$       &P,Q,D,C,B,D2,C2,B2,V,W,ITS,INFO)
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
!!$          &D(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C(((ll-1)*3*N+1):(ll*3*N)),&
!!$          &B(((ll-1)*3*N+1):(ll*3*N)),V)
!!$     EIGS(((ll-1)*N+1):ll*N) = V(1:N,1)
!!$  end do
!!$  do ll = 2,k
!!$     do ii=1,N
!!$        !print*, ii, ll, EIGS((ll-1)*N+ii)         
!!$        EIGS(ii) = EIGS(ii)*EIGS((ll-1)*N+ii)     
!!$     end do
!!$  end do

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
  deallocate(Q,D,C,B,P)
  deallocate(D2,C2,B2,V,W,EIGS,REIGS,REIGS2)
  deallocate(WORK,RWORK)

end program example_z_uprkfact_randommatr
