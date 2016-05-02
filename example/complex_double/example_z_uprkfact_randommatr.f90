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
  integer, parameter :: dd = 30
  integer, parameter :: k = 20
  logical, parameter :: output=.FALSE.
  !logical, parameter :: output=.TRUE.
  integer :: N = dd*k
  integer :: ii, jj, ll, row, col, INFO, lwork
  complex(8), allocatable :: MA(:,:),MB(:,:), EIGS(:), REIGS(:)
  complex(8), allocatable :: REIGS2(:), V(:,:),W(:,:), WORK(:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:), D(:), C(:), B(:)
  real(8), allocatable :: D2(:),C2(:),B2(:), RWORK(:)
  logical, allocatable :: P(:)
  real(8) :: bulge(3), h, maxerr
  complex(8) :: hc, Gf(2,2)

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
  allocate(WORK(5*N),RWORK(2*N),REIGS2(N))

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

  lwork = 5*N
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

  ! split MA into k upper unitary plus rank 1 upper Hessenberg matrices with one
  ! non-zero, non-unity column in the last column (stored in MB)
  !do ii = 1,k
  !   do jj = ii+1,k
  !      MB(n,jj) = MB(1,jj)/MB(1,ii);
  !   end do
  !   do ll = ii+1,k
  !      do jj = 1,n-1
  !         MB(jj,ll) = MB(jj+1,ll)!-MB(jj+1,ii)*MB(n,ll);
  !      end do
  !   end do
  !end do
  do ii = 2,k
     do jj = 1,(N-ii+1)
        MB(jj,ii) = MB(jj+ii-1,ii)
     end do
     do jj = (N-ii+2),N
        MB(jj,ii) = cmplx(0d0,0d0,kind=8)
     end do
  end do

  
  ! fill P
  P = .FALSE.

  ! initialize Q
  Q = 0d0
  do ii = 1,k*(N-1)
    Q(3*ii) = 1d0
  end do

  if (output) then     
     print*, "updated MB"
     do ii=1,N
        print*, MB(ii,:)
     end do
  end if

  ! initalize unitary-plus-rank-1 upper triangulars
  do ll = 1,k
     hc = MB(1,ll)
     MB(1:N-1,ll) = MB(2:N,ll) 
     MB(N,ll) = hc*(-1)**(N+1)
     call z_upr1_factoridpspike(.TRUE.,N,MB(1:N,ll),&
          &D(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C(((ll-1)*3*N+1):(ll*3*N)),&
          &B(((ll-1)*3*N+1):(ll*3*N)),INFO)
  end do


  if (output) then
     ! check factorization
     ! Form upper triangulars and multiply them together
     V = cmplx(0d0,0d0,kind=8)
     W = cmplx(0d0,0d0,kind=8)
     do ll = 1,N
        V(ll,ll) = cmplx(1d0,0d0,kind=8)
     end do
     do ll = 1,k
        call z_upr1fact_extracttri(.false.,N,&
             &D(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C(((ll-1)*3*N+1):(ll*3*N)),&
             &B(((ll-1)*3*N+1):(ll*3*N)),W)
        ! Apply Q
        do jj = N-1,1,-1
           ! Apply Q(jj)
           Gf(1,1) = cmplx(Q((ll-1)*(N-1)+3*jj-2),Q((ll-1)*(N-1)+3*jj-1),kind=8)
           Gf(2,1) = cmplx(Q((ll-1)*(N-1)+3*jj),0d0,kind=8)
           Gf(1,2) = -Gf(2,1)
           Gf(2,2) = conjg(Gf(1,1))
           
           W((jj):(jj+1),:) = matmul(Gf,W((jj):(jj+1),:))
        end do
        
        ! print*, "W", ll
        ! do jj = 1,N
        !    print*, jj, W(jj,:)
        ! end do
        
        V = matmul(V,W)
     end do
     
     print*, "V"
     do jj = 1,N
        print*, jj, V(jj,:)
     end do
  end if 

  ! reduce to Hessenberg * Triangular * ... * Triangular form
  ! remove rotations from Q_ii starting with k
  do ii = k,2,-1
     ! remove rotation jj
     do jj = 1,(N-1)
        col = ii-1 ! rotation right of col th upper triangular
        !print*, col
        !print*, "bulge", (3*col*(N-1)+3*(jj-1)+1),(3*col*(N-1)+3*jj)
        bulge = Q((3*col*(N-1)+3*(jj-1)+1):(3*col*(N-1)+3*jj))
        do row = jj,N-2 ! first row of the rotation
           ! pass through upper triangular
           call z_uprkfact_rot3through1tri(.FALSE.,N,k,col,D,C,B,bulge,row)
           !call z_upr1fact_rot3throughtri(.FALSE.,&
           !     &D(((col-1)*2*(N+1)+(row-1)*2+1):((col-1)*2*(N+1)+row*2)),&
           !     &C(((col-1)*3*N+1):(col*3*N)),&
           !     &B(((col-1)*3*N+1):(col*3*N)),bulge)
           ! turnover through Q
           !print*, col
           !print*, "turnover", 3*(col-1)*(N-1)+3*(row-1)+1, 3*(col-1)*(N-1)+3*row,&
           !     & 3*(col-1)*(N-1)+3*(row-1)+4, 3*(col-1)*(N-1)+3*row+3
           call z_rot3_turnover(&
                &Q((3*(col-1)*(N-1)+3*(row-1)+1):(3*(col-1)*(N-1)+3*row)),&
                &Q((3*(col-1)*(N-1)+3*(row-1)+4):(3*(col-1)*(N-1)+3*row+3)),bulge)
           ! update col
           col = col - 1
           ! update col if left most sequence was passed
           if (col==0) then
              ! update eigenvectors 
              !!!!!!!!!!!!!!!!!!!!!!!!!
              ! not implemented
              !!!!!!!!!!!!!!!!!!!!!!!!!
              ! pass through triangulars without rotations in front
              do col = k,ii+1,-1
                 !print*, "just triangle col", col, "row", row, "ii", ii, "jj", jj
                 call z_uprkfact_rot3through1tri(.FALSE.,N,k,col,D,C,B,bulge,row+1)
              end do
              !print*, "after update col", col, "ii", ii, "jj", jj
           end if
        end do
        !print*, "fusion row", row, "N", N, "col", col
        ! fuse rotation into current rotation
        call z_uprkfact_rot3through1tri(.FALSE.,N,k,col,D,C,B,bulge,row)
        !call z_upr1fact_mergebulge(.FALSE.,N,.FALSE.,&
        !     &Q(((col-1)*3*(N-1)+1):(col*3*(N-1))),&
        !     &D(((col-1)*2*(N+1)+1):(col*2*(N+1))),bulge)
        call z_upr1fact_mergebulge(.FALSE.,N,.FALSE.,&
             &Q(((col-1)*3*(N-1)+1):(col*3*(N-1))),&
             &D(((col-1)*2*(N+1)+1):((col-1)*2*(N+1)+2*N)),bulge)
     end do
  end do

  if (output) then
     ! check factorization
     ! Form upper triangulars and multiply them together
     V = cmplx(0d0,0d0,kind=8)
     W = cmplx(0d0,0d0,kind=8)
     do ll = 1,N
        V(ll,ll) = cmplx(1d0,0d0,kind=8)
     end do
     do ll = 1,k
        call z_upr1fact_extracttri(.false.,N,&
             &D(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C(((ll-1)*3*N+1):(ll*3*N)),&
             &B(((ll-1)*3*N+1):(ll*3*N)),W)
        V = matmul(V,W)
     end do
     ! Apply Q
     do jj = N-1,1,-1
        ! Apply Q(jj)
        Gf(1,1) = cmplx(Q(3*jj-2),Q(3*jj-1),kind=8)
        Gf(2,1) = cmplx(Q(3*jj),0d0,kind=8)
        Gf(1,2) = -Gf(2,1)
        Gf(2,2) = conjg(Gf(1,1))
        
        V((jj):(jj+1),:) = matmul(Gf,V((jj):(jj+1),:))
     end do
     
     do jj = 1,N
        print*, jj, V(jj,:)
     end do

     call zgeev('N','N', N, V, N, REIGS2, W, N, W, N, WORK, lwork, RWORK, INFO)
     
     print*, "Eigenvalues after reduction to upper Hessenberg times triangular form"
     do jj = 1,N
        print*, REIGS(jj),REIGS2(jj)
     end do
     
  end if

  

  call z_uprkfact_twistedqz(.FALSE.,.FALSE.,.FALSE.,l_upr1fact_upperhess,N,k,&
       &P,Q,D,C,B,D2,C2,B2,V,W,ITS,INFO)

  do ll = 1,k
     call z_upr1fact_extracttri(.TRUE.,N,&
          &D(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C(((ll-1)*3*N+1):(ll*3*N)),&
          &B(((ll-1)*3*N+1):(ll*3*N)),V)
     EIGS(((ll-1)*N+1):ll*N) = V(1:N,1)
  end do
  do ll = 2,k
     do ii=1,N
        !print*, ii, ll, EIGS((ll-1)*N+ii)         
        EIGS(ii) = EIGS(ii)*EIGS((ll-1)*N+ii)     
     end do
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
  print*, "Runtime structured QR solver ",(dble(c_stop-c_start)/dble(c_rate))
  print*, "Runtime unstructured QR solver (LAPACK) ",(dble(c_stop2-c_start2)/dble(c_rate))

  

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
