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
  integer, parameter :: dd = 2
  integer, parameter :: k = 2
  integer :: N = dd*k
  integer :: ii, jj, ll, row, col, INFO
  complex(8), allocatable :: MA(:,:),MB(:,:), EIGS(:), REIGS(:),V(:,:),W(:,:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:), D(:), C(:), B(:)
  real(8), allocatable :: D2(:),C2(:),B2(:)
  logical, allocatable :: P(:)
  real(8) :: bulge(3), h
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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print banner
  print*,""
  print*,"example_z_uprkfact_randommatr:"
  print*,""
  
  allocate(MA(N,k),MB(N,k),rev(N),iev(N),ITS(N))

  allocate(Q(3*k*(N-1)),D(2*k*(N+1)),C(3*k*N),B(3*k*N),P(N-2))
  allocate(D2(2*k*(N+1)),C2(3*k*N),B2(3*k*N),V(N,N),W(N,N),EIGS(k*N),REIGS(N))

  call u_fixedseed_initialize(INFO)  
  call z_2Darray_random_normal(N,k,MA)


  ! make P0 lower triangular
  do ii=1,k-1
     do jj=ii+1,k
        MA(ii,jj)=cmplx(0d0,0d0,kind=8)
     end do
  end do
  MB = MA
  
  if ((dd==5).AND.(k==3)) then
     REIGS(1)  = cmplx(   2.237009273533627, - 2.916844980725226,kind=8)
     REIGS(2)  = cmplx(   2.796026659062130, - 1.152854347475400,kind=8)
     REIGS(3)  = cmplx(  -1.983058630091377, + 0.445143066361035,kind=8)
     REIGS(4)  = cmplx(  -1.191503202332451, - 0.842860669031502,kind=8)
     REIGS(5)  = cmplx(  -1.116234358401707, + 0.032598704706064,kind=8)
     REIGS(6)  = cmplx(  -0.734519195470930, - 0.619406960832573,kind=8)
     REIGS(7)  = cmplx(  -0.789064990432255, + 0.336502040275259,kind=8)
     REIGS(8)  = cmplx(   0.079436398189726, - 0.887867204932278,kind=8)
     REIGS(9)  = cmplx(   1.284438277711756, - 0.538201425121114,kind=8)
     REIGS(10) = cmplx(  -0.097451888579936, + 0.937452099676855,kind=8)
     REIGS(11) = cmplx(   0.856100145047217, + 0.816998912576558,kind=8)
     REIGS(12) = cmplx(   0.211733089095752, + 0.620582573065486,kind=8)
     REIGS(13) = cmplx(   1.056657627605067, + 0.315369650520287,kind=8)
     REIGS(14) = cmplx(   0.788587798555078, - 0.277842853565542,kind=8)
     REIGS(15) = cmplx(   0.374889615453730, - 0.239004891801360,kind=8)

     REIGS(01) = cmplx(  2.239774176279825d+00, - 2.910134671502209d+00,kind=8)
     REIGS(02) = cmplx(  2.784970929932215d+00, - 1.159503309215850d+00,kind=8)
     REIGS(03) = cmplx( -1.981205838915479d+00, + 4.891168311691939d-01,kind=8)
     REIGS(04) = cmplx( -1.262281333096642d+00, - 7.968143904565121d-01,kind=8)
     REIGS(05) = cmplx( -1.013843864858387d+00, + 6.072713557298853d-02,kind=8)
     REIGS(06) = cmplx( -7.566230668121174d-01, - 3.616859783263355d-01,kind=8)
     REIGS(07) = cmplx(  1.107989755971363d-02, + 9.310046448410413d-01,kind=8)
     REIGS(08) = cmplx(  8.170499806189339d-01, + 8.695299819599946d-01,kind=8)
     REIGS(09) = cmplx( -3.628293635224717d-01, + 6.501088596278545d-02,kind=8)
     REIGS(10) = cmplx( -1.133802105277646d-01, + 3.818352120826138d-01,kind=8)
     REIGS(11) = cmplx(  5.812212052488740d-02, - 7.940806993059377d-01,kind=8)
     REIGS(12) = cmplx(  1.293719264757970d+00, - 3.053653134775179d-01,kind=8)
     REIGS(13) = cmplx(  9.938843106681681d-01, + 1.822749609686237d-01,kind=8)
     REIGS(14) = cmplx(  6.943311112663597d-01, - 4.219604333521870d-01,kind=8)
     REIGS(15) = cmplx(  3.702785050702181d-01, - 2.001911432241459d-01,kind=8)
     do ii=1,N
        print*, MB(ii,1), MB(ii,2), MB(ii,3)
     end do
  end if
  if ((dd==2).AND.(k==2)) then
     REIGS(1) = cmplx( 5.850853916553567d-01, + 1.960620778250882d+00,kind=8)
     REIGS(2) = cmplx(-1.643275827808684d+00, + 9.249101310198671d-01,kind=8)
     REIGS(3) = cmplx(-1.665510376790535d-01, - 3.317227393069380d-01,kind=8)
     REIGS(4) = cmplx( 8.478546194285020d-01, - 5.774811996649782d-01,kind=8)
  end if

  if ((dd==2)) then
     if (k==1) then
        do ii=1,N
           print*, MB(ii,1)
        end do
     end if
     if (k==2) then
        do ii=1,N
           print*, MB(ii,1), MB(ii,2)
        end do
     end if
  end if

  if ((dd==5).AND.(k==1)) then
     REIGS(1) = cmplx(  1.249115269686787d0, - 0.994793987843677d0,kind=8)
     REIGS(2) = cmplx( -0.219189348479744d0, - 0.319360123437537d0,kind=8)
     REIGS(3) = cmplx( -1.244189839128987d0, + 0.399951930506511d0,kind=8)
     REIGS(4) = cmplx(  0.580709691243725d0, + 0.669645471844493d0,kind=8)
     REIGS(5) = cmplx( -0.092395247556379d0, + 1.167340958969769d0,kind=8)
     do ii=1,N
        print*, MB(ii,1)
     end do
  end if

  if ((dd==4).AND.(k==1)) then
     REIGS(1) = cmplx( -4.506367401502382d-01, + 8.258089507059133d-01,kind=8)
     REIGS(2) = cmplx(  5.665606540229318d-01, + 5.700065546795512d-01,kind=8)
     REIGS(3) = cmplx( -2.314228515628051d-01, - 3.194639432989243d-01,kind=8)
     REIGS(4) = cmplx(  7.234246907448657d-01, - 2.718456363764743d+00,kind=8)
     do ii=1,N
        print*, MB(ii,1)
     end do
  end if
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

  print*, "updated MB"
  do ii=1,N
     print*, MB(ii,:)
  end do

  ! initalize unitary-plus-rank-1 upper triangulars
  do ll = 1,k
     hc = MB(1,ll)
     MB(1:N-1,ll) = MB(2:N,ll) 
     MB(N,ll) = hc*(-1)**(N+1)
     call z_upr1_factoridpspike(.TRUE.,N,MB(1:N,ll),&
          &D(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C(((ll-1)*3*N+1):(ll*3*N)),&
          &B(((ll-1)*3*N+1):(ll*3*N)),INFO)
  end do


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
  !print*, "V"
  !do jj = 1,N
  !   print*, jj, V(jj,:)
  !end do

  !print*, MA-V(:,N-k+1:N)
  

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
           print*, col
           print*, "turnover", 3*(col-1)*(N-1)+3*(row-1)+1, 3*(col-1)*(N-1)+3*row,&
                & 3*(col-1)*(N-1)+3*(row-1)+4, 3*(col-1)*(N-1)+3*row+3
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
                 call z_uprkfact_rot3through1tri(.FALSE.,N,k,col,D,C,B,bulge,row)
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
        print*, ii, ll, EIGS((ll-1)*N+ii)         
        EIGS(ii) = EIGS(ii)*EIGS((ll-1)*N+ii)     
     end do
  end do

  do ii = 1,N
    jj = 1
    h = abs(EIGS(ii)-REIGS(jj))
    do ll = 2,N
       if (h>abs(EIGS(ii)-REIGS(ll))) then
          jj = ll
          h = abs(EIGS(ii)-REIGS(jj))
       end if
    end do
    print*, ii, EIGS(ii), REIGS(jj), abs(EIGS(ii)-REIGS(jj))
  end do
  

  print*, "The Eigenvalues"
  do ii = 1,N
    print*, ii,  REIGS(ii)
  end do

  print*,""
  deallocate(MA,MB,rev,iev,ITS)
  deallocate(Q,D,C,B,P)
  deallocate(D2,C2,B2,V,W,EIGS,REIGS)

end program example_z_uprkfact_randommatr
