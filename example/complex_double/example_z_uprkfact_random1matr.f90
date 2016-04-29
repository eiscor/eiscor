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
program example_z_uprkfact_randommatr2


  implicit none
  
  ! compute variables
  integer, parameter :: dd = 5
  integer, parameter :: k = 1
  integer :: N = dd*k
  integer :: ii, jj, ll, row, col, INFO
  complex(8), allocatable :: MA(:,:),MB(:,:), EIGS(:), REIGS(:),V(:,:),W(:,:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:), D(:), C(:), B(:)
  real(8), allocatable :: D2(:),C2(:),B2(:)
  logical, allocatable :: P(:)
  real(8) :: bulge(3), h
  complex(8) :: Gf(2,2), hc

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
     do ii=1,N
        print*, MB(ii,1), MB(ii,2), MB(ii,3)
     end do
  end if
  if ((dd==2).AND.(k==2)) then
     REIGS(1) = cmplx(-1.898605271153636, + 0.896922815943011,kind=8)
     REIGS(2) = cmplx( 0.578128456518313, + 2.195763258929064,kind=8)
     REIGS(3) = cmplx(-0.054246276942936, - 0.703642061460719,kind=8)
     REIGS(4) = cmplx( 0.997836237174379, - 0.412717043112522,kind=8)
     do ii=1,N
        print*, MB(ii,1), MB(ii,2)
     end do
  end if

  if ((dd==5).AND.(k==1)) then
     REIGS(1) = cmplx(  1.249115269686787, - 0.994793987843677,kind=8)
     REIGS(2) = cmplx( -0.219189348479744, - 0.319360123437537,kind=8)
     REIGS(3) = cmplx( -1.244189839128987, + 0.399951930506511,kind=8)
     REIGS(4) = cmplx(  0.580709691243725, + 0.669645471844493,kind=8)
     REIGS(5) = cmplx( -0.092395247556379, + 1.167340958969769,kind=8)
     do ii=1,N
        print*, MB(ii,1)
     end do
  end if
  ! split MA into k upper unitary plus rank 1 upper Hessenberg matrices with one
  ! non-zero, non-unity column in the last column (stored in MB)
  do ii = 1,k
     do jj = ii+1,k
        MB(n,jj) = MB(1,jj)/MB(1,ii);
     end do
     do ll = ii+1,k
        do jj = 1,n-1
           MB(jj,ll) = MB(jj,ll)-MB(jj,ii)*MB(n,ll);
        end do
     end do
  end do
  
  ! fill P
  P = .FALSE.

  ! initialize Q
  Q = 0d0
  do ii = 1,k*(N-1)
    Q(3*ii) = 1d0
  end do

  ! initalize unitary-plus-rank-1 upper triangulars
  do ll = 1,k
     hc = MB(1,ll)
     MB(1:N-1,ll) = MB(2:N,ll) 
     MB(N,ll) = hc
     call z_upr1_factoridpspike(.TRUE.,N,MB(1:N,ll),&
          &D(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C(((ll-1)*3*N+1):(ll*3*N)),&
          &B(((ll-1)*3*N+1):(ll*3*N)),INFO)
  end do
  

  ! reduce to Hessenberg * Triangular * ... * Triangular form
  ! remove rotations from Q_ii starting with k
  do ii = k,2,-1
     ! remove rotation jj
     do jj = 1,(N-1)
        col = ii-1 ! rotation right of col th upper triangular
        bulge = Q(3*(ii-1)*(N-1)+3*(jj-1)+1:3*(ii-1)*(N-1)+3*jj)
        do row = jj,N-2 ! first row of the rotation
           ! pass through upper triangular
           call z_uprkfact_rot3through1tri(.FALSE.,N,k,col,D,C,B,bulge,row)
           !call z_upr1fact_rot3throughtri(.FALSE.,&
           !     &D(((col-1)*2*(N+1)+(row-1)*2+1):((col-1)*2*(N+1)+row*2)),&
           !     &C(((col-1)*3*N+1):(col*3*N)),&
           !     &B(((col-1)*3*N+1):(col*3*N)),bulge)
           ! turnover through Q
           call z_rot3_turnover(&
                &Q(3*(col-1)*(N-1)+3*(row-1)+1:3*(col-1)*(N-1)+3*row),&
                &Q(3*(col-1)*(N-1)+3*(row-1)+1:3*(col-1)*(N-1)+3*row),bulge)
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
        call z_upr1fact_mergebulge(.FALSE.,N,.FALSE.,&
             &Q(((col-1)*3*(N-1)+1):(col*3*(N-1))),&
             &D(((col-1)*2*(N+1)+1):(col*2*(N+1))),bulge)
     end do
  end do


  do ll = 1,k
     call z_upr1fact_extracttri(.false.,N,&
          &D(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C(((ll-1)*3*N+1):(ll*3*N)),&
          &B(((ll-1)*3*N+1):(ll*3*N)),V)
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

  end do

  call z_upr1fact_twistedqz(.FALSE.,.FALSE.,.FALSE.,l_upr1fact_upperhess,N,&
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
  

  print*,""
  deallocate(MA,MB,rev,iev,ITS)
  deallocate(Q,D,C,B,P)
  deallocate(D2,C2,B2,V,W,EIGS,REIGS)

end program example_z_uprkfact_randommatr2
