#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_spr1_factor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine performs a Cayley transformation of the symmetric
! tridiagonal plus spike matrix to a unitary matrix plus spike
! (a descending, (I + ue_n^H), and a 
! ascending sequence of core transformations). The matrix is then
! transformed to upper Hessenberg form by core chasing.
!
! The Cayley (Moebius) transform -(z-i)/(z+i) maps the real line to
! the unit circle, in particular the interval [-1,1] is mapped to
! [-pi/2, pi/2]. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvectors
!                    .FALSE.: no schurvectors
!
!  ID              LOGICAL
!                    .TRUE.: initialize to Z to identity
!                    .FALSE.: assume Z is already initialized
!
!  SCA             LOGICAL
!                    .TRUE.: scale the matrix to 
!                            have eigenvalues in [-1,1]
!                    .FALSE.: do not scale the matrix
!                  !! CAUTION: Not scaling the matrix can 
!                              result in inaccurate eigenvalues !! 
! 
!  N               INTEGER
!                    dimension of matrix
!
!  D               REAL(8) array of dimension (N)  
!                    diagonal entries of T
!                    on exit: if SCA=.TRUE., D and E have been scaled
!
!  E               REAL(8) array of dimension (N-1)
!                    subdiagonal entries of T
!                    on exit: if SCA=.TRUE., D and E have been scaled
!
!  U               COMPLEX(8) array of dimension (N)
!                    the entries of the spike in the last column
!
!  M               INTEGER
!                    leading dimension of Z
!
! OUTPUT VARIABLES:
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for Givens rotations
!
!  QD              REAL(8) array of dimension (2*N+2)
!                    array of generators for complex diagonal matrix
!
!  QC, QB          REAL(8) array of dimension (3*N)
!                    array of generators for Givens rotations
!
!  SCALE           REAL(8) 
!                    scaling factor for back transform
!                    
!  Z               COMPLEX(8) array of dimension (M,N)
!                    similarity transformation of the reduction to 
!                    unitary Hessenberg form
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. and ID = .TRUE. initializes Z to I 
!                    if VEC = .TRUE. and ID = .FALSE. assumes Z initialized
!
!  INFO            INTEGER
!                    INFO = 1 implies scaling failed
!                    INFO = 0 implies successful computation
!                    INFO = -4 implies N is invalid
!                    INFO = -5 implies D is invalid
!                    INFO = -6 implies E is invalid
!                    INFO = -10 implies M is invalid
!                    INFO = -11 implies Z is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_spr1_factor2(VEC,ID,SCA,N,D,E,U,Q,QD,QC,QB,SCALE,M,Z,INFO)

  implicit none
  
  ! input variables
  logical, parameter :: DEBUGOUT = .FALSE.
  !logical, parameter :: DEBUGOUT = .TRUE.
  logical, intent(in) :: VEC, ID, SCA
  integer, intent(in) :: N, M
  real(8), intent(inout) :: D(N), E(N-1), Q(3*N-3), QD(2*N+2), SCALE
  real(8), intent(inout) :: QB(3*N), QC(3*N)
  complex(8), intent(inout) :: Z(M,N), U(N)
  integer, intent(inout) :: INFO

  ! compute variables
  integer :: ii, jj, ind1, ind2
  logical :: flg
  real(8) :: cr, ci, s, nrm, bulge(3), a, b, hb(3), LD(2*N), beta, phr, phi
  complex(8) :: eu, d1, t1(2,2), mu, temp, H(N,N),HQ(N,N), t(2,2)
  complex(8) :: WORK(5*N)
  real(8) :: RWORK(2*N)
  !print*, "START d_spr1_factor"
  ! initialize INFO
  INFO = 0

  !print*, "MITTE 1 d_spr1_factor"

  ! check N
  if (N < 1) then
    INFO = -4
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 1",INFO,INFO)
    end if
    return
  end if

  ! check D
  call d_1Darray_check(N,D,flg)
  if (.NOT.flg) then
    INFO = -5
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,INFO)
    end if
    return
  end if

  ! check E
  call d_1Darray_check(N-1,E,flg)
  if (.NOT.flg) then
    INFO = -6
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"E is invalid",INFO,INFO)
    end if
    return
  end if

  ! check M
  if (VEC.AND.(M < 1)) then
    INFO = -9
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"M must be at least 1",INFO,INFO)
    end if
    return
  end if
  
  ! check Z
  if (VEC.AND..NOT.ID) then
    call d_2Darray_check(M,N,Z,flg)
    if (.NOT.flg) then
      INFO = -10
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"Z is invalid",INFO,INFO)
      end if
      return
    end if
  end if   

  ! initialize Z
  if (VEC.AND.ID) then
     Z = cmplx(0d0,0d0,kind=8)
     do ii=1,min(M,N)
        Z(ii,ii) = cmplx(1d0,0d0,kind=8)
     end do
  end if

!!$  ! scaling
!!$  if (SCA) then
!!$     ! .TRUE. : use Newton correction
!!$     call d_symtrid_specint(.TRUE.,N,D,E,a,b,INFO)
!!$
!!$     if (INFO.NE.0) then 
!!$        ! print error in debug mode
!!$        if (DEBUG) then
!!$           call u_infocode_check(__FILE__,__LINE__,"d_symtrid_specint failed",INFO,INFO)
!!$        end if
!!$        INFO = 1
!!$
!!$     end if
!!$     
!!$     SCALE = max(abs(a),abs(b))
!!$     
!!$     do ii=1,N
!!$        D(ii) = D(ii)/SCALE
!!$     end do
!!$     do ii=1,N-1
!!$        E(ii) = E(ii)/SCALE
!!$     end do
!!$
!!$  end if

  !print*, "U",U

  !print*, "MITTE 2 d_spr1_factor"

  ! Cayley transform -(T-iI)(T+iI)              
  ! implicitly form B = T - i I  and T + i I = conjg(B)
  eu = E(1)
  d1 = cmplx(D(1),-1d0,kind=8)

  ! QR decomposition QR = (T-iI)
  do ii=1,N-1
     ! compute new rotation
     call z_rot3_vec3gen(dble(d1),aimag(d1),E(ii),cr,ci,s,nrm)
     
     !print*, nrm
     
     Q(3*ii-2) = cr
     Q(3*ii-1) = ci
     Q(3*ii)   = s

     ! update Q^H u
     t1(1,1) = cmplx(cr,-ci,kind=8)
     t1(2,1) = cmplx(-s,0d0,kind=8)
     t1(1,2) = -t1(2,1)
     t1(2,2) = conjg(t1(1,1))
     U(ii:(ii+1)) = matmul(t1,U(ii:(ii+1)))


     ! update banded matrix (ignoring the second superdiagonal and the upper part of the matrix)
     d1 = -s*eu + cmplx(cr,ci,kind=8)*cmplx(D(ii+1),-1d0,kind=8)
     if (ii<N-1) then
        eu = cmplx(cr*E(ii+1),ci*E(ii+1),kind=8)
     end if
  end do

  !print*, "MITTE 3 d_spr1_factor"

  ! form diagonal R conjg(R)^-1 (which is a diagonal matrix) 
  do ii=1,N+1
     QD(2*ii-1) = -1d0
     QD(2*ii)   = 0d0
  end do
  call d_rot2_vec2gen(dble(d1),aimag(d1),cr,s,nrm)
  !print*, "cr", cr, "s", s
  !call d_rot2_vec2gen(-cr*cr+s*s,-2*cr*s,QD(2*N-1),QD(2*N),nrm)             
  call d_rot2_vec2gen(-cr,-s,QD(2*N-1),QD(2*N),nrm)             
  ! update U
  U(N) = cmplx(cr,-s,kind=8)*U(N)
  mu = conjg(U(N))
  !mu = U(N)
  !print*, d1
  d1 = conjg(cmplx(cr,-s,kind=8)*d1)
 
!!$  print*, "Q'*U"
!!$  do ii=1,N
!!$     print*, ii,U(ii)
!!$  end do
!!$  print*, "rnn (d1)", d1, nrm
!!$  print*, "mu", mu
!!$  print*, "U-conjg(U))/(d1+mu)"
  ! update U
  do ii=1,N
     !U(ii) = 2d0*cmplx(0d0,aimag(U(ii)),kind=8)/(conjg(d1)+mu)
     !U(ii) = (U(ii)-conjg(U(ii)))/(conjg(d1)+mu)
     U(ii) = (U(ii)-conjg(U(ii)))/(d1+mu)!*cmplx(cr,+s,kind=8)
     !print*, ii,U(ii)
  end do
!!$  !U(N) = cmplx(cr,-s,kind=8)*U(N)
!!$  !U(N) = (U(N)-conjg(U(N)))/(d1+mu)
!!$  !print*, N,U(N)

  
!!$  print*, "Q ", Q
!!$  print*, "QD", QD

!!$  ! check Q, form HQ
!!$  HQ = 0d0
!!$  do ii=1,N
!!$     HQ(ii,ii) = 1d0
!!$  end do
!!$  ! Apply QD from the left
!!$  do ii=1,N
!!$     HQ(ii,:) = cmplx(QD(2*ii-1),QD(2*ii),kind=8)*HQ(ii,:)
!!$  end do
!!$  ! Apply Q from the left
!!$  do ii=N-1,1,-1
!!$     t(1,1) = cmplx(Q(3*ii-2),Q(3*ii-1),kind=8)
!!$     t(2,1) = cmplx(Q(3*ii),0d0,kind=8)
!!$     t(1,2) = -t(2,1)
!!$     t(2,2) = conjg(t(1,1))
!!$     HQ(ii:(ii+1),:) = matmul(t,HQ(ii:(ii+1),:))
!!$  end do 
!!$  ! plot H
!!$  print*, "Q"
!!$  do ii=1,N
!!$     print*, ii, HQ(ii,1:N)
!!$  end do

  if (DEBUGOUT) then
     
     ! check H, form H
     H = 0d0
     do ii=1,N
        H(ii,ii) = 1d0
     end do
     do ii=1,N
        H(ii,N) = H(ii,N)+U(ii)
        print*, "H U", ii, H(ii,N), U(ii)
     end do
!!$  print*, "RiR"
!!$  do ii=1,N
!!$     print*, ii, H(ii,1:N)
!!$  end do
     
     
     ! Apply QD from the left
     do ii=1,N
        H(ii,:) = cmplx(QD(2*ii-1),QD(2*ii),kind=8)*H(ii,:)
     end do
     ! Apply Q from the left
     do ii=N-1,1,-1
        t(1,1) = cmplx(Q(3*ii-2),Q(3*ii-1),kind=8)
        t(2,1) = cmplx(Q(3*ii),0d0,kind=8)
        t(1,2) = -t(2,1)
        t(2,2) = conjg(t(1,1))
        H(ii:(ii+1),:) = matmul(t,H(ii:(ii+1),:))
     end do
     
!!$  print*, "Q*RiR"
!!$  do ii=1,N
!!$     print*, ii, H(ii,1:N)
!!$  end do
     ! Apply QD^T from the right
     do ii=1,N
        H(:,ii) = -H(:,ii)*cmplx(QD(2*ii-1),QD(2*ii),kind=8)
     end do
     ! Apply Q^T from the right
     do ii=N-1,1,-1
        t(1,1) = cmplx(Q(3*ii-2),Q(3*ii-1),kind=8)
        t(2,1) = cmplx(-Q(3*ii),0d0,kind=8)
        t(1,2) = -t(2,1)
        t(2,2) = conjg(t(1,1))
        H(:,ii:(ii+1)) = matmul(H(:,ii:(ii+1)),t)
     end do
     ! plot H
!!$  print*, "H"
!!$  do ii=1,N
!!$     print*, ii, H(ii,1:N)
!!$  end do
     
     print*, "First ZGEEV"
     call zgeev('N','N', N, H, N, HQ, Z, 1, Z, 1, WORK, 5*N, RWORK, INFO)
     
     ! EIGENVALUES ARE CORRECT HERE
  
     do ii=1,N
        print*, ii, HQ(ii,1), cmplx(0d0,1d0,kind=8)*(cmplx(1d0,0d0,kind=8)-HQ(ii,1))/&
             &(cmplx(1d0,0d0,kind=8)+HQ(ii,1))
     end do
  end if










  !U = -conjg(U)
  !U = -U
  !U = conjg(U)

  QC = 0d0
  QB = 0d0

  ! roll up U
  ! compute the phase of last coefficient
  call d_rot2_vec2gen(dble(U(N))+1,aimag(U(N)),phr,phi,beta)
 
  if (DEBUGOUT) then

     print*,""
     print*,"U(N):",U(N)
     print*,"beta:",beta
     print*,""
     
     
     print*, "QD",QD
     print*, "ph", phr, phi
     
  end if

  ! store in QD
  t1(1,1) = cmplx(QD(2*N-1),QD(2*N),kind=8)*cmplx(phr,phi,kind=8)
  call d_rot2_vec2gen(dble(t1(1,1)),aimag(t1(1,1)),QD(2*N-1),QD(2*N),nrm)
  t1(1,1) = cmplx(QD(2*N+1),QD(2*N+2),kind=8)*cmplx(phr,-phi,kind=8)
  call d_rot2_vec2gen(dble(t1(1,1)),aimag(t1(1,1)),QD(2*N+1),QD(2*N+2),nrm)

  if (DEBUGOUT) then
     print*, "QD",QD
  end if

  ! initialize bottom of QC
  call d_rot2_vec2gen(beta,-1d0,QC(3*N-2),QC(3*N),nrm)

  !!$ call z_rot3_vec3gen(dble(U(N)),aimag(U(N)),beta,QC(3*N-2),QC(3*N-1),QC(3*N),nrm)

  ! initialize bottom of QB
  QB(3*N-2) = QC(3*N)
  QB(3*N) = QC(3*N-2)

  ! roll up U into QB and QC
  temp = cmplx(nrm,0d0,kind=8)
  do ii = 1,(N-1)
    ! compute new QC
    call z_rot3_vec4gen(dble(U(N-ii)),aimag(U(N-ii)),dble(temp),aimag(temp) &
         &,QC(3*(N-ii)-2),QC(3*(N-ii)-1),QC(3*(N-ii)),nrm)

    ! store new QB
    QB(3*(N-ii)-2) = QC(3*(N-ii)-2)
    QB(3*(N-ii)-1) = -QC(3*(N-ii)-1)
    QB(3*(N-ii)) = -QC(3*(N-ii))

    ! update last entry of U using QC
    temp = cmplx(QC(3*(N-ii)-2),-QC(3*(N-ii)-1),kind=8)*U(N-ii) &
         &+ QC(3*(N-ii))*temp

  end do







  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! chasing upward
  ! initialize LD
  do ii=1,N
     LD(2*ii-1)= 1d0
     LD(2*ii) = 0d0
  end do
  
  ! fuse Q^T into Q
  do jj=1,N-1
     !print*, "MITTE 4a d_spr1_factor",jj
     bulge(1) = Q(3*jj-2)
     bulge(2) = Q(3*jj-1)
     bulge(3) = -Q(3*jj)  ! the bulge is Q^T
     
     ! update eigenvectors
     if (VEC) then
        t1(1,1) = cmplx(bulge(1),-bulge(2),kind=8)
        t1(2,1) = cmplx(-bulge(3),0d0,kind=8)
        t1(1,2) = -t1(2,1)
        t1(2,2) = conjg(t1(1,1))
        Z(:,(jj):(jj+1)) = matmul(Z(:,(jj):(jj+1)),t1)
     end if
     
     ! set indices
     ind1 = 2*jj - 1
     ind2 = ind1 + 3
     
     ! through diag 
     call z_rot3_swapdiag(.TRUE.,LD(ind1:ind2),bulge)
     
     ! main chasing loop
     do ii=jj,2,-1
        
        !print*, "MITTE 4b d_spr1_factor",jj,ii
        ! set indices
        ind1 = 3*(ii-2) + 1
        ind2 = ind1+2
        
        ! through Q
        call z_rot3_turnover(bulge,Q(ind1:ind2),Q((ind1+3):(ind2+3)))
        hb = bulge
        bulge = Q(ind1:ind2)
        Q(ind1:ind2) = Q((ind1+3):(ind2+3))
        Q((ind1+3):(ind2+3)) = hb
        
        ! set indices
        ind1 = 2*(ii-2) + 1
        ind2 = 3*(ii-2) + 1
        
        ! through diag 
        !call z_rot3_swapdiag(.TRUE.,QD(ind1:ind2),bulge)
        call z_upr1fact_rot3throughtri(.TRUE.,QD(ind1:ind1+3),QC(ind2:ind2+5),QB(ind2:ind2+5),bulge)
        
        !print*, "MITTE 4c d_spr1_factor",jj,ii
        
        ! update U
!!$        t1(1,1) = cmplx(bulge(1),bulge(2),kind=8)
!!$        t1(2,1) = cmplx(bulge(3),0d0,kind=8)
!!$        t1(1,2) = -t1(2,1)
!!$        t1(2,2) = conjg(t1(1,1))
!!$        U((ii-1):(ii)) = matmul(t1,U((ii-1):(ii)))        

!!$        print*,"update U", U(ii-1), U(ii)
!!$        print*, cmplx(bulge(1),bulge(2))*U(ii-1), cmplx(bulge(3),0d0,kind=8)*U(ii)
!!$        print*, cmplx(bulge(1),-bulge(2))*U(ii), -cmplx(bulge(3),0d0,kind=8)*U(ii-1)
!!$        t1(1,1) = cmplx(bulge(1),bulge(2))*U(ii-1) - &
!!$             &cmplx(bulge(3),0d0,kind=8)*U(ii)
!!$        U(ii) = cmplx(bulge(1),-bulge(2))*U(ii) + &
!!$             &cmplx(bulge(3),0d0,kind=8)*U(ii-1)
!!$        U(ii-1) = t1(1,1)
!!$        print*,"update U", U(ii-1), U(ii)
!!$
!!$

        ! update eigenvectors
        if (VEC) then
           t1(1,1) = cmplx(bulge(1),-bulge(2),kind=8)
           t1(2,1) = cmplx(-bulge(3),0d0,kind=8)
           t1(1,2) = -t1(2,1)
           t1(2,2) = conjg(t1(1,1))
           Z(:,(ii-1):(ii)) = matmul(Z(:,(ii-1):(ii)),t1)
        end if
        
        ! through diag 
        call z_rot3_swapdiag(.TRUE.,LD(ind1:ind2),bulge)
        
     end do
     !print*, "MITTE 4d d_spr1_factor",jj
     
     ! fusion at bottom
     call z_unifact_mergebulge(.TRUE.,Q(1:3),QD(1:4),bulge)
     
     !print*, "MITTE 4e d_spr1_factor",jj
     ! update LD(1:2)
     t1(1,1) = cmplx(LD(1),LD(2),kind=8)*cmplx(bulge(1),bulge(2),kind=8)
     call d_rot2_vec2gen(dble(t1(1,1)),aimag(t1(1,1)),LD(1),LD(2),nrm)
     
     !print*, "MITTE 4f d_spr1_factor",jj
     ! update LD(3:4)
     t1(1,1) = cmplx(LD(3),LD(4),kind=8)*cmplx(bulge(1),-bulge(2),kind=8)
     call d_rot2_vec2gen(dble(t1(1,1)),aimag(t1(1,1)),LD(3),LD(4),nrm)
     !print*, "MITTE 4g d_spr1_factor",jj
     
  end do

  !print*, "LD(N)", LD(2*N-1), LD(2*N)
  t1(1,1) = cmplx( LD(2*N-1), LD(2*N), kind=8) * cmplx(cr,s,kind=8);
  if (VEC) then
     t1(1,1) = cmplx(-cr,s,kind=8)
     Z(:,N) = Z(:,N)*t1(1,1)
  end if
  call d_rot2_vec2gen(dble(t1(1,1)),aimag(t1(1,1)),LD(2*N-1),LD(2*N),nrm)
  !print*, "LD(N)", LD(2*N-1), LD(2*N)


  ! merge LD in QD
  do ii=1,N-1
     ind1 = 2*ii-1
     ind2 = 3*ii-2
     call z_rot3_swapdiag(.FALSE.,LD(ind1:(ind1+3)),Q(ind2:(ind2+2)))
  end do
  do ii=1,N
     t1(1,1) = cmplx(LD(2*ii-1),LD(2*ii),kind=8)
     ! update U
     !U(ii) = conjg(t1(1,1))*U(ii)
     ! update eigenvectors
     !if (VEC) then
     !   Z(:,ii) = Z(:,ii)*t1(1,1)
     !end if
     
     t1(1,1) = cmplx(QD(2*ii-1),QD(2*ii),kind=8)*t1(1,1)
     call d_rot2_vec2gen(dble(t1(1,1)),aimag(t1(1,1)),QD(2*ii-1),QD(2*ii),nrm)
     
  end do









  if (DEBUGOUT) then


  ! check H, form H
  ! compute matrix H
  H = 0d0
  call z_upr1fact_extracttri(.FALSE.,N,QD,QC,QB,H)
  ! plot H
  do ii=1,N
     print*, ii, H(ii,:)
  end do

  do ii=N-1,1,-1
     t(1,1) = cmplx(Q(3*ii-2),Q(3*ii-1),kind=8)
     t(2,1) = cmplx(Q(3*ii),0d0,kind=8)
     t(1,2) = -t(2,1)
     t(2,2) = conjg(t(1,1))
     H(ii:(ii+1),:) = matmul(t,H(ii:(ii+1),:))
  end do
  print*, "reconstructed H"
  ! plot H
  do ii=1,N
     print*, ii, H(ii,:)
  end do

  print*, "Second ZGEEV"
  call zgeev('N','N', N, H, N, HQ, Z, 1, Z, 1, WORK, 5*N, RWORK, INFO)
  
  ! EIGENVALUES ARE PERTURBED HERE
  
  do ii=1,N
     print*, ii, HQ(ii,1), cmplx(0d0,1d0,kind=8)*(cmplx(1d0,0d0,kind=8)-HQ(ii,1))/&
          &(cmplx(1d0,0d0,kind=8)+HQ(ii,1))
  end do




end if
  


end subroutine d_spr1_factor2
