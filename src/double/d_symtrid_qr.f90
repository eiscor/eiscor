#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_symtrid_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the real Schur factorization of a real 
! symmetric tridiagonal matrix T.
!
! This routine performs a Moebius transformation of the symmetric
! tridiagonal matrix to a unitary matrix (a descending and a 
! ascending sequence of core transformations). The unitary matrix
! is transformed to upper Hessenberg form and then passed to 
! z_unifact_qr.
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
!  N               INTEGER
!                    dimension of matrix
!
!  D               REAL(8) array of dimension (N)  
!                    diagonal entries of T
!                    on exit D contains the eigenvalues of T
!  E               REAL(8) array of dimension (N)
!                    subdiagonal entries of T
!
!  WORK            REAL(8) array of dimension (14*N)
!                    work space for eigensolver
!
!  M               INTEGER
!                    leading dimension of Z
!
! OUTPUT VARIABLES:
!
!  Z               COMPLEX(8) array of dimension (M,N)
!                    components of schurvectors
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. and ID = .TRUE. initializes Z to I 
!                    if VEC = .TRUE. and ID = .FALSE. assumes Z initialized
!
!  ITS             INTEGER array of dimension (N-1)
!                    Contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO = 2 implies  failed
!                    INFO = 1 implies  failed
!                    INFO = 0 implies successful computation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_symtrid_qr(VEC,ID,N,D,E,WORK,M,Z,ITS,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, ID
  integer, intent(in) :: N, M
  real(8), intent(inout) :: WORK(14*N)
  integer, intent(inout) :: ITS(N-1), INFO
  real(8), intent(inout) :: D(N), E(N)
  complex(8), intent(inout) :: Z(M,N)
  ! compute variables
  integer :: ii, jj, ind1, ind2
  real(8) :: cr,ci,s,nrm

  real(8) ::   tha, thb, thg, thd, ag, ad
  complex(8) :: ca, cb, cg, cd, temp
  
  complex(8) :: D1(N), EU(N)
  real(8) :: QB(3*N-3),QC(3*N-3), DR(2*N), bulge(3), pi = EISCOR_DBL_PI
 
  !!$  complex(8) :: Z2(M,N), Z2H(M,N), H(N,N), H2(N,N)
 
   
  complex(8) :: block(2,2), t1(2,2), t2(2,2)
  
  ! initialize INFO
  INFO = 0

  if (VEC.AND.ID) then
     print*, "Reset Z"
    Z = cmplx(0d0,0d0,kind=8)
    do ii=1,min(M,N)
      Z(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if
  
  ! parameters Moebius transformation
  ca = cmplx(1d0,0d0,kind=8)
  cb = cmplx(0d0,1d0,kind=8)
  cg = cmplx(1d0,0d0,kind=8)
  cd = cmplx(0d0,-1d0,kind=8)

  ! form B = a T + b I
  ! form C = g T + d I
  ! ii = 1,...,N-1
  do ii=1,N-1
     !D1(ii) = ca * D(ii) + cb
     !D2(ii) = cg * D(ii) + cd
     ! real
     D1(ii) = cmplx(D(ii),1d0,kind=8)
     !D2(ii) = cmplx(D(ii),-1d0,kind=8)
     ! D2 = conjg(D1)
     
     ! EL == E
     ! EL(ii) = 1d0 * E(ii) ! now real
     EU(ii) = ca * E(ii)
     ! ca == cg == 1
     !EL(ii) = cg * E(ii)
     !EU(ii) = cg * E(ii)
  end do
  ! ii = N
  D1(N) = cmplx(D(N),1d0,kind=8)
  !D2(N) = cmplx(D(N),-1d0,kind=8)

  ! QR decomposition B = QB Rb
  do ii=1,N-1
     ! compute new rotation
     !call z_rot3_vec4gen(dble(D1(ii)),aimag(D1(ii)),dble(EL(ii)),aimag(EL(ii)),cr,ci,s,nrm)
     call z_rot3_vec3gen(dble(D1(ii)),aimag(D1(ii)),E(ii),cr,ci,s,nrm)
     !print*, ii, cr, ci, s
     
     QB(3*ii-2) = cr
     QB(3*ii-1) = ci
     QB(3*ii)   = s
     QC(3*ii-2) = cr
     QC(3*ii-1) = ci
     QC(3*ii)   = -s
     ! update banded matrix (ignoring the second superdiagonal and the upper part of the matrix)
     ! D1(ii) = cmplx(cr,-ci,kind=8)*D1(ii) + s*EL(ii)
     D1(ii+1) = -s*EU(ii) + cmplx(cr,ci,kind=8)*D1(ii+1)
     if (ii<N-1) then
        EU(ii+1) = cmplx(cr,ci,kind=8)*EU(ii+1)
     end if
  end do

  ! form diagonal Rb (Rc)^-1 (which is a diagonal matrix) 
  do ii=1,N
     if (ii .LT. N) then 
        ! D1(ii) = cmplx(1d0,0d0,kind=8)
        DR(2*ii-1) = 1d0
        DR(2*ii)   = 0d0
     else
        D1(ii) = D1(ii)/conjg(D1(ii))
        DR(2*ii-1) = dble(D1(ii))
        DR(2*ii)   = aimag(D1(ii))
        call d_rot2_vec2gen(DR(2*ii-1),DR(2*ii),DR(2*ii-1),DR(2*ii),nrm)             
     end if
     !print*, DR(2*ii-1), DR(2*ii), DR(2*ii-1)**2 + DR(2*ii)**2
  end do

  ! fuse QC into QB
  do jj=N-1,1,-1
     bulge(1) = QC(3*jj-2)
     bulge(2) = QC(3*jj-1)
     bulge(3) = QC(3*jj)

     ! main chasing loop
     do ii=jj,(N-2)
       
        ! set indices
        ind1 = 2*(ii-1) + 1
        ind2 = ind1+3
        
        ! through diag
        call z_rot3_swapdiag(.FALSE.,DR(ind1:ind2),bulge)
        
        ! set indices
        ind1 = 3*(ii-1) + 1
        ind2 = ind1+2
        
        ! through QB
        call z_rot3_turnover(QB(ind1:ind2),QB((ind1+3):(ind2+3)),bulge)

        ! update eigenvectors
        if (VEC) then
           t1(1,1) = cmplx(bulge(1),bulge(2),kind=8)
           t1(2,1) = cmplx(bulge(3),0d0,kind=8)
           t1(1,2) = -t1(2,1)
           t1(2,2) = conjg(t1(1,1))
           Z(:,(ii+1):(ii+2)) = matmul(Z(:,(ii+1):(ii+2)),t1)
        end if
        
     end do
     
!!$     ! update eigenvectors
!!$     if (VEC) then
!!$        t1(1,1) = cmplx(bulge(1),bulge(2),kind=8)
!!$        t1(2,1) = cmplx(bulge(3),0d0,kind=8)
!!$        t1(1,2) = -t1(2,1)
!!$        t1(2,2) = conjg(t1(1,1))
!!$        Z(:,(N-1):N) = matmul(Z(:,(N-1):N),t1)
!!$     end if
     
     ! fusion at bottom
     call z_unifact_mergebulge(.FALSE.,QB((3*N-5):(3*N-3)),DR((2*N-3):(2*N)),bulge)
     !print*, jj, DR(2*N-5:2*N)

  end do
  
  ! save H
!!$  H = 0d0
!!$  do ii=1,N
!!$     H(ii,ii) = cmplx(DR(2*ii-1),DR(2*ii),kind=8)
!!$  end do
!!$  do ii=N-1,1,-1
!!$     t1(1,1) = cmplx(QB(3*ii-2),QB(3*ii-1),kind=8)
!!$     t1(2,1) = cmplx(QB(3*ii),0d0,kind=8)
!!$     t1(1,2) = -t1(2,1)
!!$     t1(2,2) = conjg(t1(1,1))
!!$     H(ii:(ii+1),:) = matmul(t1,H(ii:(ii+1),:))    
!!$  end do
    
  ! compute eigenvalues
  call z_unifact_qr(VEC,.FALSE.,N,QB,DR,M,Z,ITS,INFO)
  !call z_unifact_qr(VEC,.TRUE.,N,QB,DR,M,Z2,ITS,INFO)

!!$  H2 = 0d0
!!$  do ii=1,N
!!$     H2(ii,ii) = cmplx(DR(2*ii-1),DR(2*ii),kind=8)
!!$     print*, H2(ii,ii)
!!$  end do
!!$  do ii=1,N
!!$     do jj=1,N
!!$        Z2H(ii,jj) = conjg(Z2(jj,ii))
!!$     end do
!!$  end do
!!$
!!$  H2(1:N,1:N) = matmul(H2(1:N,1:N),Z2H(1:N,1:N))
!!$  H2(1:N,1:N) = matmul(Z2(1:N,1:N),H2(1:N,1:N))
!!$  
!!$  print*, ""
!!$  print*, "H - H2"
!!$  do ii=1,N
!!$     do jj=1,N
!!$        print*, ii,jj, H(ii,jj)-H2(ii,jj), H(ii,jj), H2(ii,jj)
!!$     end do
!!$  end do


  ! check info
  if (INFO.NE.0) then 
     ! print error in debug mode
     call u_infocode_check(__FILE__,__LINE__,"z_unifact_qr failed",INFO,INFO)
     if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_unifact_qr failed",INFO,INFO)
     end if
     INFO = 1
     !return
  end if

  ! ca = 1
  ! cb = i
  ! cg = 1
  ! cd = -i
  
  ! back transformation
  do ii=1,N
     !print*, DR(2*ii-1),DR(2*ii)
     D1(ii) = cmplx(DR(2*ii-1),DR(2*ii),kind=8)
     D1(ii) = (cd*D1(ii)-cb)/(-cg*D1(ii)+ca)
     ! print*, D1(ii)
     D(ii) = dble(D1(ii))
  end do

  
  
end subroutine d_symtrid_qr
