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
! This routine performs a Cayley transformation of the symmetric
! tridiagonal matrix to a unitary matrix (a descending and a 
! ascending sequence of core transformations). The unitary matrix
! is transformed to upper Hessenberg form by core chasing and then 
! passed to z_unifact_qr.
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
  ! WORK(1:3N-3) = QB
  ! WORK(3N+1:5N) = DR
  integer, intent(inout) :: ITS(N-1), INFO
  real(8), intent(inout) :: D(N), E(N)
  complex(8), intent(inout) :: Z(M,N)
  ! compute variables
  integer :: ii, jj, ind1, ind2
  real(8) :: cr,ci,s,nrm
  
  complex(8) :: eu, d1
  real(8) :: bulge(3)
   
  complex(8) :: block(2,2), t1(2,2), t2(2,2)
  character(len=1024) :: filename
  
  ! initialize INFO
  INFO = 0

  if (VEC.AND.ID) then
    Z = cmplx(0d0,0d0,kind=8)
    do ii=1,min(M,N)
      Z(ii,ii) = cmplx(1d0,0d0,kind=8)
   end do
  end if
  

  ! Cayley transform -(T-iI)(T+iI)              
  ! implicitly form B = T - i I  and C = T + i I = conjg(B)
  eu = E(1)
  d1 = cmplx(D(1),-1d0,kind=8)

  ! QR decomposition B = QB Rb
  do ii=1,N-1
     ! compute new rotation
     !call z_rot3_vec4gen(dble(D1(ii)),aimag(D1(ii)),dble(EL(ii)),aimag(EL(ii)),cr,ci,s,nrm)
     !call z_rot3_vec3gen(dble(D1(ii)),aimag(D1(ii)),E(ii),cr,ci,s,nrm)
     call z_rot3_vec3gen(dble(d1),aimag(d1),E(ii),cr,ci,s,nrm)
     
     WORK(3*ii-2) = cr
     WORK(3*ii-1) = ci
     WORK(3*ii)   = s
     ! update banded matrix (ignoring the second superdiagonal and the upper part of the matrix)
     ! D1(ii) = cmplx(cr,-ci,kind=8)*D1(ii) + s*EL(ii)
     d1 = -s*eu + cmplx(cr,ci,kind=8)*cmplx(D(ii+1),-1d0,kind=8)
     if (ii<N-1) then
        eu = cmplx(cr*E(ii+1),ci*E(ii+1),kind=8)
     end if
  end do

  ! form diagonal Rb (Rc)^-1 (which is a diagonal matrix) 
  do ii=1,N-1
     WORK(3*N + 2*ii-1) = -1d0
     WORK(3*N + 2*ii)   = 0d0
  end do
  ! ii = N
  d1 = -d1/conjg(d1)
  call d_rot2_vec2gen(dble(d1),aimag(d1),WORK(3*N + 2*N-1),WORK(3*N + 2*N),nrm)             

  ! fuse QC into QB
  do jj=N-1,1,-1
     bulge(1) = WORK(3*jj-2)
     bulge(2) = WORK(3*jj-1)
     bulge(3) = -WORK(3*jj)  ! QC is linked to QB

     ! main chasing loop
     do ii=jj,(N-2)
       
        ! set indices
        ind1 = 3*N + 2*(ii-1) + 1  ! DR in WORK
        ind2 = 3*N + ind1+3
        
        ! through diag 
        call z_rot3_swapdiag(.FALSE.,WORK(ind1:ind2),bulge)
        
        ! set indices
        ind1 = 3*(ii-1) + 1
        ind2 = ind1+2
        
        ! through QB
        call z_rot3_turnover(WORK(ind1:ind2),WORK((ind1+3):(ind2+3)),bulge)

        ! update eigenvectors
        if (VEC) then
           t1(1,1) = cmplx(bulge(1),bulge(2),kind=8)
           t1(2,1) = cmplx(bulge(3),0d0,kind=8)
           t1(1,2) = -t1(2,1)
           t1(2,2) = conjg(t1(1,1))
           Z(:,(ii+1):(ii+2)) = matmul(Z(:,(ii+1):(ii+2)),t1)
        end if

     end do
          
     ! set indices
     ind1 = 3*N + 2*(N-2) + 1  ! DR in WORK
     ind2 = 3*N + ind1+3
     
     ! through diag
     call z_rot3_swapdiag(.FALSE.,WORK(ind1:ind2),bulge)
     
     ! fusion at bottom
     call z_unifact_mergebulge(.FALSE.,WORK((3*N-5):(3*N-3)),WORK((3*N + 2*N-3):(3*N + 2*N)),bulge)

  end do

!!$  do ii=1,N-1
!!$     print*, ii, WORK(3*ii-2), WORK(3*ii-1), WORK(3*ii)
!!$     print*, ii, WORK(3*N+2*ii-1), WORK(3*N+2*ii)
!!$  end do

!!$  ! write data to file
!!$  write (filename, "(A2,I4,A4)") "QD",N,".txt"
!!$  open (unit=7, file=filename, status='unknown', position="rewind")
!!$  write (7,*) "new file"
!!$  do ii=1,N-1
!!$     write (7,*) WORK(3*ii-2), WORK(3*ii-1), WORK(3*ii)
!!$  end do
!!$  do ii=1,N
!!$     write (7,*) WORK(3*N+2*ii-1), WORK(3*N+2*ii)
!!$  end do
!!$  close(7)
  
  ! compute eigenvalues
  call z_unifact_qr(VEC,.FALSE.,N,WORK(1:3*N-3),WORK((3*N+1):(5*N-1)),M,Z,ITS,INFO)

  ! check info
  if (INFO.NE.0) then 
     ! print error in debug mode
     if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_unifact_qr failed",INFO,INFO)
     end if
     INFO = 1
     ! since some of the eigenvalues have been found, the back transform is performed for all
     ! return
  end if
  
  ! back transformation
  do ii=1,N
     D(ii) = WORK(3*N+2*ii)/(1d0+WORK(3*N+2*ii-1))
     !print*, WORK(3*N+2*ii-1), WORK(3*N+2*ii), D(ii)
  end do
  
end subroutine d_symtrid_qr
