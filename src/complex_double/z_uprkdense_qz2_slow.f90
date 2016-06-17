#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkdense_qz2_slow.f90
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This functions computes the eigenvalues of a matrix polynomial 
! B_1 \lambda^d + (B_2 + A_1) \lambda^d-1 + ... + (A_(d-1) + B_d \lambda + A_d
! of degree d with B_1 lower triangular and A_d upper triangular 
!
! If not QZ then B_1 = I and B_i = 0 for i>1.
!
! After factorizing the problem in d upper Hessenberg matrices for A and
! d upper triangular matrices for B, the problem is transformed to a product 
! eigenvalue problem as a product of one Hessenberg matrix and 2*d-1 upper 
! trinagular matrices. This is done by similarity transformations.
!
! The product eigenvalue problen in Hessenberg, triagular, ..., triangular form
! is then solved with z_uprkfact_twistedqz.
!
! Finally the eigenvalues are retrieved
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  QZ              LOGICAL
!                    .TRUE.: second triangular factor is assumed nonzero
!                    .FALSE.: second triangular factor is assumed to be identity
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvector
!                    .FALSE.: no schurvectors
!
!  N               INTEGER
!                    dimension of matrix (K*d)
!
!  K               INTEGER
!                    rank, i.e., number of upper triangulars
!
!  A, B            COMPLEX(8) array of dimension (N,K)
!                    coefficients matrix polynomial
!                    the entries will be overwritten on exit
! 
! OUTPUT VARIABLES:
!
!  EIGSA, EIGSB    COMPLEX(8) array of dimension (N)
!                    the quotient EIGSA(ii)/EIGSB(ii) is the eigenvalue
!                    like in LAPACK
!
!  V               COMPLEX(8) array of dimension (N,N)
!                    right schur vectors
!
!  W               COMPLEX(8) array of dimension (N,N)
!                    left schur vectors
!
!  INFO            INTEGER
!                    TBA   
!                    INFO = 2 implies z_uprkdense_twistedqz failed
!                    INFO = 1 implies z_uprkdense_factor failed
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies N, Q, D1, C1, B1, D2, C2 or B2 is invalid
!                    INFO = -9 implies V is invalid
!                    INFO = -10 implies W is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkdense_qz2_slow(QZ,VEC,N,K,A,B,EIGSA,EIGSB,V,W,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: QZ, VEC
  integer, intent(in) :: N,K
  complex(8), intent(in) :: A(N,K), B(N,K)
  complex(8), intent(inout) :: V(N,N), W(N,N)
  integer, intent(inout) :: INFO

  ! compute variables
  logical :: P(N-2)
  real(8) :: Q(3*K*(N-1)), D1(2*K*(N+1)), C1(3*N*K), B1(3*N*K)
  real(8) :: D2(2*K*(N+1)), C2(3*N*K), B2(3*N*K)
  complex(8) :: MV(N,N), MW(N,N)
  real(8) :: bulge(3)
  complex(8) :: hc, Gf(2,2)
  complex(8) :: EIGSA(K*N), EIGSB(K*N), WORK(5*N)
  real(8) :: RWORK(2*N)
  integer :: ii, jj, ll, row, col, lwork
  integer :: ITS(N)
  ! interface
  interface
    function l_upr1fact_upperhess(m,flags)
      logical :: l_upr1fact_upperhess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_upperhess
  end interface

  ! initialize info
  INFO = 0
  
!!$  ! check N
!!$  if (N < 2) then
!!$    INFO = -1
!!$    ! print error in debug mode
!!$    if (DEBUG) then
!!$      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
!!$    end if
!!$    return
!!$  end if
!!$  
!!$  ! check H
!!$  call z_2Darray_check(N,N,H,flg)
!!$  if (.NOT.flg) then
!!$    INFO = -2
!!$    ! print error in debug mode
!!$    if (DEBUG) then
!!$      call u_infocode_check(__FILE__,__LINE__,"H is invalid",INFO,INFO)
!!$    end if
!!$    return
!!$  end if  
  
  call z_uprkdense_factor2_slow(QZ,VEC,.TRUE.,N,k,A,B,P,Q,&
       &D1,C1,B1,D2,C2,B2,V,W,INFO)
  if (INFO.NE.0) then
     print*, "Info code from z_uprkdense_factor: ", INFO
     INFO = 1
  end if


  call z_uprkfact_twistedqz(QZ,VEC,.FALSE.,l_upr1fact_upperhess,N,k,&
       &P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
  if (INFO.NE.0) then
     print*, "Info code from z_uprkfact_twistedqz: ", INFO
     INFO = 2
  end if

  do ll = 1,k
     call z_upr1fact_extracttri(.TRUE.,N,&
          &D1(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C1(((ll-1)*3*N+1):(ll*3*N)),&
          &B1(((ll-1)*3*N+1):(ll*3*N)),V)
     EIGSA(((ll-1)*N+1):ll*N) = V(1:N,1)
  end do
  do ll = 2,k
     do ii=1,N
        EIGSA(ii) = EIGSA(ii)*EIGSA((ll-1)*N+ii)     
     end do
  end do


  if (QZ) then
     do ll = 1,k
        call z_upr1fact_extracttri(.TRUE.,N,&
             &D2(((ll-1)*2*(N+1)+1):(ll*2*(N+1))),C2(((ll-1)*3*N+1):(ll*3*N)),&
             &B2(((ll-1)*3*N+1):(ll*3*N)),V)
        EIGSB(((ll-1)*N+1):ll*N) = V(1:N,1)
     end do
     do ll = 2,k
        do ii=1,N
           EIGSB(ii) = EIGSB(ii)*EIGSB((ll-1)*N+ii)     
        end do
     end do
  end if
  
end subroutine z_uprkdense_qz2_slow
