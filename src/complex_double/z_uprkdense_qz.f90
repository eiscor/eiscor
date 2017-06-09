#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkdense_qz.f90
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
!  M               INTEGER
!                    leading dimension of V and W
! 
! OUTPUT VARIABLES:
!
!  EIGSA, EIGSB    COMPLEX(8) array of dimension (N)
!                    the quotient EIGSA(ii)/EIGSB(ii) is the eigenvalue
!                    like in LAPACK
!
!  V               COMPLEX(8) array of dimension (N,N)
!                    right schur vectors
!                    ignored if VEC == .FALSE.
!
!  W               COMPLEX(8) array of dimension (N,N)
!                    left schur vectors
!                    ignored if VEC == .FALSE.
!
!  S               COMPLEX(8) array of dimension (N,N)
!                    Upper triangular matrix which is the constant
!                    term of the generalized Schur form of the pencil. 
!                    ignored if VEC == .FALSE.
!
!  T               COMPLEX(8) array of dimension (N,N)
!                    Upper triangular matrix which is the linear
!                    term of the generalized Schur form of the pencil. 
!                    ignored if VEC == .FALSE.
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
subroutine z_uprkdense_qz(VEC,N,K,A,B,M,EIGSA,EIGSB,V,W,S,T,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: N, K, M 
  complex(8), intent(in) :: A(N,K), B(N,K)
  complex(8), intent(inout) :: V(M,N), W(M,N), S(N,N), T(N,N)
  integer, intent(inout) :: INFO
  complex(8), intent(inout) :: EIGSA(N), EIGSB(N)

  ! compute variables
  logical :: P(N-2)
  real(8) :: Q(3*K*(N-1)), D1(2*K*N), C1(3*N*K), B1(3*N*K)
  real(8) :: D2(2*K*N), C2(3*N*K), B2(3*N*K)
  complex(8) :: HS(N,N)
  integer :: ii, ll
  integer :: ITS(N)
  ! interface
  interface
    function l_upr1fact_hess(m,flags)
      logical :: l_upr1fact_hess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_hess
  end interface

  ! initialize info
  INFO = 0

  !print*, "z_uprkdense_qz",N, K

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
  
  call z_uprk_compress(.TRUE.,VEC,.TRUE.,N,K,A,B,P,Q,&
       &D1,C1,B1,D2,C2,B2,V,W,INFO)
  if (INFO.NE.0) then
     print*, "Info code from z_uprkdense_factor: ", INFO
     INFO = 1
  end if

  !print*, "Q",Q
  !print*, "B",A
  !print*, "A",B
  !print*, "P",P
  !print*, "D1",D1
  !print*, "C1",C1
  !print*, "B1",B1
  !print*, "D2",D2
  !print*, "C2",C2
  !print*, "B2",B2
  

  call z_uprkfpen_qz(VEC,.FALSE.,l_upr1fact_hess,N,k,&
       &P,Q,D1,C1,B1,D2,C2,B2,M,V,W,ITS,INFO)
  if (INFO.NE.0) then
     print*, "Info code from z_uprkfpen_qz: ", INFO
     INFO = 2
     return
  end if

  if (VEC) then
     call z_upr1utri_decompress(.FALSE.,N,&
          &D1(1:2*N),C1(1:(3*N)),&
          &B1(1:(3*N)),S)
     do ll = 2,k
        call z_upr1utri_decompress(.FALSE.,N,&
             &D1(((ll-1)*2*N+1):(ll*2*N)),C1(((ll-1)*3*N+1):(ll*3*N)),&
             &B1(((ll-1)*3*N+1):(ll*3*N)),HS)
        S = matmul(S,HS)

     end do
     do ii = 1, N
        EIGSA(ii) = S(ii,ii)
     end do
  else
     call z_upr1utri_decompress(.TRUE.,N,&
          &D1(1:(2*N)),C1(1:(3*N)),B1(1:(3*N)),HS)
     do ii = 1, N
        EIGSA(ii) = HS(ii,1)
     end do
     do ll = 2,k
        call z_upr1utri_decompress(.TRUE.,N,&
             &D1(((ll-1)*2*N+1):(ll*2*N)),C1(((ll-1)*3*N+1):(ll*3*N)),&
             &B1(((ll-1)*3*N+1):(ll*3*N)),HS)
        do ii = 1, N
           EIGSA(ii) = EIGSA(ii)*HS(ii,1)
        end do
     end do
  end if

  if (VEC) then
     call z_upr1utri_decompress(.FALSE.,N,&
          &D2(1:(2*N)),C2(1:(3*N)),&
          &B2(1:(3*N)),T)
     do ll = 2,k
        call z_upr1utri_decompress(.FALSE.,N,&
             &D2(((ll-1)*2*N+1):(ll*2*N)),C2(((ll-1)*3*N+1):(ll*3*N)),&
             &B2(((ll-1)*3*N+1):(ll*3*N)),HS)
        T = matmul(T,HS)
     end do
     do ii = 1, N
        EIGSB(ii) = T(ii,ii)
     end do
  else
     call z_upr1utri_decompress(.TRUE.,N,&
          &D2(1:(2*N)),C2(1:(3*N)),B2(1:(3*N)),HS)
     do ii = 1, N
        EIGSB(ii) = HS(ii,1)
     end do
     do ll = 2,k
        call z_upr1utri_decompress(.TRUE.,N,&
             &D2(((ll-1)*2*N+1):(ll*2*N)),C2(((ll-1)*3*N+1):(ll*3*N)),&
             &B2(((ll-1)*3*N+1):(ll*3*N)),HS)
        do ii = 1, N
           EIGSB(ii) = EIGSB(ii)*HS(ii,1)
        end do
     end do
  end if

  !print*, "end z_uprkdense_qz",N, K
  
end subroutine z_uprkdense_qz
