#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_urffact_eigenvectors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine diagonalizes a unitary upper hessenberg matrix that is stored as
! a product of N Givens rotations, without computing square roots.
!
! | u1       -v1 |
! | v1  conj(u1) | | u2       -v2 | 
!                  | v2  conj(u2) | | u3       -v3 | | 1   0 |
!                                   | v3  conj(u3) | | 0  u4 | 
!                                                                     
! The square root free algorithm only requires the storage of the vi^2,
! so the arrays U and VV contain the following:
!
!  U(i) = ui
! VV(i) = vi^2
!
! The input must satisfy the following:
!
!  |U(i)|^2 + VV(i) = 1
!             VV(N) = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  U               COMPLEX(8) array of dimension N
!                    array of complex generators for Givens rotations
!                    on output contains eigenvalues
!
!  VV              REAL(8) array of dimension N
!                    array of real generators (squared) for Givens rotations
!
! OUTPUT VARIABLES:
!
!  ITS             INTEGER array of dimension N-1
!                    contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO = 1 implies no convergence
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies N is invalid
!                    INFO = -2 implies U or VV is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_urffact_eigenvectors(N,U,VV,K,E,Z,CWORK,WORK)

  implicit none
  
  ! input variables
  integer, parameter :: nits = 2
  integer, intent(in) :: N, K
  complex(8), intent(in) :: E(N)
  complex(8), intent(inout) :: U(N), Z(N,K), CWORK(N,nits*K)
  real(8), intent(inout) :: VV(N), WORK(N,nits*K)
  
  ! compute variables
  integer :: ii, jj, kk, cind
  real(8) :: vvt, cc, ss, cr, ci, nrm
  complex(8) :: ut, w, p, rho, A(2,2)
  

  ! initialize Z
  Z = cmplx(0d0,0d0,kind=8)
  do kk = 1,K
    Z((N-K)+kk,kk) = cmplx(1d0,0d0,kind=8)
  end do
  
  ! loop for number of eigenvalues
  do kk = 1,K

    ! loop for inverse iterations
    do jj=1,nits
   
      ! set cind
      cind = nits*(kk-1) + jj 

      ! initialize singlestep
      rho = E(kk)
      w = -rho
      cc = 1d0
      ss = 0d0
      p = w
    
      ! core chasing loop
      do ii=0,(N-1)
    
        ! store rotation
        p = p*(cmplx(1d0,0d0,kind=8) + conjg(w)*U(ii+1))
        call z_rot3_vec3gen(dble(p),aimag(p),sqrt(VV(ii+1)),cr,ci,WORK(ii+1,cind),nrm)
        p = cmplx(cr,ci,kind=8)
        CWORK(ii+1,cind) = p
        
        ! set ut and vvt
        ut = U(ii+1)
        vvt = VV(ii+1)
    
        ! turnover
        call z_rfr3_turnover(w,cc,ss,ut,vvt,rho)
    
        ! store ut and vvt
        if ( ii > 0 ) then
          U(ii) = ut
          VV(ii) = vvt
        end if
    
      end do
!print*,""
!print*,"U and VV"
!do ii = 1,N
!print*,U(ii),VV(ii)
!end do
  
    end do

  end do

!print*,""
!print*,"CWORK and WORK"
!do ii = 1,N
!print*,CWORK(ii,1),WORK(ii,1),abs(CWORK(ii,1))**2+WORK(ii,1)**2,CWORK(ii,2),WORK(ii,2),abs(CWORK(ii,2))**2+WORK(ii,2)**2
!end do
  

  ! loop for number of eigenvalues
  do kk = 1,K

    ! loop for inverse iterations
    do jj=1,nits
   
      ! set cind
      cind = nits*K + 1 - (nits*(kk-1) + jj)

      ! update vectors
      Z(N,:) = -CWORK(N,cind)*conjg(rho)*Z(N,:)
      do ii = N-1,1,-1
    
        ! 2x2 matrix
        A(1,1) = CWORK(ii,cind) 
        A(2,2) = conjg(A(1,1))
        A(2,1) = WORK(ii,cind) 
        A(1,2) = -A(2,1)
        A(:,1) = -conjg(rho)*A(:,1)
        A(:,2) = -rho*A(:,2)
     
        ! update Z
        Z(ii:ii+1,:) = matmul(A,Z(ii:ii+1,:))
    
      end do
  
    end do

  end do

end subroutine z_urffact_eigenvectors
