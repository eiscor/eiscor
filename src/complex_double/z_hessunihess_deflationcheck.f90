#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_hessunihess_deflationcheck 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks for deflations in a dense Hessenberg-unitary
! Hessenberg pencil. When a deflation occurs the corresponding 
! rotation in the unitary part is set to the identity matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: compute Schurvectors
!                    .FALSE.: ignore V,W, no Schurvectors
!
!  N               INTEGER
!                    dimension of matrix
!
!  STR             INTEGER
!                    index first active rotation
!
!  STP             INTEGER
!                    index last active rotation
!
!  A               COMPLEX(8) array of dimension (N,N)
!                    array of the dense Hessnberg matrix
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for unitary Hessenberg matrix
!
!  M               INTEGER
!                    leading dimension of V and W
!
!  V,W             COMPLEX(8) array of dimension (M,N)
!                    array of eigenvectors
!
! OUTPUT VARIABLES:
!
!  ZERO            INTEGER
!                     on output contains index of newest deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_hessunihess_deflationcheck(VEC,N,STR,STP,A,Q,M,V,W,ZERO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N, M, STR, STP
  logical, intent(in) :: VEC
  real(8), intent(inout) :: Q(3*(N-1))
  complex(8), intent(inout) :: A(N,N), V(M,N), W(M,N)
  integer, intent(inout) :: ZERO

  ! compute variables
  integer :: ii, down
  real(8), parameter :: tol = EISCOR_DBL_EPS
  real(8) :: qr, qi, nrm, s 
  
  ! check for deflation
  do ii=STP,STR,-1
         
    
    ! deflate if subdiagonal is small enough
    nrm = abs(Q(3*ii)) + abs(A(ii+1,ii))
    if ((ii .EQ. STR).AND.(nrm .GE. 2*tol)) then
      ! check for parallelity
      ! first column
 
    end if
    if ((ii .EQ. STR).AND.(nrm .GE. 2*tol)) then
      ! check for parallelity
      ! last row

    end if
    if (nrm .LT. 2*tol) then

      ! set ZERO
      ZERO = ii-STR+1
      
      ! extract diagonal
      qr = Q(3*ii-2)
      qi = Q(3*ii-1)
                
      ! set rotation to identity
      Q(3*ii-2) = 1d0
      Q(3*ii-1) = 0d0
      Q(3*ii) = 0d0
      
      ! scale column of A
      call z_hess_unimodscale(.FALSE.,0,ii,A,cmplx(qr,qi,kind=8))
      
      ! scale row of A
      call z_hess_unimodscale(.TRUE.,0,ii+1,A,cmplx(qr,-qi,kind=8))

      ! update left schurvectors
      if (VEC) then
        V(:,ii) = V(:,ii)*cmplx(qr,qi,kind=8)
        W(:,ii+1) = W(:,ii+1)*cmplx(qr,-qi,kind=8)
      end if

      ! exit loop 
      exit
        
    end if
    
  end do
  
end subroutine z_hessunihess_deflationcheck
