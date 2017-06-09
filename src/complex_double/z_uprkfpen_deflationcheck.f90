#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkfpen_deflationcheck 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks for deflations in a factored unitary plus rank
! one (uprkfpen) matrix. When a deflation occurs the corresponding 
! rotation in the unitary part is set to the identity matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VECL            LOGICAL
!                    .TRUE.: compute schurvector
!                    .FALSE.: no schurvectors
!                    right Schur vectors are not updated
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D1,D2           REAL(8) arrays of dimension (2*N)
!                    arrays of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1,C2,B1,B2     REAL(8) arrays of dimension (3*N)
!                    arrays of generators for unitary plus rank one
!                    upper-trinagular matrices
!
!  M               INTEGER
!                    leading dimension of V and W
!
!  V,W             COMPLEX(8) array of dimension (M,N)
!                    right and left schurvectors 
!
! OUTPUT VARIABLES:
!
!  ZERO            INTEGER
!                     on output contains index of newest deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkfpen_deflationcheck(VECL,N,K,STR,STP,P,Q,D1,C1,B1,D2,C2,B2,M,V,W,ZERO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N, M, K, STR, STP
  logical, intent(in) :: VECL, P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*N*K), C1(3*N*K), B1(3*N*K)
  real(8), intent(inout) :: D2(2*N*K), C2(3*N*K), B2(3*N*K)
  complex(8), intent(inout) :: V(M,N), W(M,N)
  integer, intent(inout) :: ZERO

  ! compute variables
  integer :: ii, jj, up, down
  real(8), parameter :: tol = EISCOR_DBL_EPS
  real(8) :: qr, qi, dr, di, cr, ci, s, nrm
  
  ZERO = 0
  ! check for deflation
  do ii=(N-STP),(N-STR)
  
    ! deflate if subdiagonal is small enough
    nrm = abs(Q(3*(N-ii)))
    print*, "deflation", N-ii
    if(nrm < tol)then

      ! set ZERO
      ZERO = N-ii
      
      ! extract diagonal
      qr = Q(3*ZERO-2)
      qi = Q(3*ZERO-1)
                
      ! set rotation to identity
      Q(3*ZERO-2) = 1d0
      Q(3*ZERO-1) = 0d0
      Q(3*ZERO) = 0d0
      
      ! upper left entry of Q(ZERO)
      ! deflation at top
      if ( ZERO.EQ.STR ) then

        ! scale row of R1
        call z_upr1utri_unimodscale(.TRUE.,D1(2*ZERO-1:2*ZERO), &
             C1(3*ZERO-2:3*ZERO),B1(3*ZERO-2:3*ZERO),cmplx(qr,qi,kind=8))
        ! .TRUE. means C1, B1 don't matter

      ! deflation in the middle
      ! P(ZERO-1) == 0
      elseif ( .NOT.P(ZERO-1) ) then

        ! scale row of R1
        call z_upr1utri_unimodscale(.TRUE.,D1(2*ZERO-1:2*ZERO), &
             C1(3*ZERO-2:3*ZERO),B1(3*ZERO-2:3*ZERO),cmplx(qr,qi,kind=8))
        ! .TRUE. means C1, B1 don't matter

      ! P(ZERO-1) == 1
      else

        ! scale row of R2
        call z_upr1utri_unimodscale(.TRUE.,D2(2*ZERO-1:2*ZERO), &
             C2(3*ZERO-2:3*ZERO),B2(3*ZERO-2:3*ZERO),cmplx(qr,-qi,kind=8))
        ! .TRUE. means C2, B2 don't matter

        ! update left schurvectors
        if (VECL) then
          W(:,ZERO) = W(:,ZERO)*cmplx(qr,qi,kind=8)
        end if

      end if
      
      ! lower entry of Q(ZERO)
      ! deflation at bottom
      if ( ZERO.EQ.STP ) then

        ! scale row of R1
        call z_upr1utri_unimodscale(.TRUE.,D1(2*ZERO+1:2*ZERO+2), &
             C1(3*ZERO+1:3*ZERO+3),B1(3*ZERO+1:3*ZERO+3),cmplx(qr,-qi,kind=8))
        ! .TRUE. means C1, B1 don't matter

      ! deflation in the middle
      ! P(ZERO) == 0
      elseif ( .NOT.P(ZERO) ) then

        ! scale row of R2
        call z_upr1utri_unimodscale(.TRUE.,D2(2*ZERO+1:2*ZERO+2), &
             C2(3*ZERO+1:3*ZERO+3),B2(3*ZERO+1:3*ZERO+3),cmplx(qr,qi,kind=8))
        ! .TRUE. means C2, B2 don't matter

        ! update left schurvectors
        if (VECL) then
          W(:,ZERO+1) = W(:,ZERO+1)*cmplx(qr,-qi,kind=8)
        end if

      ! P(ZERO) == 1
      else

        ! scale row of R1
        call z_upr1utri_unimodscale(.TRUE.,D1(2*ZERO+1:2*ZERO+2), &
             C1(3*ZERO+1:3*ZERO+3),B1(3*ZERO+1:3*ZERO+3),cmplx(qr,-qi,kind=8))
        ! .TRUE. means C1, B1 don't matter

      end if

      ZERO = ZERO - STR +1
      ! exit loop 
      exit
        
    end if
    
  end do

  print*, "deflation ZERO", ZERO
  
end subroutine z_uprkfpen_deflationcheck
