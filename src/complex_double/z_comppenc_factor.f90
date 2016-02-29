#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_comppenc_factor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  QZ              LOGICAL
!                    .TRUE.: second triangular factor is used
!                    .FALSE.: second triangular factor is assumed 
!                             to be identity
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  V               COMPLEX(8) array of dimension (N)
!                    coefficients for left triangular factor
!
!  W               COMPLEX(8) array of dimension (N)
!                    coefficients for right traingular factor
!
! OUTPUT VARIABLES:
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D1,D2           REAL(8) arrays of dimension (2*(N+1))
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1,B1,C2,B2     REAL(8) arrays of dimension (3*N)
!                    array of generators for upper-triangular parts
!                    of the pencil
!
!  INFO           INTEGER
!                   INFO = 1 implies no convergence 
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_comppenc_factor(QZ,N,P,V,W,Q,D1,C1,B1,D2,C2,B2,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: QZ
  integer, intent(in) :: N
  logical, intent(in) :: P(N-2)
  complex(8), intent(inout) :: V(N), W(N)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*(N+1)), C1(3*N), B1(3*N)
  real(8), intent(inout) :: D2(2*(N+1)), C2(3*N), B2(3*N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii
  !real(8) :: phr, phi, nrm, beta
  complex(8) :: temp
  
  ! initialize info
  INFO = 0
  
  ! initialize Q
  Q = 0d0
  do ii = 1,(N-1)
    Q(3*ii) = 1d0
  end do

  ! shuffle V
  do ii = 1,(N-2)

    ! permute if necessary
    if (P(N-ii-1)) then

      temp = (-1d0)**(N-ii-1)*V(N-ii)
      V(2:(N-ii)) = V(1:(N-ii-1))
      V(N-ii) = temp    

    end if
 
  end do

  call z_upr1_factoridpspike(.TRUE.,N,V,D1,C1,B1,INFO)

  ! shuffle W
  if (QZ) then

    do ii = 1,(N-2)

      ! permute if necessary
      if (P(N-ii-1)) then

        temp = (-1d0)**(N-ii-1)*W(N-ii)
        W(2:(N-ii)) = W(1:(N-ii-1))
        W(N-ii) = temp    

      end if
 
    end do

    call z_upr1_factoridpspike(.TRUE.,N,W,D2,C2,B2,INFO)

  end if

end subroutine z_comppenc_factor
