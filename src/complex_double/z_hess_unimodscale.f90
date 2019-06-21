#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_hess_unimodscale
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine scales either a row or column of a unitary plus rank
! one upper-triangular matrix (upr1utri) by a complex unimodular scalar.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  ROW             LOGICAL
!                    .TRUE.: scale row
!                    .FALSE.: scale column
!
!  N               INTEGER
!                    dimension of matrix
!
!  N               INTEGER
!                    index of row/column to be scaled
!
!  A               COMPLEX(8) array of dimension (N,N)
!                    array of the dense Hessenberg matrix
!
!  SCL             COMPLEX(8) 
!                    scalar, assumed unimodular
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_hess_unimodscale(ROW,N,ii,A,SCL)

  implicit none
  
  ! input variables
  logical, intent(in) :: ROW
  integer, intent(in) :: N, ii
  complex(8), intent(inout) :: A(N,N)
  complex(8), intent(in) :: SCL
 
  ! compute variables
  real(8) :: nrm
  complex(8) :: temp
 
  ! update A
  if (.NOT.ROW) then

    if (ii .LT. N) then
      A(1:ii+1,ii) = A(1:ii+1,ii) * SCL
    else
      A(1:N,ii) = A(1:N,ii) * SCL
    end if

  else

    if (ii .GT. 1) then
      A(ii,ii-1:N) = A(ii,ii-1:N) * SCL
    else
      A(ii,1:N) = A(ii,1:N) * SCL
    end if
    
  end if

end subroutine z_hess_unimodscale
