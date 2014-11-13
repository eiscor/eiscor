#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZTTEEV (Zomplex Two by Two Eigenvalues and EigenVectors)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Schur decomposition of a general
! 2x2 COMPLEX(kind=8) matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  H              COMPLEX(8) array of dimension (2,2)
!                   Contains the 2x2 matrix.
!
! OUTPUT VARIABLES:
!
!  E              COMPLEX(8) array of dimension (2)
!                   On exit contains the eigenvalues of H.
!
!  Z              COMPLEX(8) array of dimension (2,2)
!                   On exit the columns of Z contain the Schur vectors 
!                   of H sorted to match the array E.
!
! INFO            INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies H is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZTTEEV(H,E,Z,INFO)
  
  implicit none
  
  ! input variables
  integer, intent(inout) :: INFO
  complex(8), intent(in) :: H(2,2)
  complex(8), intent(inout) :: E(2), Z(2,2)
  
  ! compute variables
  integer :: ii, id
  real(8) :: temp, nrm1, nrm2
  complex(8) :: WORK(4)
  complex(8) :: trace, detm, disc
  
  ! initialize info
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
  
    ! check H
    call ZARACH2(2,2,H,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"H is invalid",INFO,-1)
      return
    end if
  
  end if
  
  ! compute intermediate values
  trace = H(1,1) + H(2,2)
  detm = H(1,1)*H(2,2) - H(2,1)*H(1,2)
  disc = sqrt((H(1,1)-H(2,2))**2 + cmplx(4d0,0d0,kind=8)*H(1,2)*H(2,1))
  
  ! compute E
  if(abs(trace+disc) > abs(trace-disc))then
     if(abs(trace+disc) == 0)then
        E(1) = cmplx(0d0,0d0,kind=8)
        E(2) = cmplx(0d0,0d0,kind=8)
     else
        E(1) = (trace+disc)/cmplx(2d0,0d0,kind=8)
        E(2) = detm/E(1)
     end if
  else
     if(abs(trace-disc) == 0)then
        E(1) = cmplx(0d0,0d0,kind=8)
        E(2) = cmplx(0d0,0d0,kind=8)
     else
        E(1) = (trace-disc)/cmplx(2d0,0d0,kind=8)
        E(2) = detm/E(1)
     end if
  end if
  
  ! compute diagonal with least cancellation
  WORK(1) = H(1,1) - E(1)
  WORK(2) = H(2,2) - E(1)
  WORK(3) = H(1,1) - E(2)
  WORK(4) = H(2,2) - E(2)
  
  id = 1
  temp = abs(WORK(1))
  do ii=1,3
    if(abs(WORK(ii+1)) > temp) then
        id = ii+1
        temp = abs(WORK(id))
    end if
  end do
  
  ! initialize first column of Z
  ! swap eigenvalues if necessary
  if (temp .EQ. 0d0) then
    Z(1,1) = cmplx(1d0,0d0,kind=8)
    Z(2,1) = cmplx(0d0,0d0,kind=8)
  else if (id .EQ. 1) then
    Z(1,1) = H(1,2)
    Z(2,1) = -WORK(1)
  else if (id .EQ. 2) then
    Z(1,1) = -WORK(2)
    Z(2,1) = H(2,1)
  else if (id .EQ. 3) then
    Z(1,1) = H(1,2)
    Z(2,1) = -WORK(3)
    trace = E(1)
    E(1) = E(2)
    E(2) = trace
  else
    Z(1,1) = -WORK(4)
    Z(2,1) = H(2,1)
    trace = E(1)
    E(1) = E(2)
    E(2) = trace
  end if
  
  ! normalize first column of Z
  nrm1 = abs(Z(1,1))
  nrm2 = abs(Z(2,1))
  if (nrm1 > nrm2) then
    nrm2 = nrm2/nrm1
    temp = nrm1*sqrt(1d0 + nrm2*nrm2)
  else
    nrm1 = nrm1/nrm2
    temp = nrm2*sqrt(1d0 + nrm1*nrm1)
  end if
  Z(1,1) = Z(1,1)/temp
  Z(2,1) = Z(2,1)/temp
  
  ! store second column of Z
  Z(1,2) = -conjg(Z(2,1))
  Z(2,2) = conjg(Z(1,1))
  
end subroutine ZTTEEV
