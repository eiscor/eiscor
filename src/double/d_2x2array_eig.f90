#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_2x2array_eig 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Schur decomposition of a general
! 2x2 real matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  H              REAL(8) array of dimension (2,2)
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
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies H is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_2x2array_eig(H,E,Z,INFO)
  
  implicit none
  
  ! input variables
  integer, intent(inout) :: INFO
  real(8), intent(in) :: H(2,2)
  complex(8), intent(inout) :: E(2), Z(2,2)
  
  ! compute variables
  integer :: ii, id
  real(8) :: trace, disc, detm, temp, nrm1, nrm2
  complex(8) :: WORK(4), swap
  
  ! initialize info
  INFO = 0
 
  ! print error in debug mode
  if (DEBUG) then
  
    ! check H
    call d_2Darray_check(2,2,H,INFO)
    call u_infocode_check(__FILE__,__LINE__,"H is invalid",INFO,-1)
    if (INFO.NE.0) then
      return
    end if 

  end if  

  ! compute intermediate values
  trace = H(1,1) + H(2,2)
  detm = H(1,1)*H(2,2) - H(2,1)*H(1,2)
  temp = (H(1,1)-H(2,2))**2 + 4d0*H(1,2)*H(2,1)
  
  ! zero roots
  if (abs(detm).EQ.0d0) then
    E(1) = cmplx(0d0,0d0,kind=8)
    E(2) = cmplx(0d0,0d0,kind=8)   
  
  ! imaginary roots
  else if (temp < 0) then
    disc = sqrt(-temp)
    E(1) = cmplx(trace,disc,kind=8)/2d0
    E(2) = cmplx(trace,-disc,kind=8)/2d0    

  ! real roots
  else
    disc = sqrt(temp)
    
    ! compute E
    if(abs(trace+disc) > abs(trace-disc))then
      E(1) = cmplx((trace+disc)/2d0,0d0,kind=8)
      E(2) = cmplx(detm,0d0,kind=8)/E(1)
    else
      E(1) = cmplx((trace-disc)/2d0,0d0,kind=8)
      E(2) = cmplx(detm,0d0,kind=8)/E(1)
    end if
    
  end if
  
  ! compute diagonal with least cancellation
  WORK(1) = H(1,1) - E(1)
  WORK(2) = H(2,2) - E(1)
  WORK(3) = H(1,1) - E(2)
  WORK(4) = H(2,2) - E(2)
  
  id = 1
  ! complex abs does not matter here
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
    swap = E(1)
    E(1) = E(2)
    E(2) = swap
  else
    Z(1,1) = -WORK(4)
    Z(2,1) = H(2,1)
    swap = E(1)
    E(1) = E(2)
    E(2) = swap
  end if

  ! normalize first column of Z
  nrm1 = abs(Z(1,1)) 
  nrm2 = abs(Z(2,1))
  ! complex abs does not matter here, since nrm2 < eps_single 
  ! => nrm2**2 < eps_double
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
  
end subroutine d_2x2array_eig
