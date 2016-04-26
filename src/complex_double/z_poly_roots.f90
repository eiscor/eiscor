#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_poly_roots 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the roots of a polynomial expressed in the 
! monomial basis using the fast algorithm described in:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    degree of the polynomial
!
!  COEFFS          COMPLEX(8) array of dimension (N+1)
!                    coefficients of polynomial ordered from highest
!                    degree coefficient to lowest degree
!
! OUTPUT VARIABLES:
!
!  ROOTS           COMPLEX(8) array of dimension (N)
!                    computed roots
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_poly_roots(N,COEFFS,ROOTS,RESIDUALS,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  complex(8), intent(in) :: COEFFS(N+1)
  complex(8), intent(inout) :: ROOTS(N)
  real(8), intent(inout) :: RESIDUALS(N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii, jj
  logical, allocatable :: P(:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:), D(:), C(:), B(:)
  complex(8), allocatable :: V(:)
  interface
    function l_upr1fact_hess(m,flags)
      logical :: l_upr1fact_hess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_hess
  end interface
  
  ! allocate memory
  allocate(P(N-2),ITS(N-1),Q(3*(N-1)),D(2*N),C(3*N),B(3*N),V(N))    

  ! initialize INFO
  INFO = 0

  ! fill P
  P = .FALSE.

  ! set V
  V(N) = ((-1d0)**(N))*COEFFS(N+1)/COEFFS(1)
  do ii=1,(N-1)
    V(ii) = -COEFFS(N+1-ii)/COEFFS(1)
  end do

  ! factor companion matrix
  call z_compmat_compress(N,P,V,Q,D,C,B)       

  ! call z_upr1fact_qr
  call z_upr1fact_qr(.FALSE.,.FALSE.,l_upr1fact_hess, & 
                     N,P,Q,D,C,B,N,V,ITS,INFO)

  ! check INFO
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__, & 
           "z_upr1fact_qr failed to compute eigenvalues",INFO,INFO)
    end if
    INFO = 1
    deallocate(P,ITS,Q,D,C,B,V)
    return
  end if

  ! extract roots
  call z_upr1utri_decompress(.TRUE.,N,D,C,B,ROOTS)
    
  ! compute residuals
  call z_poly_residuals(N,COEFFS,ROOTS,0,RESIDUALS)

  ! free memory
  deallocate(P,ITS,Q,D,C,B,V)

end subroutine z_poly_roots
