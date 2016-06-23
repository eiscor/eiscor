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
  real(8) :: scl
  logical, allocatable :: P(:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:), D1(:), C1(:), B1(:)
  real(8), allocatable :: D2(:), C2(:), B2(:)
  complex(8), allocatable :: V(:), W(:)
  complex(8), allocatable :: H(:,:), T(:,:)
  interface
    function l_upr1fact_hess(m,flags)
      logical :: l_upr1fact_hess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_hess
  end interface
  
  ! allocate memory
  allocate(P(N-2),ITS(N-1),Q(3*(N-1)),D1(2*N),C1(3*N),B1(3*N))    
  allocate(D2(2*N),C2(3*N),B2(3*N),V(N),W(N))    
  allocate(H(N,N),T(N,N))    

  ! initialize INFO
  INFO = 0

  ! fill P
  P = .FALSE.

  ! fill V and W
  scl = maxval(abs(COEFFS))
  V(N) = ((-1d0)**(N))*COEFFS(N+1)/scl
  do ii=1,(N-1)
    V(ii) = -COEFFS(N+1-ii)/scl
  end do
  W = cmplx(0d0,0d0,kind=8)
  W(N) = COEFFS(1)/scl

  ! compress companion pencil
  call z_comppen_compress(N,P,V,W,Q,D1,C1,B1,D2,C2,B2) 

  ! decompress 
  call z_upr1fpen_decompress(N,P,Q,D1,C1,B1,D2,C2,B2,H,T)

  ! print pencil
print*,""
print*,"H"
do ii=1,N
write(*,2000) H(ii,:)
end do
print*,""
print*,"T"
do ii=1,N
write(*,2000) T(ii,:)
end do
print*,""

  ! free memory
  deallocate(P,Q,D1,C1,B1,D2,C2,B2,V,W,H,T,ITS)

!  ! set V
!  V(N) = ((-1d0)**(N))*COEFFS(N+1)/COEFFS(1)
!  do ii=1,(N-1)
!    V(ii) = -COEFFS(N+1-ii)/COEFFS(1)
!  end do
!
!  ! factor companion matrix
!  call z_compmat_compress(N,P,V,Q,D,C,B)       
!
!  ! call z_upr1fact_qr
!  call z_upr1fact_qr(.FALSE.,.FALSE.,l_upr1fact_hess, & 
!                     N,P,Q,D,C,B,N,V,ITS,INFO)
!
!  ! check INFO
!  if (INFO.NE.0) then
!    ! print error message in debug mode
!    if (DEBUG) then
!      call u_infocode_check(__FILE__,__LINE__, & 
!           "z_upr1fact_qr failed to compute eigenvalues",INFO,INFO)
!    end if
!    INFO = 1
!    deallocate(P,ITS,Q,D,C,B,V)
!    return
!  end if
!
!  ! extract roots
!  call z_upr1utri_decompress(.TRUE.,N,D,C,B,ROOTS)
!    
!  ! compute residuals
!  call z_poly_residuals(N,COEFFS,ROOTS,0,RESIDUALS)
!
!  ! free memory
!  deallocate(P,ITS,Q,D,C,B,V)

2000 format(9(es11.3,es10.3))

end subroutine z_poly_roots
