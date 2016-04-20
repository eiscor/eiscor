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
subroutine z_poly_roots(N,COEFFS,ROOTS,RESIDUALS)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  complex(8), intent(in) :: COEFFS(N+1)
  complex(8), intent(inout) :: ROOTS(N)
  real(8), intent(inout) :: RESIDUALS(N)
  
  ! compute variables
  integer :: ii, jj, INFO
  real(8), parameter :: pi = 3.14159265358979323846264338327950d0
  real(8) :: scl
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

  ! fill P
  P = .FALSE.

  ! set V
  do ii = 1,N
    V(ii) = -COEFFS(N+2-ii)/COEFFS(1)
  end do

  ! factor companion matrix
  call z_compmat_compress(N,P,V,Q,D,C,B,INFO)       
    
  ! call z_upr1fact_twistedqz
  call z_upr1fact_qr(.FALSE.,.FALSE.,l_upr1fact_hess, & 
                     N,P,Q,D,C,B,N,V,ITS,INFO)

  ! check INFO
  if (INFO.NE.0) then
    print*,""
    print*,"INFO:",INFO
    print*,""
    deallocate(P,ITS,Q,D,C,B,V)
    return  
  end if

  ! extract roots
  call z_upr1utri_decompress(.TRUE.,N,D,C,B,ROOTS)
    
  ! compute residuals
  call z_poly_residuals(N,COEFFS,ROOTS,0,RESIDUALS)
    
print*,""
print*,""
print*,"Inside Roots"
print*,""
print*,"INFO:",INFO
print*,""
print*,"Roots residuals and iterations"
call z_scalar_argument(dble(ROOTS(1)),aimag(ROOTS(1)),scl,INFO)
jj = nint(N*scl/pi/2d0)
print*,jj,ROOTS(1),RESIDUALS(1)
do ii=1,(N-1)
call z_scalar_argument(dble(ROOTS(ii+1)),aimag(ROOTS(ii+1)),scl,INFO)
jj = nint(N*scl/pi/2d0)
print*,jj,ROOTS(ii+1),RESIDUALS(ii+1),ITS(ii)
end do
print*,""

  ! free memory
    deallocate(P,ITS,Q,D,C,B,V)

end subroutine z_poly_roots
