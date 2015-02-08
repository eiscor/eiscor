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
subroutine z_poly_roots(N,COEFFS,ROOTS)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  complex(8), intent(in) :: COEFFS(N+1)
  complex(8), intent(inout) :: ROOTS(N)
  
  ! compute variables
  integer :: ii, INFO
  logical :: QZ, VEC, ID
  logical, allocatable :: P(:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:),D(:),C(:),B(:)
  complex(8) :: V
  interface
    function l_upr1fact_upperhess(m,flags)
      logical :: l_upr1fact_upperhess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_upperhess
  end interface
  
  ! set variables
  QZ = .FALSE.
  VEC = .FALSE.
  ID = .FALSE.

  ! allocate memory
  allocate(P(N-2),ITS(N-1),Q(3*(N-1)),D(2*(N+1)),C(3*N),B(3*N))    

  ! factor companion matrix
  call z_poly_factorcomp(QZ,VEC,ID,N,COEFFS,P,Q,D,C,B,D,C,B,V,V)  
    
  ! call z_upr1fact_twistedqz
  call z_upr1fact_twistedqz(QZ,VEC,ID,l_upr1fact_upperhess,N,P,Q,D,C,B,D,C,B,V,V,ITS,INFO)
    
  ! extract roots
  call z_upr1fact_extracttri(.TRUE.,N,D,C,B,ROOTS)
    
  ! free memory
  deallocate(P,ITS,Q,D,C,B)

end subroutine z_poly_roots
