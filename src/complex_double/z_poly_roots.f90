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
  integer :: ii, INFO
  real(8) :: scl
  logical, allocatable :: P(:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:),D1(:),C1(:),B1(:)
  real(8), allocatable :: D2(:),C2(:),B2(:)
  complex(8), allocatable :: V(:),W(:),T(:,:)
  interface
    function l_upr1fact_hess(m,flags)
      logical :: l_upr1fact_hess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_hess
  end interface
  interface
    function l_upr1fact_inversehess(m,flags)
      logical :: l_upr1fact_inversehess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_inversehess
  end interface
  interface
    function l_upr1fact_cmv(m,flags)
      logical :: l_upr1fact_cmv
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_cmv
  end interface
  interface
    function l_upr1fact_random(m,flags)
      logical :: l_upr1fact_random
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_random
  end interface
  
  ! allocate memory
  allocate(P(N-2),ITS(N-1),Q(3*(N-1)),D1(2*(N+1)),C1(3*N),B1(3*N))    
  allocate(V(N),W(N),D2(2*(N+1)),C2(3*N),B2(3*N),T(N,N))    

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

  ! factor companion matrix
  call z_comppenc_factor(.TRUE.,N,P,V,W,Q,D1,C1,B1,D2,C2,B2,INFO)  
    
  ! call z_upr1fact_twistedqz
!  call z_upr1fact_twistedqz(.TRUE.,.FALSE.,.FALSE.,l_upr1fact_hess,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
!  call z_upr1fact_twistedqz(.TRUE.,.FALSE.,.FALSE.,l_upr1fact_inversehess,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
!  call z_upr1fact_twistedqz(.TRUE.,.FALSE.,.FALSE.,l_upr1fact_cmv,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
  call z_upr1fact_twistedqz(.TRUE.,.FALSE.,.FALSE.,l_upr1fact_random,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
!  call z_upr1fact_twistedqz(.FALSE.,.FALSE.,.FALSE.,l_upr1fact_hess,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)

  ! check INFO
  if (INFO.NE.0) then
    print*,""
    print*,"INFO:",INFO
    print*,""
    deallocate(P,ITS,Q,D1,C1,B1,D2,C2,B2,V,W,T)
    return  
  end if

  ! extract roots
  call z_upr1fact_extracttri(.TRUE.,N,D1,C1,B1,V)
  call z_upr1fact_extracttri(.TRUE.,N,D2,C2,B2,W)
  do ii=1,N
    ROOTS(ii) = V(ii)/W(ii)
!    ROOTS(ii) = T(ii,ii)
  end do
    
  ! compute residuals
  call z_poly_residuals(N,COEFFS,ROOTS,0,RESIDUALS)
    
print*,""
print*,""
print*,"Inside Roots"
print*,""
print*,"INFO:",INFO
print*,""
print*,"Roots residuals and iterations"
print*,ROOTS(1),RESIDUALS(1)
do ii=1,(N-1)
print*,ROOTS(ii+1),RESIDUALS(ii+1),ITS(ii)
end do
print*,""

  ! free memory
  deallocate(P,ITS,Q,D1,C1,B1,D2,C2,B2,V,W,T)

end subroutine z_poly_roots
