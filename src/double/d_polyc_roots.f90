#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_polyc_roots 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the roots of a polynomial expressed in the 
! Chebyshev basis using a Cayley transformation and the upr1
! QR algorithm. The highest degree coefficient is 1d0.
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
subroutine d_polyc_roots(N,COEFFS,ROOTS,SCA)

  implicit none
  
  ! input variables
  logical, intent(in) :: SCA
  integer, intent(in) :: N
  real(8), intent(in) :: COEFFS(N)
  complex(8), intent(inout) :: ROOTS(N)
!  real(8), intent(inout) :: RESIDUALS(N)
  
  ! compute variables
  integer :: ii, INFO
  logical, allocatable :: P(:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:),D1(:),C1(:),B1(:),D(:),E(:),Z(:,:)
  real(8), allocatable :: DWORK(:), D2(:),C2(:),B2(:)
  real(8) :: a,b,scale, norm
  complex(8), allocatable :: U(:)
  complex(8) :: V,W

!!$  complex(8) :: WORK(5*N),t(2,2)
!!$  real(8) :: RWORK(2*N)
!!$  complex(8), allocatable :: H(:,:)

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
  allocate(D(N),E(N-1),U(N),Z(1,N),DWORK(2*N))
  allocate(D2(2*(N+1)),C2(3*N),B2(3*N))    
  !allocate(H(N,N))

  ! fill P
  P = .FALSE.

  ! set up D, E, and U
  D = 0d0
  E = 5d-1
  E(1) = sqrt(2d0)/2d0
  U(1) = cmplx(-COEFFS(N)*sqrt(2d0)/2d0,0d0,kind=8)
  do ii=2,N
     U(ii) = cmplx(-COEFFS(N-ii+1)/2d0,0d0,kind=8)
  end do
  
  ! factorize \Phi(T + Ue^H) and reduce it to Hessenberg form
  !call d_spr1_factor(.FALSE.,.FALSE.,.TRUE.,N,D,E,U,Q,D1,C1,B1,scale,1,Z,DWORK,INFO)
  call d_spr1_factor(.FALSE.,.FALSE.,SCA,N,D,E,U,Q,D1,C1,B1,scale,1,Z,DWORK,INFO)
  !call d_spr1_factor2(.FALSE.,.FALSE.,.FALSE.,N,D,E,U,Q,D1,C1,B1,scale,1,Z,INFO)

  ! check info
  if (INFO.NE.0) then 
     ! print error in debug mode
     if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"d_symtrid_factor failed",INFO,INFO)
     end if
     INFO = 1
     ! no eigenvalues found, no back transform necessary
     return
  end if
  
  ! compute the roots with upr1fact Hessenberg QR
  !call z_upr1fact_twistedqz(.FALSE.,.FALSE.,.FALSE.,l_upr1fact_hess,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
  call z_upr1fact_qr(.FALSE.,.FALSE.,l_upr1fact_hess,N,P,Q,D1,C1,B1,N,V,ITS,INFO)

  
  ! check INFO
  ! if INFO == 1 some eigenvalues might have been found
  if ((INFO.NE.0).AND.(INFO.NE.1)) then
    print*,""
    print*,"INFO:",INFO
    print*,""
    INFO = 2
    deallocate(P,ITS,Q,D1,C1,B1,D,E,U,Z)!,H)
    deallocate(DWORK,D2,C2,B2)    
    return  
  end if

  ! extract roots
  call z_upr1utri_decompress(.TRUE.,N,D1,C1,B1,ROOTS)
  !call z_upr1fact_extracttri(.TRUE.,N,D1,C1,B1,ROOTS)
  
  ! back transformation
  do ii=1,N
!!$     if (N.LE.16) then
!!$        print*, ii,"ROOTS(ii)", ROOTS(ii),"1/ROOTS(ii)", 1d0/ROOTS(ii),"nrm",scale           
!!$     end if
!     if (abs(1d0-scale)<EISCOR_DBL_EPS) then
!        ROOTS(ii) = aimag(ROOTS(ii))/(1d0+dble(ROOTS(ii)))*scale
!     else
        ROOTS(ii) = cmplx(0d0,1d0,kind=8)*(cmplx(1d0,0d0,kind=8)-ROOTS(ii))/(cmplx(1d0,0d0,kind=8)+ROOTS(ii))*scale
 !    end if

  end do
    
  ! compute residuals
  ! call z_poly_residuals(N,COEFFS,ROOTS,0,RESIDUALS)

  deallocate(P,ITS,Q,D1,C1,B1,D,E,U,Z)!,H)
  deallocate(DWORK,D2,C2,B2)    

    
end subroutine d_polyc_roots
