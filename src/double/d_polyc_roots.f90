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
subroutine d_polyc_roots(N,COEFFS,ROOTS,RESIDUALS)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(in) :: COEFFS(N)
  complex(8), intent(inout) :: ROOTS(N)
  real(8), intent(inout) :: RESIDUALS(N)
  
  ! compute variables
  integer :: ii, INFO
  real(8) :: scl
  logical, allocatable :: P(:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:),D1(:),C1(:),B1(:),D(:),E(:),Z(:,:)
  real(8) :: D2,C2,B2,scale
  complex(8), allocatable :: U(:), H(:,:)
  complex(8) :: V,W, t(2,2)

  complex(8) :: WORK(5*N)
  real(8) :: RWORK(2*N)

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
  allocate(D(N),E(N-1),U(N),Z(1,N),H(N,N))

  ! fill P
  P = .FALSE.

  ! set up D, E, and U
  D = 0d0
  E = 5d-1
  E(1) = sqrt(2d0)/2d0
  U(1) = -sqrt(2d0)/COEFFS(N-1+1)/2d0
  do ii=2,N
     U(ii) = -COEFFS(N-ii+1)/2d0
  end do
  
  ! factorize \Phi(T + Ue^H) and reduce it to Hessenberg form
  call d_spr1_factor(.FALSE.,.FALSE.,.FALSE.,N,D,E,U,Q,D1,C1,B1,scale,1,Z,INFO)

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
  
!!$  ! compute matrix H
!!$  H = 0d0
!!$  call z_upr1fact_extracttri(.FALSE.,N,D1,C1,B1,H)
!!$  ! plot H
!!$  do ii=1,N
!!$     print*, ii, H(ii,:)
!!$  end do
!!$
!!$  do ii=N-1,1,-1
!!$     t(1,1) = cmplx(Q(3*ii-2),Q(3*ii-1),kind=8)
!!$     t(2,1) = cmplx(Q(3*ii),0d0,kind=8)
!!$     t(1,2) = -t(2,1)
!!$     t(2,2) = conjg(t(1,1))
!!$     H(ii:(ii+1),:) = matmul(t,H(ii:(ii+1),:))
!!$  end do
!!$  print*, "reconstructed H"
!!$  ! plot H
!!$  do ii=1,N
!!$     print*, ii, H(ii,:)
!!$  end do
!!$
!!$  call zgeev('N','N', N, H, N, ROOTS, Z, 1, Z, 1, WORK, 5*N, RWORK, INFO)
!!$  
!!$  do ii=1,N
!!$     print*, ii, ROOTS(ii), cmplx(0d0,1d0,kind=8)*(cmplx(1d0,0d0,kind=8)-ROOTS(ii))/&
!!$          &(cmplx(1d0,0d0,kind=8)+ROOTS(ii))
!!$  end do
!!$  print*, "START z_upr1fact_twistedqz"

  ! compute the roots with upr1fact Hessenberg QR
  call z_upr1fact_twistedqz(.FALSE.,.FALSE.,.FALSE.,l_upr1fact_hess,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
  ! call z_upr1fact_twistedqz(.TRUE.,.FALSE.,.FALSE.,l_upr1fact_inversehess,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
  ! call z_upr1fact_twistedqz(.TRUE.,.FALSE.,.FALSE.,l_upr1fact_cmv,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
  ! call z_upr1fact_twistedqz(.TRUE.,.FALSE.,.FALSE.,l_upr1fact_random,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)
  ! call z_upr1fact_twistedqz(.FALSE.,.FALSE.,.FALSE.,l_upr1fact_hess,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)

  ! check INFO
  if (INFO.NE.0) then
    print*,""
    print*,"INFO:",INFO
    print*,""
    deallocate(P,ITS,Q,D1,C1,B1,D,E,U,Z,H)
    return  
  end if

  ! extract roots
  call z_upr1fact_extracttri(.TRUE.,N,D1,C1,B1,ROOTS)
  ! back transformation
  do ii=1,N
     !ROOTS(ii) = aimag(ROOTS(ii))/(1d0+dble(ROOTS(ii)))
     ROOTS(ii) = cmplx(0d0,1d0,kind=8)*(1-ROOTS(ii))/(1+ROOTS(ii))
  end do
    
  ! compute residuals
  ! call z_poly_residuals(N,COEFFS,ROOTS,0,RESIDUALS)

  deallocate(P,ITS,Q,D1,C1,B1,D,E,U,Z,H)

    
end subroutine d_polyc_roots
