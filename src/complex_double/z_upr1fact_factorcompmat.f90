#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_factorcompmat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a compressed QZ factorization for the 
! companion matrix of a polynomial expressed in the monomial basis. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  A1          COMPLEX(8) array of dimension (N)
!                   coefficients of polynomial, assumed to be of degree
!                   exactly N, have leading coefficient 1 and have no 
!                   zero roots
!
!  Q               REAL(8) array of dimension (3*N)
!                    array of generators for first sequence of Givens' 
!                    rotations
!
!  D1, D2          REAL(8) array of dimension (2*(N+1))
!                    arrays of generators for complex diagonal matrices
!
!  C1, B1, C2, B2  REAL(8) array of dimension (3*N)
!                    arrays of generators for second sequence of Givens' 
!                    rotations
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                   INFO = 0 implies successful computation.
!                   INFO negative implies that an input variable has
!                   an improper value, i.e. INFO=-2 => N is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_factorcompmat(ALG,COMPZ,N,A1,A2,Q,D1,C1,B1,D2,C2,B2,V,W,INFO)

  implicit none
  
  ! input variables
  character(2), intent(in) :: ALG
  character, intent(in) :: COMPZ
  integer, intent(in) :: N
  complex(8), intent(in) :: A1(N), A2(N)
  real(8), intent(inout) :: Q(3*N), D1(2*(N+1)), D2(2*(N+1))
  real(8), intent(inout) :: C1(3*N), B1(3*N), C2(3*N), B2(3*N)
  complex(8), intent(inout) :: V(N,N), W(N,N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii, ind
  real(8) :: phr, phi, nrm
  complex(8) :: t1, t2
  
  ! initialize INFO
  INFO = 0
  
  ! check input
  if (DEBUG) then

  end if      
  
  ! initialize V
  if (COMPZ.EQ.'I') then
    V = cmplx(0d0,0d0,kind=8)
    do ii=1,N
      V(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if
  
  ! QR factorization
  if (ALG.EQ.'QR') then
  
    ! compute the phase of A1(N)
    call d_rot2_vec2gen(dble(A1(N)),aimag(A1(N)),phr,phi,nrm,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! update V
    if (COMPZ.NE.'N') then
      V(:,N) = cmplx(phr,-phi,kind=8)*V(:,N)
    end if
    
    ! set Q
    do ii=1,(N-1)
      ind = 3*(ii-1)
      Q(ind+1) = 0d0
      Q(ind+2) = 0d0
      Q(ind+3) = 1d0
    end do
    ind = 3*(N-1)
    Q(ind+1) = 1d0
    Q(ind+2) = 0d0
    Q(ind+3) = 0d0
    
    ! set D1
    do ii=1,N+1
      ind = 2*(ii-1)
      D1(ind+1) = 1d0
      D1(ind+2) = 0d0
    end do
    ind = 2*(N-2)
    D1(ind+1) = phr
    D1(ind+2) = phi
    
    ! initialize B1 and C1
    t1 = cmplx(nrm*(-1d0)**(N),0d0,kind=8)
    t2 = cmplx((-1d0)**(N-1),0d0,kind=8)
    
    ind = 3*(N-1)
    call d_rot2_vec2gen(dble(t1),dble(t2),C1(ind+1),C1(ind+3),nrm,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    C1(ind+2) = 0d0
    B1(ind+1) = C1(ind+3)*(-1d0)**(N)
    B1(ind+2) = 0d0
    B1(ind+3) = C1(ind+1)*(-1d0)**(N)
    
    do ii=2,N
      ind = 3*(N-ii+1)
      t2 = cmplx(C1(ind+1),-C1(ind+2),kind=8)*t1 + cmplx(C1(ind+3),0d0,kind=8)*t2
      t1 = -A1(ii-1)*cmplx(phr,-phi,kind=8)
      ind = 3*(N-ii)
      call z_rot3_vec4gen(dble(t1),aimag(t1),dble(t2),aimag(t2),C1(ind+1),C1(ind+2),C1(ind+3),nrm,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec4gen",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if      
      
      B1(ind+1) = C1(ind+1)
      B1(ind+2) = -C1(ind+2)
      B1(ind+3) = -C1(ind+3)
    end do
  
  ! QZ factorization
  else
  
    ! initialize W
    if (COMPZ.EQ.'I') then
      W = cmplx(0d0,0d0,kind=8)
      do ii=1,N
        W(ii,ii) = cmplx(1d0,0d0,kind=8)
      end do
    end if  
  
  end if

end subroutine z_upr1fact_factorcompmat
