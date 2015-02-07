!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_poly_factorcomp 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a compressed QR factorization for the 
! companion matrix of a polynomial expressed in the monomial basis.
! This factorization is stored as a product of three sequences of 
! Givens' rotations and a complex unimodular diagonal matrix. 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  QZ              LOGICAL
!                    .TRUE.: second triangular factor is assumed nonzero
!                    .FALSE.: second triangular factor is assumed to be identity
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvector
!                    .FALSE.: no schurvectors
!
!  ID              LOGICAL
!                    .TRUE.: initialize V and W to I
!                    .FALSE.: assume V and W already initialized
!
!  N               INTEGER
!                    dimension of matrix
!
!  COEFFS          COMPLEX(8) array of dimension (N)
!                   coefficients of polynomial, assumed to be of degree
!
! OUTPUT VARIABLES:
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of Givens' 
!
!  D1,D2           REAL(8) arrays of dimension (2*(N+1))
!                    array of generators for complex diagonal matrix
!
!  C1,B1,C2,B2     REAL(8) arrays of dimension (3*N)
!                    array of generators for second sequence of Givens' 
!
!  V,W             COMPLEX(8) array of dimension (N,N)
!                    array of generators for third sequence of Givens' 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_poly_factorcomp(QZ,VEC,ID,N,COEFFS,P,Q,D1,C1,B1,D2,C2,B2,V,W)

  implicit none
  
  ! input variables
  logical, intent(in) :: QZ, VEC, ID
  integer, intent(in) :: N
  complex(8), intent(in) :: COEFFS(N+1)
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*(N+1)), C1(3*N), B1(3*N)
  real(8), intent(inout) :: D2(2*(N+1)), C2(3*N), B2(3*N)
  complex(8), intent(inout) :: V(N,N), W(N,N)
  
  ! compute variables
  integer :: ii, ind
  real(8) :: phr, phi, nrm
  complex(8) :: t1, t2, phase
 

  ! QZ factorization
  if (QZ) then 



  ! QR factorization
  else

    ! set P
    P = .FALSE.

    ! initialize V
    if (VEC.AND.ID) then
      V = cmplx(0d0,0d0,kind=8)
      do ii=1,N
        V(ii,ii) = cmplx(1d0,0d0,kind=8)
      end do
    end if

    ! compute the phase of constant coefficient
    call d_rot2_vec2gen(dble(COEFFS(N+1)/COEFFS(1)),aimag(COEFFS(N+1)/COEFFS(1)) &
    ,phr,phi,nrm)
  
    ! update V
    if (VEC) then
      V(:,N) = cmplx(phr,phi,kind=8)*V(:,N)
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
    do ii=1,(N+1)
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
    call d_rot2_vec2gen(nrm*(-1d0)**(N),(-1d0)**(N-1),C1(ind+1),C1(ind+3),nrm)
    C1(ind+2) = 0d0
    B1(ind+1) = C1(ind+3)*(-1d0)**(N)
    B1(ind+2) = 0d0
    B1(ind+3) = C1(ind+1)*(-1d0)**(N)
  
    do ii=2,N
      ind = 3*(N-ii+1)
      t2 = cmplx(C1(ind+1),-C1(ind+2),kind=8)*t1 + cmplx(C1(ind+3),0d0,kind=8)*t2
      t1 = -COEFFS(ii)*conjg(phase)/COEFFS(1)
      ind = 3*(N-ii)
      call z_rot3_vec4gen(dble(t1),aimag(t1),dble(t2),aimag(t2),C1(ind+1),C1(ind+2),C1(ind+3),nrm)
      B1(ind+1) = C1(ind+1)
      B1(ind+2) = -C1(ind+2)
      B1(ind+3) = -C1(ind+3)
    end do

  end if

end subroutine z_poly_factorcomp
