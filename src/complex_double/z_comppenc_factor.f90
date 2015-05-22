#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_comppenc_factor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  QZ              LOGICAL
!                    .TRUE.: second triangular factor is used
!                    .FALSE.: second triangular factor is assumed 
!                             to be identity
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  V               COMPLEX(8) array of dimension (N)
!                    coefficients for left triangular factor
!
!  W               COMPLEX(8) array of dimension (N)
!                    coefficients for right traingular factor
!
! OUTPUT VARIABLES:
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D1,D2           REAL(8) arrays of dimension (2*(N+1))
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1,B1,C2,B2     REAL(8) arrays of dimension (3*N)
!                    array of generators for upper-triangular parts
!                    of the pencil
!
!  INFO           INTEGER
!                   INFO = 1 implies no convergence 
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_comppenc_factor(QZ,N,P,V,W,Q,D1,C1,B1,D2,C2,B2,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: QZ
  integer, intent(in) :: N
  logical, intent(in) :: P(N-2)
  complex(8), intent(inout) :: V(N), W(N)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*(N+1)), C1(3*N), B1(3*N)
  real(8), intent(inout) :: D2(2*(N+1)), C2(3*N), B2(3*N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii
  real(8) :: phr, phi, nrm, beta
  complex(8) :: temp, block(2,2), vec(2,1), h(2,2)
  
  ! initialize info
  INFO = 0
  
  ! initialize Q
  Q = 0d0
  do ii = 1,(N-1)
    Q(3*ii) = 1d0
  end do

  ! shuffle V
  do ii = 1,(N-2)

    ! permute if necessary
    if (P(N-ii-1)) then

      temp = (-1d0)**(N-ii-1)*V(N-ii)
      V(2:(N-ii)) = V(1:(N-ii-1))
      V(N-ii) = temp    

    end if
 
  end do

  ! compress first triangle
  D1 = 0d0
  do ii = 1,N
    D1(2*ii-1) = 1d0
  end do
  B1 = 0d0
  C1 = 0d0

  ! compute the phase of last coefficient
  call d_rot2_vec2gen(dble(V(N)),aimag(V(N)),phr,phi,beta)
 
  ! compute square root of phase
  temp = sqrt(cmplx(phr,phi,kind=8))  
  call d_rot2_vec2gen(dble(temp),aimag(temp),phr,phi,nrm)
  temp = cmplx(-phr,-phi,kind=8)

  ! initialize bottom of B1
  call z_rot3_vec3gen(-beta*phr,-beta*phi,1d0 &
  ,B1(3*N-2),B1(3*N-1),B1(3*N),nrm)

  ! initialize bottom of C1
  C1(3*N-2) = B1(3*N-2)
  C1(3*N-1) = -B1(3*N-1)
  C1(3*N) = -B1(3*N)

  ! update bottom of B1
  call z_rot3_vec3gen(C1(3*N)*phr,C1(3*N)*phi,-beta/nrm &
  ,B1(3*N-2),B1(3*N-1),B1(3*N),nrm)

  ! roll up V into B1 and C1
  do ii = 1,(N-1)

    ! update last entry of V using C1
    temp = cmplx(C1(3*(N-ii+1)-2),C1(3*(N-ii+1)-1),kind=8)*V(N-ii+1) &
    - C1(3*(N-ii+1))*temp

    ! compute new B1
    call z_rot3_vec4gen(dble(V(N-ii)),aimag(V(N-ii)),dble(temp),aimag(temp) &
    ,B1(3*(N-ii)-2),B1(3*(N-ii)-1),B1(3*(N-ii)),nrm)

    ! store new C1
    C1(3*(N-ii)-2) = B1(3*(N-ii)-2)
    C1(3*(N-ii)-1) = -B1(3*(N-ii)-1)
    C1(3*(N-ii)) = -B1(3*(N-ii))

  end do

  ! shuffle W
  if (QZ) then

    do ii = 1,(N-2)

      ! permute if necessary
      if (P(N-ii-1)) then

        temp = (-1d0)**(N-ii-1)*W(N-ii)
        W(2:(N-ii)) = W(1:(N-ii-1))
        W(N-ii) = temp    

      end if
 
    end do

    ! compress second triangle
    D2 = 0d0
    do ii = 1,N
      D2(2*ii-1) = 1d0
    end do
    B2 = 0d0
    C2 = 0d0

    ! compute the phase of last coefficient
    call d_rot2_vec2gen(dble(W(N)),aimag(W(N)),phr,phi,nrm)

    ! compute square root of phase
    temp = sqrt(cmplx(phr,phi,kind=8))  
    call d_rot2_vec2gen(dble(temp),aimag(temp),phr,phi,nrm)

    ! initialize bottom of B2
    call z_rot3_vec4gen(dble(W(N)),aimag(W(N)),-phr,-phi &
    ,B2(3*N-2),B2(3*N-1),B2(3*N),nrm)

    ! initialize bottom of C2
    C2(3*N-2) = B2(3*N-2)
    C2(3*N-1) = -B2(3*N-1)
    C2(3*N) = -B2(3*N)

    ! update bottom of B2
    call z_rot3_vec4gen(C2(3*N)*phr,C2(3*N)*phi,C2(3*N-2)*phr+C2(3*N-1)*phi &
    ,C2(3*N-2)*phi-C2(3*N-1)*phr,B2(3*N-2),B2(3*N-1),B2(3*N),nrm)

    ! roll up W into B2 and C2
    temp = cmplx(-phr,-phi,kind=8)
    do ii = 1,(N-1)

      ! update last entry of W using C2
      temp = cmplx(C2(3*(N-ii+1)-2),C2(3*(N-ii+1)-1),kind=8)*W(N-ii+1) &
      - C2(3*(N-ii+1))*temp

      ! compute new B2
      call z_rot3_vec4gen(dble(W(N-ii)),aimag(W(N-ii)),dble(temp),aimag(temp) &
      ,B2(3*(N-ii)-2),B2(3*(N-ii)-1),B2(3*(N-ii)),nrm)

      ! store new C2
      C2(3*(N-ii)-2) = B2(3*(N-ii)-2)
      C2(3*(N-ii)-1) = -B2(3*(N-ii)-1)
      C2(3*(N-ii)) = -B2(3*(N-ii))

    end do

  end if

end subroutine z_comppenc_factor
