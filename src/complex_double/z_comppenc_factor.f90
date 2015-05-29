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
!  P               LOGICAL array of dimension (N)
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
!  Q               REAL(8) array of dimension (4*(N+1))
!                    array of generators for first sequence of rotations
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
subroutine z_comppenc_factor(QZ,N,P,V,W,Q,C1,B1,C2,B2,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: QZ
  integer, intent(in) :: N
  logical, intent(inout) :: P(N)
  integer, intent(inout) :: INFO
  real(8), intent(inout) :: Q(4*(N+1)), C1(3*N), B1(3*N), C2(3*N), B2(3*N)
  complex(8), intent(inout) :: V(N), W(N)
  
  ! compute variables
  logical :: tP(2)
  integer :: ii
  real(8) :: phr, phi, nrm, beta, tQ(8)
  complex(8) :: temp
  
  ! initialize info
  INFO = 0
 
  ! fix P(1)
  P(1) = P(2)

  ! initialize Q
  Q = 0d0
  Q(1) = 1d0
  do ii = 2,N
    Q(4*ii-1) = 1d0
  end do

  ! shuffle V
  do ii = 1,(N-2)

    ! permute if necessary
    if (P(N-ii)) then

      temp = (-1d0)**(N-ii-1)*V(N-ii)
      V(2:(N-ii)) = V(1:(N-ii-1))
      V(N-ii) = temp    

    end if
 
  end do

  ! compress first triangle
  B1 = 0d0
  C1 = 0d0

  ! compute the phase of last coefficient
  call d_rot2_vec2gen(dble(V(N)),aimag(V(N)),phr,phi,beta)
 
  ! store in Q
  Q(4*N+1) = phr
  Q(4*N+2) = phi

  ! move to the correct side if necessary
  P(N) = .FALSE.
  call z_upr1fact_correctend(.FALSE.,P((N-1):N),Q((4*N-3):(4*N+4)))

  ! initialize bottom of C1
  call d_rot2_vec2gen(beta,-1d0,C1(3*N-2),C1(3*N),nrm)

  ! initialize bottom of B1
  B1(3*N-2) = C1(3*N)
  B1(3*N) = C1(3*N-2)

  ! roll up V into B1 and C1
  temp = cmplx(nrm,0d0,kind=8)
  do ii = 1,(N-1)

    ! compute new C1
    call z_rot3_vec4gen(dble(V(N-ii)),aimag(V(N-ii)),dble(temp),aimag(temp) &
    ,C1(3*(N-ii)-2),C1(3*(N-ii)-1),C1(3*(N-ii)),nrm)

    ! store new B1
    B1(3*(N-ii)-2) = C1(3*(N-ii)-2)
    B1(3*(N-ii)-1) = -C1(3*(N-ii)-1)
    B1(3*(N-ii)) = -C1(3*(N-ii))

    ! update last entry of V using C1
    temp = cmplx(C1(3*(N-ii)-2),-C1(3*(N-ii)-1),kind=8)*V(N-ii) &
    + C1(3*(N-ii))*temp

  end do

  ! shuffle W
  if (QZ) then

    do ii = 1,(N-2)

      ! permute if necessary
      if (P(N-ii)) then

        temp = (-1d0)**(N-ii-1)*W(N-ii)
        W(2:(N-ii)) = W(1:(N-ii-1))
        W(N-ii) = temp    

      end if
 
    end do

    ! compress second triangle
    B2 = 0d0
    C2 = 0d0

    ! compute the phase of last coefficient
    call d_rot2_vec2gen(dble(W(N)),aimag(W(N)),phr,phi,beta)
 
    ! store in tQ
    tQ = 0d0
    tQ(5) = phr
    tQ(6) = -phi

    ! merge into Q
    if (P(N)) then

      ! update Q
      temp = cmplx(tQ(5),tQ(6),kind=8)*cmplx(Q(4*N+1),Q(4*N+2),kind=8)
      call d_rot2_vec2gen(dble(temp),aimag(temp),Q(4*N+1),Q(4*N+2),nrm)

    else
      ! move to the correct side if necessary
      tP(1) = P(N-1); tP(2) = .TRUE.
      tQ(1:4) = Q((4*N-3):(4*N))
      call z_upr1fact_correctend(.FALSE.,tP,tQ)
      P((N-1):N) = tP
      Q((4*N-3):(4*N)) = tQ(1:4) 

      ! update Q
      temp = cmplx(tQ(5),tQ(6),kind=8)*cmplx(Q(4*N+1),Q(4*N+2),kind=8)
      call d_rot2_vec2gen(dble(temp),aimag(temp),Q(4*N+1),Q(4*N+2),nrm)
  
    end if

    ! initialize bottom of C2
    call d_rot2_vec2gen(beta,-1d0,C2(3*N-2),C2(3*N),nrm)

    ! initialize bottom of B2
    B2(3*N-2) = C2(3*N)
    B2(3*N) = C2(3*N-2)

    ! roll up W into B2 and C2
    temp = cmplx(nrm,0d0,kind=8)
    do ii = 1,(N-1)

      ! compute new C2
      call z_rot3_vec4gen(dble(W(N-ii)),aimag(W(N-ii)),dble(temp),aimag(temp) &
      ,C2(3*(N-ii)-2),C2(3*(N-ii)-1),C2(3*(N-ii)),nrm)

      ! store new B2
      B2(3*(N-ii)-2) = C2(3*(N-ii)-2)
      B2(3*(N-ii)-1) = -C2(3*(N-ii)-1)
      B2(3*(N-ii)) = -C2(3*(N-ii))

      ! update last entry of W using C2
      temp = cmplx(C2(3*(N-ii)-2),-C2(3*(N-ii)-1),kind=8)*W(N-ii) &
      + C2(3*(N-ii))*temp

    end do

  end if

end subroutine z_comppenc_factor
