#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_idpspike_compress
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine factors a identity matrix plus a spike in the last
! column into D C ( B + u e_n^H), where D is a (n+1 x n+1)unitary 
! diagonal matrix, C an ascending sequence of n rotations and B a 
! descending sequence. Since the whole matrix can be represented by D,
! C, and B the updated u is not stored.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Input: U
! 1       u 
!   1     u
!     1   u
!       1 u
!         u
!
! The vector u is not reshuffeled.
! 
!
! Output: D, C, B
!
! D          C [ B         a 00001 ]
!  D        C  [  B        0       ]
!   D      C   [   B    +  0       ]
!    D    C    [    B      0       ]
!     D  C     [     B     0       ]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  ID              LOGICAL
!                    .TRUE.: initialize D to be the indentity
!                    .FALSE.: D is already initialize
!
!  N               INTEGER
!                    dimension of matrix
!
!  M               INTEGER
!                    column of the spike
!
!  U               COMPLEX(8) array of dimension (N)
!                    coefficients for left triangular factor
!
! OUTPUT VARIABLES:
!
!  D               REAL(8) arrays of dimension (2*(N+1))
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C,B            REAL(8) arrays of dimension (3*N)
!                    array of generators for upper-triangular parts
!                    of the pencil
!
!  INFO           INTEGER
!                   INFO = 0 implies successful computation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_idpspike_compress(ID,N,M,U,D,C,B,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: ID
  integer, intent(in) :: N, M
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: D(2*N), C(3*N), B(3*N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii
  real(8) :: phr, phi, nrm, beta
  complex(8) :: temp
  
  ! initialize info
  INFO = 0

  ! initialize D
  if (ID) then
     D = 0d0
     do ii = 1,N
        D(2*ii-1) = 1d0
     end do
  end if
  
  ! initialize B and C
  B = 0d0
  C = 0d0

  ! compute the phase of last coefficient
  call d_rot2_vec2gen(dble(U(M)),aimag(U(M)),phr,phi,beta)
 
  ! store in D
  if (ID) then
     D(2*M-1) = phr
     D(2*M) = phi
     !D(2*N+1) = phr
     !D(2*N+2) = -phi
  else
     temp = cmplx(D(2*M-1),D(2*M),kind=8)*cmplx(phr,phi,kind=8)
     call d_rot2_vec2gen(dble(temp),aimag(temp),D(2*M-1),D(2*M),nrm)
     !temp = cmplx(D(2*N+1),D(2*N+2),kind=8)*cmplx(phr,-phi,kind=8)
     !call d_rot2_vec2gen(dble(temp),aimag(temp),D(2*M+1),D(2*M+2),nrm)
  end if

  if (M.EQ.N) then
     ! initialize bottom of C
     call d_rot2_vec2gen(beta,-1d0,C(3*N-2),C(3*N),nrm)
     ! initialize bottom of B
     B(3*N-2) = C(3*N)
     B(3*N) = C(3*N-2)
  else
     do ii=N,M+1,-1
        C(3*ii) = 1d0
        B(3*ii) = -1d0
     end do
     ! initialize bottom of C
     call d_rot2_vec2gen(beta,-1d0,C(3*M-2),C(3*M),nrm)
     ! initialize bottom of B
     B(3*M-2) = C(3*M)
     B(3*M) = C(3*M-2)     
  end if


  ! roll up U into B and C
  temp = cmplx(nrm,0d0,kind=8)
  do ii = 1+N-M,(N-1)

    ! compute new C
    call z_rot3_vec4gen(dble(U(N-ii)),aimag(U(N-ii)),dble(temp),aimag(temp) &
         &,C(3*(N-ii)-2),C(3*(N-ii)-1),C(3*(N-ii)),nrm)

    ! store new B
    B(3*(N-ii)-2) = C(3*(N-ii)-2)
    B(3*(N-ii)-1) = -C(3*(N-ii)-1)
    B(3*(N-ii)) = -C(3*(N-ii))

    ! update last entry of U using C
    temp = cmplx(C(3*(N-ii)-2),-C(3*(N-ii)-1),kind=8)*U(N-ii) &
         &+ C(3*(N-ii))*temp

  end do

end subroutine z_idpspike_compress
