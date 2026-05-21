#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_singlestep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift
! algorithm on a unitary upper hessenberg matrix that is stored as a
! product of givens rotations and a complex diagonal matrix.
!
! This revised version avoids small MATMUL calls in the eigenvector
! updates.  Instead, it applies each 2 x 2 transformation to the
! corresponding two columns of Z by explicit loops.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: update eigenvectors
!                    .FALSE.: no eigenvectors
!
!  N               INTEGER
!                    dimension of matrix, must be >= 2
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  M               INTEGER
!                    leading dimension of Z
!
!  Z               COMPLEX(8) array of dimension (M,N)
!                    if VEC = .TRUE. updated
!                    if VEC = .FALSE. unused
!
!  ITCNT           INTEGER
!                    contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_singlestep(VEC,N,Q,D,M,Z,ITCNT)

  implicit none

  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: N, M
  integer, intent(inout) :: ITCNT
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N)
  complex(8), intent(inout) :: Z(M,N)

  ! compute variables
  integer :: ii, jj, ind1, ind2
  real(8) :: s1, s2, ar, ai, br, bi, nrm
  real(8) :: bulge(3), binv(3), qt(6)
  complex(8) :: shift
  complex(8) :: block(2,2), t1(2,2), t2(2,2)

  ! eigenvector update variables
  complex(8) :: z1, z2
  complex(8) :: q11, q12, q21, q22

  ! compute a nonzero shift
  ! random shift
  if (mod(ITCNT+1,11) == 0) then

    call random_number(s1)
    call random_number(s2)
    shift = cmplx(s1,s2,kind=8)

  ! Wilkinson shift
  else

    ! set generators if N == 2
    if (N.EQ.2) then
      qt(1) = 1d0
      qt(2:3) = 0d0
      qt(4:6) = Q

    ! else if N > 2
    else
      qt = Q((3*N-8):(3*N-3))
      nrm = sqrt(qt(1)**2+qt(2)**2)

      if (nrm.EQ.0d0) then
        qt(1) = 1d0
        qt(2:3) = 0d0
      else
        qt(3) = 0d0
        qt(1:2) = qt(1:2)/nrm
      end if

    end if

    ! get 2 x 2 block
    call z_unifact_2x2diagblock(.FALSE.,qt,D((2*N-3):(2*N)),block)

    ! compute eigenvalues and eigenvectors
    t1 = block
    call z_2x2array_eig(.FALSE.,t1,t1,t2,t2)

    ! choose Wilkinson shift
    ! complex abs does not matter here
    if (abs(block(2,2)-t1(1,1)) < abs(block(2,2)-t1(2,2))) then
      shift = t1(1,1)
    else
      shift = t1(2,2)
    end if

  end if

  ! project shift onto unit circle
  ar = dble(shift)
  ai = aimag(shift)
  call d_rot2_vec2gen(ar,ai,br,bi,nrm)
  shift = cmplx(br,bi,kind=8)

  ! set generators if N == 2
  if (N.EQ.2) then
    qt(1:3) = Q
    qt(4) = 1d0
    qt(5:6) = 0d0

  ! else if N > 2
  else
    qt = Q(1:6)
  end if

  ! build bulge
  call z_unifact_buildbulge(qt,D(1:4),shift,bulge)

  ! update eigenvectors
  if (VEC) then

    q11 = cmplx(bulge(1),bulge(2),kind=8)
    q21 = cmplx(bulge(3),0d0,kind=8)
    q12 = -q21
    q22 = conjg(q11)

    do jj = 1,M

      z1 = Z(jj,1)
      z2 = Z(jj,2)

      Z(jj,1) = z1*q11 + z2*q21
      Z(jj,2) = z1*q12 + z2*q22

    end do

  end if

  ! bulge inverse
  binv(1) = bulge(1)
  binv(2) = -bulge(2)
  binv(3) = -bulge(3)

  ! fusion at top
  call z_unifact_mergebulge(.TRUE.,Q(1:3),D(1:4),binv)

  ! bulge through diag
  call z_rot3_swapdiag(D(1:4),bulge)

  ! update eigenvectors
  if (VEC) then

    q11 = cmplx(binv(1),binv(2),kind=8)
    q22 = cmplx(binv(1),-binv(2),kind=8)

    do jj = 1,M
      Z(jj,1) = Z(jj,1)*q11
      Z(jj,2) = Z(jj,2)*q22
    end do

  end if

  ! update D(1:2)
  t1(1,1) = cmplx(D(1),D(2),kind=8)*cmplx(binv(1),binv(2),kind=8)
  call d_rot2_vec2gen(dble(t1(1,1)),aimag(t1(1,1)),D(1),D(2),nrm)

  ! update D(3:4)
  t1(1,1) = cmplx(D(3),D(4),kind=8)*cmplx(binv(1),-binv(2),kind=8)
  call d_rot2_vec2gen(dble(t1(1,1)),aimag(t1(1,1)),D(3),D(4),nrm)

  ! main chasing loop
  do ii = 1,(N-2)

    ! set indices
    ind1 = 3*(ii-1) + 1
    ind2 = ind1+2

    ! through Q
    call z_rot3_turnover(Q(ind1:ind2),Q((ind1+3):(ind2+3)),bulge)

    ! update eigenvectors
    if (VEC) then

      q11 = cmplx(bulge(1),bulge(2),kind=8)
      q21 = cmplx(bulge(3),0d0,kind=8)
      q12 = -q21
      q22 = conjg(q11)

      do jj = 1,M

        z1 = Z(jj,ii+1)
        z2 = Z(jj,ii+2)

        Z(jj,ii+1) = z1*q11 + z2*q21
        Z(jj,ii+2) = z1*q12 + z2*q22

      end do

    end if

    ! set indices
    ind1 = 2*ii + 1
    ind2 = ind1+3

    ! through diag
    call z_rot3_swapdiag(D(ind1:ind2),bulge)

  end do

  ! fusion at bottom
  call z_unifact_mergebulge(.FALSE.,Q((3*N-5):(3*N-3)),D((2*N-3):(2*N)),bulge)

end subroutine z_unifact_singlestep
