#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_singlestep 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift 
! algorithm on a upr1 pencil. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  QZ              LOGICAL
!                    .TRUE. second triangular factor is assumed nonzero
!                    .FALSE. second triangular factor is assumed to be identity
!
!  VEC             LOGICAL
!                    .TRUE. update schurvectors
!                    .FALSE. no schurvectors
!
!  FUN             LOGICAL FUNCTION FUN(N,P)
!                    takes integer N and logical array P of 
!                    dimension N-2 and outputs a logical 
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (4*(N+1))
!                    array of generators for givens rotations
!
!  C1,B1,C2,B2     REAL(8) arrays of dimension (3*N)
!                    array of generators for upper-triangular parts of the pencil
!                    If QZ = .FALSE., C2 and B2 are unused.
!
!  M               INTEGER
!                    leading dimesnion of V and W
!
!  V               COMPLEX(8) array of dimension (M,N)
!                    right schur vectors
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. update V to store right schurvectors 
!
!  W               COMPLEX(8) array of dimension (M,N)
!                    left schur vectors
!                    if QZ = .FALSE. unused
!                    if VEC = .FALSE. unused
!                    if VEC = .TRUE. update W to store left schurvectors
!
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_singlestep(QZ,VEC,FUN,N,P,Q,C1,B1,C2,B2,M,V,W,ITCNT)

  implicit none
  
  ! input variables
  logical, intent(in) :: QZ, VEC
  integer, intent(in) :: M, N
  logical, intent(inout) :: P(N)
  real(8), intent(inout) :: Q(4*(N+1)), C1(3*N), B1(3*N), C2(3*N), B2(3*N)
  complex(8), intent(inout) :: V(M,N), W(M,N)
  integer, intent(in) :: ITCNT
  interface
    function FUN(m,flags)
      logical :: FUN
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function FUN
  end interface
  
  ! compute variables
  integer :: ii, ir1, ir2, id1, id2
  logical :: final_flag
  real(8) :: G1(4), G2(4), G3(4)
  complex(8) :: shift
  logical :: tP(3)
  real(8) :: tQ(12), tC1(9), tB1(9), tC2(9), tB2(9)
  
  ! compute final_flag
  final_flag = FUN(N,P(2:(N-1)))

  ! correct top end of Q
  call z_upr1fact_correctend(.TRUE.,P(1:2),Q(1:8))
  
  ! compute shift
  ! random shift
  if ((mod(ITCNT+1,16) == 0).OR.(ITCNT.LT.0)) then
    call random_number(G1(1))
    call random_number(G1(2))
    shift = cmplx(G1(1),G1(2),kind=8)
          
  ! wilkinson shift
  else

    tC2 = 0d0
    tB2 = 0d0
  
    ! N < 3
    if (N.LT.3) then

      tP(1) = .FALSE.; tP(2:3) = P((N-1):N)
      tQ = 0d0; tQ(1) = 1d0; tQ(5:12) = Q(5:12)
      tC1 = 0d0; tC1(3) = 1d0; tC1(4:9) = C1
      tB1 = 0d0; tB1(3) = -1d0; tB1(4:9) = B1

      if (QZ) then
        tC2 = 0d0; tC2(3) = 1d0; tC2(4:9) = C2
        tB2 = 0d0; tB2(3) = -1d0; tB2(4:9) = B2
      end if

    ! N > 2
    else
  
      tP = P((N-2):N)
      tQ = Q((4*(N+1)-11):(4*(N+1)))
      tC1 = C1((3*N-8):(3*N))
      tB1 = B1((3*N-8):(3*N))

      if (QZ) then
        tC2 = C2((3*N-8):(3*N))
        tB2 = B2((3*N-8):(3*N))
      end if

    end if

    ! compute shift
    call z_upr1fact_singleshift(QZ,tP,tQ,tC1,tB1,tC2,tB2,shift)
  
  end if

  ! build bulge
  call z_upr1fact_buildbulge(QZ,P(1:2),Q(1:8),C1(1:6),B1(1:6) &
  ,C2(1:6),B2(1:6),shift,G1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! iteration for QZ
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (QZ) then
  
    ! initialize chase
    ! inverse hess
    if(P(2)) then

! this warning is temporary and should be removed when the code is
! finalized
print*,"" 
print*,"" 
print*,"Inside z_upr1fact_singlestep."
print*,"  P(2) == .TRUE. is not supported!" 
print*,"" 
print*,"" 

    ! hess
    else

      ! G2 = G1 = G1^-1
      G1(2) = -G1(2); G1(3) = -G1(3)
      G2 = G1

      ! merge G1 with Q
      call z_upr1fact_mergebulge(.TRUE.,P(1:2),Q(1:8),G1)

      ! pass G2 from left to right through R2
      call z_upr1fact_rot3throughtri(.TRUE.,C2(1:6),B2(1:6),G2)

      ! G2 = G2^-1
      G2(2) = -G2(2); G2(3) = -G2(3)

      ! pass G2 from right to left through R1
      call z_upr1fact_rot3throughtri(.FALSE.,C1(1:6),B1(1:6),G2)

      ! prepare G1, G2 and G3 for turnover
      G1 = Q(5:8)
      G3 = G2
      G2 = Q(9:12)      

    end if

    ! chase bulge
    do ii=1,(N-2)

! this warning is temporary and should be removed when the code is
! finalized
print*,"" 
print*,"" 
print*,"Inside z_upr1fact_singlestep."
print*,"  N > 2 is not supported!" 
print*,"" 
print*,"" 

    end do  
  
    ! set new position flag
    P(N-1) = final_flag

    ! finalize chase
    ! inverse hess
    if(P(N-1)) then

! this warning is temporary and should be removed when the code is
! finalized
print*,"" 
print*,"" 
print*,"Inside z_upr1fact_singlestep."
print*,"  P(N-1) == .TRUE. is not supported!" 
print*,"" 
print*,"" 

    ! hess
    else

      ! correct top end of Q
      call z_upr1fact_correctend(.FALSE.,P((N-1):N),Q((4*(N+1)-7):(4*(N+1))))
  
      ! merge G3 with Q
      call z_upr1fact_mergebulge(.FALSE.,P((N-1):N),Q((4*(N+1)-7):(4*(N+1))),G3)

    end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! iteration for QR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  else
 
! this warning is temporary and should be removed when the code is
! finalized
print*,"" 
print*,"" 
print*,"Inside z_upr1fact_singlestep."
print*,"  QZ == .FALSE. is not supported!" 
print*,"" 
print*,"" 
    
  end if
  
end subroutine z_upr1fact_singlestep
