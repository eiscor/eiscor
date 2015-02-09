#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_twistedqz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generalized Schur decomposition of an 
! extended upper-hessenberg, upper-triangular pencil. Both the hessenberg
! and triangular matrices are the sum of a unitary matrix and a rank 
! one matrix. These matrices are stored in 5 sequences of rotations 
! and 2 unimodular diagonal matrices.
!
! The hessenberg part is stored as H = Q*D1*C1*B1
! The triangular part is stored as S = D2*C2*B2
!
! The matrices V and W are the left and right Schur vectors respectively.
! Namely, V*(H,S)W is upper-triangular.
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
!  FUN             LOGICAL FUNCTION FUN(N,P)
!                    takes integer N and logical array P of 
!                    dimension N-2 and outputs a logical 
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
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
! OUTPUT VARIABLES:
!
!  V              COMPLEX(8) array of dimension (N,N)
!                   right schur vectors
!
!  W              COMPLEX(8) array of dimension (N,N)
!                   left schur vectors
!
!  ITS            INTEGER array of dimension (N-1)
!                   Contains the number of iterations per deflation
!
!  INFO           INTEGER
!                   INFO = 1 implies no convergence 
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies QZ, N, Q, D or R is invalid
!                   INFO = -2 implies VEC is invalid
!                   INFO = -9 implies V is invalid
!                   INFO = -10 implies W is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_twistedqz(QZ,VEC,ID,FUN,N,P,Q,D1,C1,B1,D2,C2,B2,V,W,ITS,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: QZ, VEC, ID
  integer, intent(in) :: N
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*(N+1)), C1(3*N), B1(3*N)
  real(8), intent(inout) :: D2(2*(N+1)), C2(3*N), B2(3*N)
  complex(8), intent(inout) :: V(N,N), W(N,N)
  integer, intent(inout) :: INFO, ITS(N-1)
  interface
    function FUN(m,flags)
      logical :: FUN
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function FUN
  end interface
  
  ! compute variables
  integer :: ii, jj, kk
  integer :: STR, STP, ZERO, ITMAX, ITCNT
  
  ! initialize info
  INFO = 0
  
  ! initialize storage
  ITS = 0
  
  if (VEC.AND.ID) then
    V = cmplx(0d0,0d0,kind=8)
    do ii=1,n
      V(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if   
  
  if (QZ.AND.VEC.AND.ID) then
    W = cmplx(0d0,0d0,kind=8)
    do ii=1,n
      W(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if   
 
  ! initialize indices
  STR = 1
  STP = N-1
  ZERO = 0
  ITMAX = 20*N
  ITCNT = 0
  
  ! iteration loop
  do kk=1,ITMAX

    ! check for completion
    if(STP <= 0)then    
      exit
    end if
   
print*,""
print*,""
print*,"before and after deflation"
print*,"(ZERO,STR,STP) =",ZERO,STR,STP
if (STP < 4) then
print*,"Q"
print*,Q(1:3)
print*,Q(4:6)
print*,Q(7:9)
print*,"D"
print*,D1(1:2)
print*,D1(3:4)
print*,D1(5:6)
print*,D1(7:8)
print*,"C"
print*,C1(1:3)
print*,C1(4:6)
print*,C1(7:9)
print*,C1(10:12)
print*,"B"
print*,B1(1:3)
print*,B1(4:6)
print*,B1(7:9)
print*,B1(10:12)
end if
    ! check for deflation
    call z_upr1fact_deflationcheck(STP-STR+2,P(STR:(STP-1)),Q((3*STR-2):(3*STP)) &
    ,D1((2*STR-1):(2*STP+2)),ZERO)
  
print*,"(ZERO,STR,STP) =",ZERO,STR,STP
if (STP < 4) then
print*,"Q"
print*,Q(1:3)
print*,Q(4:6)
print*,Q(7:9)
print*,"D"
print*,D1(1:2)
print*,D1(3:4)
print*,D1(5:6)
print*,D1(7:8)
print*,"C"
print*,C1(1:3)
print*,C1(4:6)
print*,C1(7:9)
print*,C1(10:12)
print*,"B"
print*,B1(1:3)
print*,B1(4:6)
print*,B1(7:9)
print*,B1(10:12)
end if
!pause

    ! update ITCNT
    if (ZERO > 0) then
      ITS(STR+STP-1) = ITCNT
      ITCNT = 0
    end if
 
    ! if 1x1 block remove and check again 
    if(STP == (STR+ZERO-1))then
    
      ! update indices
      STP = STP - 1
      ZERO = 0
      STR = 1
    
    ! if 2x2 block remove and check again
    else if(STP == (STR+ZERO))then
    !else if(.FALSE.) then
    
print*,""
print*,"before 2x2 deflation"
print*,"(ZERO,STR,STP) =",ZERO,STR,STP
 
      ! call 2x2 deflation
      call z_upr1fact_2x2deflation(QZ,VEC,Q((3*STP-2):(3*STP)),D1((2*STP-1):(2*STP+2)),C1((3*STP-2):(3*STP+3)) &
      ,B1((3*STP-2):(3*STP+3)),D2((2*STP-1):(2*STP+2)),C2((3*STP-2):(3*STP+3)),B2((3*STP-2):(3*STP+3)),N &
      ,V(:,STP:(STP+1)),W(:,STP:(STP+1)))
    
      ! update indices
!      STP = STP - 2
!      ZERO = 0
!      STR = 1
    
    ! if greater than 2x2 chase a bulge
    else

      ! check STR
      if (STR <= ZERO) then
        STR = STR+ZERO
        ZERO = 0
      end if

      ! perform singleshift iteration
      call z_upr1fact_singlestep(QZ,VEC,FUN,STP-STR+2,P(STR:(STP-1)),Q((3*STR-2):(3*STP)),D1((2*STR-1):(2*STP+2)) &
      ,C1((3*STR-2):(3*STP+3)),B1((3*STR-2):(3*STP+3)),D2((2*STR-1):(2*STP+2)),C2((3*STR-2):(3*STP+3)) &
      ,B2((3*STR-2):(3*STP+3)),N,V(:,STR:(STP+1)),W(:,STR:(STP+1)),ITCNT)
     
      ! update indices
      ITCNT = ITCNT + 1
 
    end if
    
    ! if ITMAX hit
    if (kk == ITMAX) then
      INFO = 1
      ITS(STR+STP-1) = ITCNT
    end if
    
  end do

end subroutine z_upr1fact_twistedqz
