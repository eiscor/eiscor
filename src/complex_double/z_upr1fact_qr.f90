#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Schur decomposition of an 
! extended upper-hessenberg matrix. Both the matrix
! is the sum of a unitary matrix and a rank 
! one matrix. These matrices are stored in 3 sequences of rotations 
! and 1 unimodular diagonal matrix.
!
! The hessenberg part is stored as H = Q*D*C*B
!
! The matrices V are the right Schur vectors.
! Namely, V*HV is upper-triangular.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
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
!  D               REAL(8) arrays of dimension (2*(N+1))
!                    array of generators for complex diagonal matrix
!                    in the upper-triangular factor
!
!  C,B             REAL(8) arrays of dimension (3*N)
!                    array of generators for upper-triangular parts
!
!  M               INTEGER
!                    leading dimension of V matrix
!
! OUTPUT VARIABLES:
!
!  V              COMPLEX(8) array of dimension (M,N)
!                   right schur vectors
!
!  ITS            INTEGER array of dimension (N-1)
!                   Contains the number of iterations per deflation
!
!  INFO           INTEGER
!                   INFO = 2 implies no convergence 
!                   INFO = 1 random seed initialization failed
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies N, Q, D, C or B is invalid
!                   INFO = -14 implies V is invalid
!                   INFO = -15 implies W is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_qr(VEC,ID,FUN,N,P,Q,D,C,B,M,V,ITS,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, ID
  integer, intent(in) :: M,N
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D(2*(N+1)), C(3*N), B(3*N)
  complex(8), intent(inout) :: V(M,N)
  integer, intent(inout) :: INFO, ITS(N-1)
  interface
    function FUN(m,flags)
      logical :: FUN
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function FUN
  end interface
  
  ! compute variables
  logical :: flg
  integer :: ii, jj, kk
  integer :: STR, STP, ZERO, ITMAX, ITCNT
  
  ! initialize info
  INFO = 0

  ! initialize random seed
  call u_randomseed_initialize(INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,&
           "Failed to initialize random seed",INFO,INFO)
    end if
    INFO = 1
    return
  end if
 
  ! check factorization
  call z_upr1fact_factorcheck(N,Q,D,C,B,INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N, Q, D, C or B is invalid",INFO,INFO)
    end if
    INFO = -1
    return
  end if
  
  ! check V
  if (VEC.AND..NOT.ID) then
    call z_2Darray_check(M,N,V,flg)
    if (.NOT.flg) then
      INFO = -14
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"V is invalid",INFO,INFO)
      end if
      return
    end if
  end if   
  
  ! initialize storage
  ITS = 0
  
  if (VEC.AND.ID) then
    V = cmplx(0d0,0d0,kind=8)
    do ii=1,min(M,N)
      V(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if   
  
  ! initialize indices
  STR = 1
  STP = N-1
  ZERO = 0
  ITMAX = 20*N
  ITCNT = -1
  
  ! iteration loop
  do kk=1,ITMAX

    ! check for completion
    if(STP <= 0)then    
      exit
    end if

!print*,""
!print*,"Inside twisted QZ"
!print*,"STR:",STR
!print*,"STP:",STP
!print*,"Q"
!do ii=1,(N-1)
!print*,Q(3*ii-2),Q(3*ii-1),Q(3*ii)
!end do
!print*,""
!print*,"D1"
!do ii=1,(N+1)
!print*,D1(2*ii-1),D1(2*ii)
!end do
!print*,""
!print*,"C1"
!do ii=1,(N)
!print*,C1(3*ii-2),C1(3*ii-1),C1(3*ii)
!end do
!print*,""
!print*,"B1"
!do ii=1,(N)
!print*,B1(3*ii-2),B1(3*ii-1),B1(3*ii)
!end do
!print*,""
!print*,"D2"
!do ii=1,(N+1)
!print*,D2(2*ii-1),D2(2*ii)
!end do
!print*,""
!print*,"C2"
!do ii=1,(N)
!print*,C2(3*ii-2),C2(3*ii-1),C2(3*ii)
!end do
!print*,""
!print*,"B2"
!do ii=1,(N)
!print*,B2(3*ii-2),B2(3*ii-1),B2(3*ii)
!end do
!print*,""
!pause

    ! check for deflation
    call z_upr1fact_deflationcheck(STP-STR+2,P(STR:(STP-1)) &
    ,Q((3*STR-2):(3*STP)),D1((2*STR-1):(2*STP+2)),ZERO)
    
    ! if 1x1 block remove and check again 
    if(STP == (STR+ZERO-1))then
    
      ! update indices
      ITS(STR+ZERO-1) = ITCNT
      ITCNT = 0
      STP = STP - 1
      ZERO = 0
      STR = 1
    
    ! if 2x2 block remove and check again
    else if(STP == (STR+ZERO))then
    
      ! check STR
      if (STR <= ZERO) then
        STR = ZERO+1
      end if

      ! call 2x2 deflation
      call z_upr1fact_2x2deflation(QZ,VEC,Q((3*STP-2):(3*STP)) &
      ,D1((2*STP-1):(2*STP+2)),C1((3*STP-2):(3*STP+3)) &
      ,B1((3*STP-2):(3*STP+3)),D2((2*STP-1):(2*STP+2)) &
      ,C2((3*STP-2):(3*STP+3)),B2((3*STP-2):(3*STP+3)) &
      ,N,V(:,STP:(STP+1)),W(:,STP:(STP+1)))
    
      ! update indices
      ITCNT = ITCNT + 1
 
      ! update indices
!      ITS(STR+STP-1) = ITCNT
!      ITCNT = 0
!      STP = STP - 2
!      ZERO = 0
!      STR = 1
    
    ! if greater than 2x2 chase a bulge
    else

      ! check STR
      if (STR <= ZERO) then
        STR = ZERO+1
      end if

      ! perform singleshift iteration
      call z_upr1fact_singlestep(QZ,VEC,FUN,STP-STR+2,P(STR:(STP-1)) &
      ,Q((3*STR-2):(3*STP)),D1((2*STR-1):(2*STP+2)) &
      ,C1((3*STR-2):(3*STP+3)),B1((3*STR-2):(3*STP+3)) &
      ,D2((2*STR-1):(2*STP+2)),C2((3*STR-2):(3*STP+3)) &
      ,B2((3*STR-2):(3*STP+3)),N,V(:,STR:(STP+1)),W(:,STR:(STP+1)),ITCNT)
     
      ! update indices
      if (ITCNT.EQ.-1) then 
        ITCNT = 1 
      else
        ITCNT = ITCNT + 1
      end if

    end if
    
    ! if ITMAX hit
    if (kk == ITMAX) then
      INFO = 1
      ITS(STR+STP-1) = ITCNT
    end if
    
  end do

end subroutine z_upr1fact_qr
