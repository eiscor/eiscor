#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_hessunihess_qz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generalized Schur decomposition of an Hessenberg-
! unitary Hessenberg pencil, where the first matrix is dense and the second
! matrix is stored in factored form.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: compute Schurvector
!                    .FALSE.: no Schurvectors
!
!  ID              LOGICAL
!                    .TRUE.: initialize V and W to I
!                    .FALSE.: assume V and W already initialized
!
!  N               INTEGER
!                    dimension of matrix
!
!  A               REAL(8) array of dimension (N,N)
!                    array of the dense Hessenberg matrix
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  M               INTEGER
!                    leading dimension of V and W
!
! OUTPUT VARIABLES:
!
!  V,W             COMPLEX(8) array of dimension (M,N)
!                    right and left Schurvectors 
!
!  ITS             INTEGER array of dimension (N-1)
!                    Contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO = 2 implies no convergence 
!                    INFO = 1 random seed initialization failed
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies N is invalid
!                    INFO = -2 implies Q is invalid
!                    INFO = -3 implies V is invalid
!                    INFO = -4 implies W is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_hessunihess_qz(VEC,ID,N,A,Q,M,V,W,ITS,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, ID
  integer, intent(in) :: N, M
  real(8), intent(inout) :: A(N,N), Q(3*(N-1))
  complex(8), intent(inout) :: V(M,N), W(M,N)
  integer, intent(inout) :: INFO, ITS(N-1)
  
  ! compute variables
  logical :: flg
  integer :: ii, kk
  integer :: STR, STP, ZERO, ITMAX, ITCNT
  
  ! initialize info
  INFO = 0
  
  ! initialize random seed
  call u_randomseed_initialize(INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__, & 
           "Failed to initialize random seed",INFO,INFO)
    end if
    INFO = 1
    return
  end if
  
  ! check N
  if (N < 2) then
    INFO = -1
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
    end if
    return
  end if
  
  ! check Q 
  call z_rot3array_check(N-1,Q,flg)
  if (.NOT.flg) then
    INFO = -2
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,INFO)
    end if
    return
  end if

  ! check V
  if (VEC.AND..NOT.ID) then
    call z_2Darray_check(M,N,V,flg)
    if (.NOT.flg) then
      INFO = -3
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"V is invalid",INFO,INFO)
      end if
      return
    end if
  end if

  ! check W
  if (VEC.AND..NOT.ID) then
    call z_2Darray_check(M,N,W,flg)
    if (.NOT.flg) then
      INFO = -4
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"W is invalid",INFO,INFO)
      end if
      return
    end if
  end if

  ! initialize storage
  ITS = 0
  
  if (VEC.AND.ID) then
    V = cmplx(0d0,0d0,kind=8)
    do ii=1,min(m,n)
      V(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if   
  
  if (VEC.AND.ID) then
    W = cmplx(0d0,0d0,kind=8)
    do ii=1,min(m,n)
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

    ! check for deflation
    !call z_upr1fpen_deflationcheck(VEC,STP-STR+2,P(STR:(STP-1)), &
    !     Q((3*STR-2):(3*STP)),D1((2*STR-1):(2*STP+2)),C1((3*STR-2):(3*STP+3)), & 
    !     B1((3*STR-2):(3*STP+3)),D2((2*STR-1):(2*STP+2)),C2((3*STR-2):(3*STP+3)), & 
    1     B2((3*STR-2):(3*STP+3)),M,V(:,STR:STP+1),W(:,STR:STP+1),ZERO)
    
    ! if 1x1 block remove and check again 
    if(STP == (STR+ZERO-1))then
    
      ! update indices
      ITS(STR+ZERO-1) = ITCNT
      ITCNT = 0
      STP = STP - 1
      ZERO = 0
      STR = 1
    
    ! if greater than 1x1 chase a bulge
    else

      ! update STR
      if (ZERO.GT.0) then
         STR = STR + ZERO
         ZERO = 0
         ITS(STR+ZERO-1) = ITCNT
         ITCNT = 0
      end if

      ! perform singleshift iteration
      !call z_upr1fpen_singlestep(VEC,FUN,STP-STR+2,P(STR:(STP-1)),Q((3*STR-2):(3*STP)) &
      !,D1((2*STR-1):(2*STP+2)),C1((3*STR-2):(3*STP+3)),B1((3*STR-2):(3*STP+3)) &
      !,D2((2*STR-1):(2*STP+2)),C2((3*STR-2):(3*STP+3)),B2((3*STR-2):(3*STP+3)) &
      !,M,V(:,STR:(STP+1)),W(:,STR:(STR+1)),ITCNT)
     
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

end subroutine z_hessunihess_qz
