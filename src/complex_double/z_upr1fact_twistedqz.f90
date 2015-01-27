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
! The matrices V and W are the right and left Schur vectors respectively.
! Namely, W*(H,S)V is upper-triangular.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  ALG             CHARACTER(2)
!                    'QR': second triangular factor is assumed to be identity
!                    'QZ': second triangular factor is assumed nonzero
!
!  COMPZ           CHARACTER
!                    'N': no schurvectors
!                    'I': schurvectors, initializing V and W to the identity
!                    'V': schurvectors, assume V and W already initialized
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  FUN             LOGICAL FUNCTION FUN(N,P)
!                    takes integer N and logical array P of 
!                    dimension N-2 and outputs a logical 
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D               REAL(8) array of dimension (2,2*(N+1))
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!                    D1 = D(1,:)
!                    D2 = D(2,:)
!
!  R               REAL(8) array of dimension (4,3*N)
!                    array of generators for upper-triangular parts
!                    of the pencil
!                    C1 = R(1,:)
!                    B1 = R(2,:)
!                    C2 = R(3,:)
!                    B2 = R(4,:)
!
! OUTPUT VARIABLES:
!
!  V              COMPLEX(8) array of dimension (N,N)
!                   right schur vectors
!                   if COMPZ = 'N' unused
!                   if COMPZ = 'I' stores right schurvectors in V 
!                   if COMPZ = 'V' update V to store right schurvectors 
!
!  W              COMPLEX(8) array of dimension (N,N)
!                   left schur vectors
!                   if COMPZ = 'N' unused
!                   if COMPZ = 'I' stores left schurvectors in W 
!                   if COMPZ = 'V' update W to store left schurvectors 
!
!  ITS            INTEGER array of dimension (N-1)
!                   Contains the number of iterations per deflation
!
!  INFO           INTEGER
!                   INFO = 1 implies no convergence 
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies ALG, N, Q, D or R is invalid
!                   INFO = -2 implies COMPZ is invalid
!                   INFO = -9 implies V is invalid
!                   INFO = -10 implies W is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_twistedqz(ALG,COMPZ,N,P,FUN,Q,D,R,V,W,ITS,INFO)

  implicit none
  
  ! input variables
  character(2), intent(in) :: ALG
  character, intent(in) :: COMPZ
  integer, intent(in) :: N
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D(2,2*(N+1)), R(4,3*N)
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
  integer :: start_index, stop_index, zero_index, it_max, it_count
  
  ! initialize info
  INFO = 0
  
  ! check factorization
  call z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"ALG, N, Q, D or R is invalid",INFO,INFO)
    end if
    INFO = -1
    return
  end if
  
  ! check COMPZ
  if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
    INFO = -2
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"COMPZ must be 'N', 'I' or 'V'",INFO,INFO)
    end if
    return
  end if
  
  ! check V
  if (COMPZ.EQ.'V') then
    call z_2Darray_check(N,N,V,INFO)
    if (INFO.NE.0) then
      INFO = -9
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"V is invalid",INFO,INFO)
      end if
      return
    end if
  end if   
  
  ! check W
  if ((ALG.EQ.'QZ').AND.(COMPZ.EQ.'V')) then
    call z_2Darray_check(N,N,W,INFO)
    if (INFO.NE.0) then
      INFO = -10
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"W is invalid",INFO,INFO)
      end if
      return
    end if
  end if
  
  ! initialize storage
  ITS = 0
  
  if (COMPZ.EQ.'I') then
    V = cmplx(0d0,0d0,kind=8)
    do ii=1,n
      V(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if   
  
  if ((ALG.EQ.'QZ').AND.(COMPZ.EQ.'I')) then
    W = cmplx(0d0,0d0,kind=8)
    do ii=1,n
      W(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if   
 
  ! initialize indices
  start_index = 1
  stop_index = N-1
  zero_index = 0
  it_max = 20*N
  it_count = 0
  
  ! iteration loop
  do kk=1,it_max
  
    ! check for completion
    if(stop_index <= 0)then    
      exit
    end if
    
    ! check for deflation
    call z_upr1fact_deflationcheck(N,start_index,stop_index,zero_index,P,Q,D(1,:),it_count,ITS,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_deflationcheck failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! if 1x1 block remove and check again 
    if(stop_index == zero_index)then
    
      ! update indices
      stop_index = stop_index - 1
      zero_index = 0
      start_index = 1
    
    ! if 2x2 block remove and check again
    else if(stop_index-1 == zero_index)then
    
      ! call 2x2 deflation
      call z_upr1fact_2x2deflation(ALG,COMPZ,N,stop_index,P,Q,D,R,V,W,INFO)
    
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_2x2deflation failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
      
      ! update indices
      stop_index = stop_index - 2
      zero_index = 0
      start_index = 1
    
    ! if greater than 2x2 chase a bulge
    else
      
      ! perform singleshift iteration
      call z_upr1fact_singlestep(ALG,COMPZ,N,start_index,stop_index,P,FUN,Q,D,R,V,W,it_count,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_upr1fact_singlestep failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
    
    end if
    
    ! if it_max hit
    if (kk == it_max) then
      INFO = 1
      ITS(stop_index) = it_count
    end if
    
  end do

end subroutine z_upr1fact_twistedqz
