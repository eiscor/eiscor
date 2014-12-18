#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_twistedqz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generalized schur decomposition of an 
! extended upper-hessenberg, upper-triangular pencil. Both the hessenberg
! and triangular matrices are the sum of a unitary matrix and a rank 
! one matrix. These matrices are stored in 5 sequences of rotations.
! 
! This factorization is described in:
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  ALG             CHARACTER
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
!  P               LOGICAL array of dimension (N-1)
!                    array of position flags for Q
!
!  FUN             LOGICAL FUNCTION FUN(N,P)
!                    takes integer N and logical array P of 
!                    dimension N-1 and outputs a logical 
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D               REAL(8) array of dimension (2,2*(N+1))
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  R               REAL(8) array of dimension (4,3*N)
!                    array of generators for upper-triangular parts
!                    of the pencil
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
!                   INFO = -1 implies COMPZ is invalid
!                   INFO = -2 implies ALG, N, Q, D or R is invalid
!                   INFO = -3 implies FUN is invalid
!                   INFO = -7 implies V is invalid
!                   INFO = -8 implies W is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_twistedqz(ALG,COMPZ,N,P,FUN,Q,D,R,V,W,ITS,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: ALG, COMPZ
  integer, intent(in) :: N
  logical, intent(inout) :: P(N-1)
  real(8), intent(inout) :: Q(3*(N-1)), D(2,2*(N+1)), R(4,3*N)
  complex(8), intent(inout) :: V(N,N), W(N,N)
  integer, intent(inout) :: INFO, ITS(N-1)
  interface
    logical function FUN(N,P)
      integer, intent(in) :: N
      logical, intent(in) :: P(N-1)
    end function FUN
  end interface
  
  ! compute variables
  integer :: ii
  
  ! initialize info
  INFO = 0
  
  ! check COMPZ
  if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
    INFO = -1
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"COMPZ must be 'N', 'I' or 'V'",INFO,INFO)
    end if
    return
  end if
  
  ! check factorization
  call z_upr1fact_factorcheck(N,P,Q,D,R,INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"ALG, N, Q, D or R is invalid",INFO,INFO)
    end if
    INFO = -2
    return
  end if
  
  ! check V
  if (COMPZ.EQ.'V') then
    call z_2Darray_check(N,N,V,INFO)
    if (INFO.NE.0) then
      INFO = -7
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
      INFO = -8
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
 
  ! initialize local variables
  
  ! interation loop

end subroutine z_upr1fact_twistedqz
