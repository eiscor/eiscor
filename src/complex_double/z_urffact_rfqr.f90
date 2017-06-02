#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_urffact_rfqr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Schur factorization of
! a unitary upper hessenberg matrix that is stored as a product of 
! N-1 Givens rotations and a complex diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for Givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
! OUTPUT VARIABLES:
!
!  ITS            INTEGER array of dimension (N-1)
!                   contains the number of iterations per deflation
!
!  INFO           INTEGER
!                   INFO = 1 implies no convergence
!                   INFO = 0 implies successful computation
!                   INFO = -3 implies N is invalid
!                   INFO = -4 implies Q is invalid
!                   INFO = -5 implies D is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_urffact_rfqr(N,U,VV,ITS,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: VV(N)
  integer, intent(inout) :: INFO, ITS(N-1)
  
  ! compute variables
  logical :: flg
  integer :: ii, jj, kk, ind1, ind2, ll, strt, k
  integer :: STR, STP, ZERO, ITMAX, ITCNT
  real(8) :: xx
  
  ! initialize info
  INFO = 0
  
  ! check factorization
!  call z_urffact_factorcheck(N,Q,D,INFO)
!  if (INFO.NE.0) then
!    INFO = INFO - 2
!    ! print error message in debug mode
!    if (DEBUG) then
!      call u_infocode_check(__FILE__,__LINE__,"N, Q, or D is invalid",INFO,INFO)
!    end if
!    return
!  end if
 
  ! initialize storage
  ITS = 0
  
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
      ! store eigenvalues in U
      do ii = 1,N-1
        U(N+1-ii) = conjg(U(N-ii))*U(N+1-ii)
        xx = dble(U(N+1-ii))**2 + aimag(U(N+1-ii))**2
        U(N+1-ii) = 5d-1*U(N+1-ii)*(3d0-xx)
      end do
      exit
    end if
    
    ! check for deflation
    call z_urffact_deflationcheck(STP-STR+1,U(STR:STP),VV(STR:STP),ZERO)
    
    if (ZERO.GT.0) then
      ITS(STR+ZERO-1) = ITS(STR+ZERO-1) + ITCNT
      ITCNT = 0
    end if
    
    ! if 1x1 block remove and check again 
    if(STP == (STR+ZERO-1))then
      ! update indices
      STP = STP - 1
      ZERO = 0
      STR = 1
    
    ! if greater than 1x1 chase a bulge
    else

      ! check ZERO
      if (ZERO.GT.0) then
        STR = STR+ZERO
      end if

      ! perform singleshift iteration
      call z_urffact_singlestep(STP-STR+2,U(STR:STP+1),VV(STR:STP+1),ITCNT)
     
      ! update indices
      ITCNT = ITCNT + 1
 
    end if
    
    ! if ITMAX hit
    if (kk == ITMAX) then
      INFO = 1
      ITS(STR+STP-1) = ITCNT
    end if
    
  end do

end subroutine z_urffact_rfqr
