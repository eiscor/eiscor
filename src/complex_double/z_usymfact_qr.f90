#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_usymfact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine diagonalizes a unitary upper hessenberg matrix that is stored as
! a product of N Givens rotations.
!
! | u1       -v1 |
! | v1  conj(u1) | | u2       -v2 | 
!                  | v2  conj(u2) | | u3       -v3 | | 1   0 |
!                                   | v3  conj(u3) | | 0  u4 |                                   
!                                                                     
! The arrays U and V contain the following:
!
!  U(i) = ui
!  V(i) = vi
!
! The input must satisfy the following:
!
!  |U(i)|^2 + V(i)^2 = 1
!               V(N) = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  U               COMPLEX(8) array of dimension N
!                    array of complex generators for Givens rotations
!                    on output contains eigenvalues
!
!  V               REAL(8) array of dimension N
!                    array of real generators for Givens rotations
!
! OUTPUT VARIABLES:
!
!  ITS             INTEGER array of dimension N-1
!                    contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO =  1 implies no convergence
!                    INFO =  0 implies successful computation
!                    INFO = -1 implies N is invalid
!                    INFO = -2 implies U or VV is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_usymfact_qr(N,U,V,ITS,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  complex(8), intent(inout) :: U(N)
  real(8), intent(inout) :: V(N)
  integer, intent(inout) :: INFO, ITS(N-1)
  
  ! compute variables
  integer :: ii, kk 
  integer :: STR, STP, ZERO, ITMAX, ITCNT
  real(8) :: xx
  complex(8) :: nu
  
  ! initialize info
  INFO = 0
  
  ! check N
  if (N < 2) then
    INFO = -1
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N must be >= 2.",INFO,INFO)
    end if
    return
  end if
  

  ! check U and V
  do ii = 1,N
    if (abs(abs(U(ii))**2+V(ii)**2-1d0) > 10d0*EISCOR_DBL_EPS) then
      INFO = -2
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"U or V is invalid.",INFO,INFO)
      end if
      return
    end if
  end do
 
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
    call z_usymfact_deflationcheck(STP-STR+1,U(STR:STP),V(STR:STP),ZERO)
    
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

      ! set nu for top deflations
      if (STR > 1) then
        nu = conjg(U(STR-1))
      else
        nu = cmplx(1d0,0d0,kind=8)
      end if

      ! perform singleshift iteration
      call z_usymfact_singlestep(STP-STR+2,U(STR:STP+1),V(STR:STP+1),nu,ITCNT)
     
      ! update indices
      ITCNT = ITCNT + 1
 
    end if
    
    ! if ITMAX hit
    if (kk == ITMAX) then
      INFO = 1
      ITS(STR+STP-1) = ITCNT
    end if
    
  end do

end subroutine z_usymfact_qr
