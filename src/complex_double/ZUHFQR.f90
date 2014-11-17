#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZUHFQR (Zomplex Unitary Hessenberg Fast QR eigensolver)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the eigenvalues and optionally eigenvectors
! of a unitary upper hessenberg matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  COMPZ           CHARACTER
!                    'N': do not compute eigenvectors
!                    'I': stores eigenvectors, initializing Z to the identity
!                    'V': stores eigenvectors, assume Z already initialized
!
!  N               INTEGER
!                    dimension of matrix
!
!  H               COMPLEX(8) array of dimension (N,N)
!                    unitary hessenberg matrix, assumed that H(ii,jj) = 0 for |ii-jj| > 0
!                    on exit contains a diagonal matrix whose entries are the 
!                    eigenvalues of H
!
!  WORK            REAL(8) array of dimension (5*N)
!                    work space for eigensolver
!
! OUTPUT VARIABLES:
!
!  Z              COMPLEX(8) array of dimension (N,N)
!                   if COMPZ = 'N' unused
!                   if COMPZ = 'I' stores eigenvectors, initializing to identity 
!                   if COMPZ = 'V' stores eigenvectors, assume Z already initialized
!
!  ITS            INTEGER array of dimension (N-1)
!                   Contains the number of iterations per deflation
!
!  INFO           INTEGER
!                   INFO = 2 implies ZUFFQR failed
!                   INFO = 1 implies ZUHRFF failed
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies COMPZ is invalid
!                   INFO = -2 implies N is invalid
!                   INFO = -3 implies H is invalid
!                   INFO = -4 implies Z is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZUHFQR(COMPZ,N,H,Z,ITS,WORK,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: COMPZ
  integer, intent(in) :: N
  real(8), intent(inout) :: WORK(5*N)
  integer, intent(inout) :: ITS(N-1), INFO
  complex(8), intent(inout) :: H(N,N), Z(N,N)
  
  ! compute variables
  integer :: ii
  
  ! initialize INFO
  INFO = 0
  
  ! check COMPZ
  if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
    INFO = -1
    ! print error in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"COMPZ must be 'N', 'I' or 'V'",INFO,INFO)
    end if
    return
  end if
  
  ! check N
  call IARNAN(N,INFO)
  if (INFO.NE.0) then
    ! print error in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"N is invalid",INFO,INFO)
    end if
    INFO = -2
    return
  end if
  call IARINF(N,INFO)
  if (INFO.NE.0) then
    ! print error in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"N is invalid",INFO,INFO)
    end if
    INFO = -2
    return
  end if
  if (N < 2) then
    INFO = -2
    ! print error in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
    end if
    return
  end if
  
  ! check H
  call ZARACH2(N,N,H,INFO)
  if (INFO.NE.0) then
    ! print error in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"H is invalid",INFO,INFO)
    end if
    INFO = -3 
    return
  end if  
  
  ! check Z
  if (COMPZ.EQ.'V') then
    call ZARACH2(N,N,Z,INFO)
    if (INFO.NE.0) then
      ! print error in debug mode
      if (DEBUG) then
        call UARERR(__FILE__,__LINE__,"Z is invalid",INFO,INFO)
      end if
      INFO = -4 
      return
    end if 
  end if
  
  ! compress H
  call ZUHRFF(N,H,WORK(1:(3*N)),WORK((3*N+1):(5*N)),INFO)

  ! check info
  if (INFO.NE.0) then 
    ! print error in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"ZUHRFF failed",INFO,INFO)
    end if 
    INFO = 1
    return
  end if
  
  ! compute eigenvalues
  call ZUFFQR(COMPZ,N,WORK(1:(3*N)),WORK((3*N+1):(5*N)),Z,ITS,INFO) 

  ! check info
  if (INFO.NE.0) then 
    ! print error in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"ZUFFQR failed",INFO,INFO)
    end if
    INFO = 2 
    return
  end if
  
  ! update H
  H = cmplx(0d0,0d0,kind=8)
  do ii=1,N
    H(ii,ii) = cmplx(WORK((3*N)+2*(ii-1)+1),WORK((3*N)+2*(ii-1)+2),kind=8)
  end do
  
end subroutine ZUHFQR
