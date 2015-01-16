#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthhess_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the eigenvalues and optionally eigenvectors
! of a real orthogonal upper hessenberg matrix.
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
!  H               Real(8) array of dimension (N,N)
!                    orthogonal hessenberg matrix, assumed that H(ii,jj) = 0 for |ii-jj| > 0
!                    on exit contains a block diagonal matrix with blocks no greater than 2x2. 
!
!  WORK            REAL(8) array of dimension (4*N)
!                    work space for eigensolver
!
! OUTPUT VARIABLES:
!
!  Z              REAL(8) array of dimension (N,N)
!                    if COMPZ = 'N' unused
!                    if COMPZ = 'I' stores eigenvectors, initializing to identity 
!                    if COMPZ = 'V' stores eigenvectors, assume Z already initialized
!
!  ITS            INTEGER array of dimension (N-1)
!                    Contains the number of iterations per deflation
!
!  INFO           INTEGER
!                    INFO = 2 implies d_orthfact_qr failed
!                    INFO = 1 implies d_orthhess_factor failed
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies COMPZ is invalid
!                    INFO = -2 implies N is invalid
!                    INFO = -3 implies H is invalid
!                    INFO = -4 implies Z is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthhess_qr(COMPZ,N,H,Z,ITS,WORK,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: COMPZ
  integer, intent(in) :: N
  real(8), intent(inout) :: WORK(4*N)
  integer, intent(inout) :: ITS(N-1), INFO
  real(8), intent(inout) :: H(N,N), Z(N,N)
  
  ! compute variables
  integer :: ii,cpair
  real(8) :: c,s
  
  ! initialize INFO
  INFO = 0
  
  ! check COMPZ
  if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
    INFO = -1 
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"COMPZ must be 'N', 'I' or 'V'",INFO,INFO)
    end if   
    return
  end if
  
  ! check N
  if (N < 2) then
    INFO = -2
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
    end if  
    return
  end if
  
  ! check H
  call d_2Darray_check(N,N,H,INFO)
  if (INFO.NE.0) then
    INFO = -3
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"H is invalid",INFO,INFO)
    end if  
    return
  end if  
  
  ! check Z
  if (COMPZ.EQ.'V') then
    call d_2Darray_check(N,N,Z,INFO)
    if (INFO.NE.0) then
      INFO = -4
      ! print error in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"Z is invalid",INFO,INFO)
      end if  
      return
    end if 
  end if
  
  ! compress H
  call d_orthhess_factor(N,H,WORK(1:(2*N)),WORK((2*N+1):(4*N)),INFO)
  if (INFO.NE.0) then
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_orthhess_factor",INFO,1)
    end if  
    return
  end if  
  
  ! compute eigenvalues
  print*, N
  call d_orthfact_qr(COMPZ,N,WORK(1:(2*N)),WORK((2*N+1):(4*N)),Z,ITS,INFO)
  if (INFO.NE.0) then
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_orthfact_qr failed",INFO,2)
    end if  
    return
  end if 
  
  ! update H
  cpair = 0
  H = 0d0
  do ii=1,N
    c = WORK((2*N)+2*(ii-1)+1)
    s = WORK((2*N)+2*(ii-1)+2)
    if (abs(s).EQ.0d0) then
      H(ii,ii) = c
    else if (cpair.EQ.0) then
      H(ii,ii) = c
      H(ii,ii+1) = s
      H(ii+1,ii) = -s
      H(ii+1,ii+1) = c
      cpair = 1
    else
      cpair = 0
    end if
  end do
  
end subroutine d_orthhess_qr
