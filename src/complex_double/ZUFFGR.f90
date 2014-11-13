#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZUFFGR (Zouble Unitary hessenberg Factored Fuse Givens Rotation)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine fuses a Givens rotation represented by 3 real numbers: 
! the real and imaginary parts of a complex cosine and a scrictly real 
! sine into the top or bottom of a unitary hessenberg matrix that is 
! stored as a product of givens rotations and a complex diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  JOB             CHARACTER
!                    'T': fuses at the Top from the left
!                    'B': fuses at the Bottom from the right
!
!  N               INTEGER
!                    dimension of matrix
!
!  STR             INTEGER
!                    index of the top most givens rotation where 
!                    fusion could happen
!
!  STP             INTEGER
!                    index of the bottom most givens rotation where 
!                    fusion could happen
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  B               REAL(8) array of dimension (3)
!                    generators for rotation that will be fused
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies JOB is invalid
!                   INFO = -2 implies N is invalid
!                   INFO = -3 implies STR is invalid
!                   INFO = -4 implies STP is invalid
!                   INFO = -7 implies B is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZUFFGR(JOB,N,STR,STP,Q,D,B,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: JOB
  integer, intent(in) :: N, STR, STP
  integer, intent(inout) :: INFO
  real(8), intent(inout) :: Q(3*(N-1)),D(2*N)
  real(8), intent(in) :: B(3)
  
  ! compute variables
  integer :: k,ii
  real(8) :: c1r, c1i, s1
  real(8) :: c2r, c2i, s2
  real(8) :: c3r, c3i, s3r, s3i
  real(8) :: d1r, d1i
  real(8) :: d2r, d2i
  real(8) :: phr, phi
  real(8) :: nrm
  
  ! initialize INFO
  INFO = 0

  ! check input in debug mode
  if (DEBUG) then
  
    ! check JOB
    if ((JOB.NE.'T').AND.(JOB.NE.'B')) then
      INFO = -1
      call UARERR(__FILE__,__LINE__,"JOB must be 'T' or 'B'",INFO,INFO)
      return
    end if
    
    ! check N
    call IARNAN(N,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"N is invalid",INFO,-2)
      return
    end if
    call IARINF(N,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"N is invalid",INFO,-2)
      return
    end if
    if (N < 2) then
      INFO = -2
      call UARERR(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
      return
    end if
    
    ! check STR
    if ((STR < 1).OR.(STR > N-1)) then
      INFO = -3
      call UARERR(__FILE__,__LINE__,"STR must 1 <= STR <= N-1",INFO,INFO)
      return
    end if 
    
    ! check STP
    if ((STP < STR).OR.(STP > N-1)) then
      INFO = -4
      call UARERR(__FILE__,__LINE__,"STP must STR <= STP <= N-1",INFO,INFO)
      return
    end if  
    
    ! check B
    call DARACH1(3,B,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"B is invalid",INFO,-7)
      return
    end if   

  end if
  
  ! fusion at bottom
  if(JOB.EQ.'B')then
  
    ! pass through diag
    call ZARGTD('R',D((2*(STP-1)+1):(2*(STP-1)+4)),B,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"ZARGTD failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! set inputs  
    c2r = B(1)
    c2i = B(2)
    s2 = B(3)
  
    ! retrieve Q  
    k = 3*(STP-1)
    c1r = Q(k+1)
    c1i = Q(k+2)
    s1 = Q(k+3)
     
    ! retrieve D
    k = 2*(STP-1)
    d1r = D(k+1)
    d1i = D(k+2)
    d2r = D(k+3)
    d2i = D(k+4)
     
    ! compute givens product
    c3r = c1r*c2r - c1i*c2i - s1*s2
    c3i = c1r*c2i + c1i*c2r
    s3r = s1*c2r + s2*c1r
    s3i = s1*c2i - s2*c1i
    
    ! compute phase
    nrm = abs(cmplx(s3r,s3i,kind=8))
    if(nrm /= 0)then
      phr = s3r/nrm
      phi = s3i/nrm
    else
      phr = 1d0
      phi = 0d0
    end if
     
    ! update Q
    c2r = c3r*phr + c3i*phi
    c2i = -c3r*phi + c3i*phr
    s2 = nrm
    nrm = sqrt(c2r*c2r + c2i*c2i + s2*s2)
    c2r = c2r/nrm
    c2i = c2i/nrm
    s2 = s2/nrm
    k = 3*(STP-1)
    Q(k+1) = c2r
    Q(k+2) = c2i
    Q(k+3) = s2
     
    ! update D
    c1r = phr*d1r - phi*d1i
    c1i = phr*d1i + phi*d1r
    nrm = sqrt(c1r*c1r + c1i*c1i)
    c1r = c1r/nrm
    c1i = c1i/nrm
    c2r = phr*d2r + phi*d2i
    c2i = phr*d2i - phi*d2r
    nrm = sqrt(c2r*c2r + c2i*c2i)
    c2r = c2r/nrm
    c2i = c2i/nrm
    k = 2*(STP-1)
    D(k+1) = c1r 
    D(k+2) = c1i
    D(k+3) = c2r
    D(k+4) = c2i
     
  ! fusion at the top
  else
  
    ! set inputs  
    c2r = B(1)
    c2i = B(2)
    s2 = B(3)
    
    ! retrieve Q  
    k = 3*(STR-1)
    c1r = Q(k+1)
    c1i = Q(k+2)
    s1 = Q(k+3)
     
    ! retrieve D
    k = 2*(STR-1)
    d1r = D(k+1)
    d1i = D(k+2)
    k = 2*(STP)
    d2r = D(k+1)
    d2i = D(k+2)
     
    ! compute givens product
    c3r = c1r*c2r - c1i*c2i - s1*s2
    c3i = c1r*c2i + c1i*c2r
    s3r = s1*c2r + s2*c1r
    s3i = -(s1*c2i - s2*c1i)
     
    ! compute phase
    nrm = abs(cmplx(s3r,s3i,kind=8))
    if(nrm /= 0)then
       phr = s3r/nrm
       phi = s3i/nrm
    else
       phr = 1d0
       phi = 0d0
    end if
     
    ! update Q
    c2r = c3r*phr + c3i*phi
    c2i = -c3r*phi + c3i*phr
    s2 = nrm
    nrm = sqrt(c2r*c2r + c2i*c2i + s2*s2)
    c2r = c2r/nrm
    c2i = c2i/nrm
    s2 = s2/nrm
    k = 3*(STR-1)
    Q(k+1) = c2r
    Q(k+2) = c2i
    Q(k+3) = s2
     
    do ii=(STR+1),STP
      k = 3*(ii-1)
      c1r = Q(k+1)
      c1i = Q(k+2)
      s1 = Q(k+3)
      nrm = c1r*phr + c1i*phi
      c1i = -c1r*phi + c1i*phr
      c1r = nrm
      nrm = sqrt(c1r*c1r + c1i*c1i + s1*s1)
      c1r = c1r/nrm
      c1i = c1i/nrm
      s1 = s1/nrm
      Q(k+1) = c1r
      Q(k+2) = c1i
      Q(k+3) = s1
    end do
     
    ! update D
    c1r = phr*d1r - phi*d1i
    c1i = phr*d1i + phi*d1r
    nrm = sqrt(c1r*c1r + c1i*c1i)
    c1r = c1r/nrm
    c1i = c1i/nrm
    c2r = phr*d2r + phi*d2i
    c2i = phr*d2i - phi*d2r
    nrm = sqrt(c2r*c2r + c2i*c2i)
    c2r = c2r/nrm
    c2i = c2i/nrm
    
    k = 2*(STR-1)
    D(k+1) = c1r 
    D(k+2) = c1i
    k = 2*(STP)
    D(k+1) = c2r
    D(k+2) = c2i
     
  end if

end subroutine ZUFFGR
