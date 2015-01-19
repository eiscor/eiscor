#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthhess_factor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine factors a real orthogonal upper hessenberg matrix into
! a sequence of Givens' rotations and a diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  H               REAL(8) array of dimension (N,N)
!                    orthogonal hessenberg matrix to be reduced
!
! OUTPUT VARIABLES:
!
!  Q               REAL(8) array of dimension (2*(N-1))
!                    array of generators for Givens' rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  INFO           INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = 1 implies d_rot2_vec2gen failed
!                    INFO = -1 implies N is invalid
!                    INFO = -2 implies H is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthhess_factor(N,H,Q,D,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(inout) :: H(N,N), Q(2*(N-1)), D(2*N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii, ind
  real(8) :: nrm, tol 
  real(8) :: c, s
  
  ! initialize info
  INFO = 0
  
  ! check N
  if (N < 2) then
    INFO = -1
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
    end if  
    return
  end if
  
  ! check H
  call d_2Darray_check(N,N,H,INFO)
  if (INFO.NE.0) then
    INFO = -2
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"H is invalid",INFO,INFO)
    end if  
    return
  end if  
  
  ! set tol
  tol = dble(max(10,N))*epsilon(1d0)
  
  ! loop for reduction
  ! ii is the column being reduced
  do ii=1,(N-1)
            
    ! reduce to block diagonal
    call d_rot2_vec2gen(H(ii,ii),H(ii+1,ii),c,s,nrm,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,1)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
     
    ! check for unitarity       
    if (abs(nrm-1d0) >= tol) then
      INFO = -3
      ! print error in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"H is not orthogonal",INFO,INFO)
      end if  
      return          
    end if
    
    ! update H
    H(ii,ii) = c*H(ii,ii) + s*H(ii+1,ii)
    H(ii+1,(ii+1):N) = -s*H(ii,(ii+1):N) + c*H(ii+1,(ii+1):N)
         
    ! store in Q
    Q(2*(ii-1)+1) = c
    Q(2*(ii-1)+2) = s
          
    ! store in D
    D(2*(ii-1)+1) = sign(1d0,H(ii,ii))
    D(2*(ii-1)+2) = 0d0
   
  end do
  
  ! check for unitarity       
  if (abs(abs(H(N,N))-1d0) >= tol) then
    INFO = -3
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"H is not orthogonal",INFO,INFO)
    end if  
    return        
  end if
          
  ! store in D
  D(2*(N-1)+1) = sign(1d0,H(N,N))
  D(2*(N-1)+2) = 0d0

end subroutine d_orthhess_factor
