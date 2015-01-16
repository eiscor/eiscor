#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_buildbulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the first transformation in a single shift
! iteration for a unitary upper hessenberg matix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  K               INTEGER
!                    index of the diagonal block
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!                    generators must be orthogonal to working precision
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
!  SHFT            COMPLEX(8) 
!                    contains the shift need for the first transformation
!
! OUTPUT VARIABLES:
!
!  B               REAL(8) array of dimension (3)
!                    on exit contains the generators for the first
!                    transformation
!
!  INFO            INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies N is invalid
!                   INFO = -2 implies K is invalid
!                   INFO = -3 implies Q is invalid
!                   INFO = -4 implies D is invalid
!                   INFO = -5 implies SHFT is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_buildbulge(N,K,Q,D,SHFT,B,INFO)
  
  implicit none
  
  ! input variables
  integer, intent(in) :: N, K
  integer, intent(inout) :: INFO
  real(8), intent(in) :: Q(3*N), D(2*N)
  complex(8), intent(in) :: SHFT
  real(8), intent(inout) :: B(3)
  
  ! compute variables
  real(8) :: nrm
  complex(8) :: block(2,2)
  
  ! initialize info
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
    
    ! check N, Q and D
    call z_unifact_factorcheck(N,Q,D,INFO)
    if (INFO.EQ.-1) then
      INFO = -1
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
      return
    end if
    if (INFO.EQ.-2) then
      INFO = -5
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,INFO)
      return
    end if
    if (INFO.EQ.-3) then
      INFO = -6
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,INFO)
      return
    end if
    
    ! check K
    if ((K < 1).OR.(K > N-1)) then
      INFO = -2
      call u_infocode_check(__FILE__,__LINE__,"K must be 1 <= K <= N-1",INFO,INFO)
      return
    end if 
  
    ! check SHFT
    call z_scalar_nancheck(SHFT,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"SHFT is invalid",INFO,-5)
      return
    end if
    call z_scalar_infcheck(SHFT,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"SHFT is invalid",INFO,-5)
      return
    end if
  
  end if
  
  ! get top block
  call z_unifact_2x2diagblock(N,K,Q,D,block,INFO) 
      
  ! check INFO in debug mode
  if (DEBUG) then
    call u_infocode_check(__FILE__,__LINE__,"z_unifact_2x2diagblock failed",INFO,INFO)
    if (INFO.NE.0) then 
      return 
    end if 
  end if
  
  ! shift first entry
  block(1,1) = block(1,1) - SHFT
  
  ! bulge
  call z_rot3_vec4gen(dble(block(1,1)),aimag(block(1,1)),dble(block(2,1)), &
    aimag(block(2,1)),B(1),B(2),B(3),nrm,INFO)
      
  ! check INFO in debug mode
  if (DEBUG) then
    call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec4gen failed",INFO,INFO)
    if (INFO.NE.0) then 
      return 
    end if 
  end if

end subroutine z_unifact_buildbulge
