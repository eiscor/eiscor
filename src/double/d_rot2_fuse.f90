#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_rot2_fuse
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the product of two real Given's rotations and 
! and stores the output in one of the input arrays.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  JOB             CHARACTER
!                    'L': store fusion in Q1
!                    'R': store fusion in Q2
!
!  Q1, Q2          REAL(8) array of dimension (2)
!                    generators for two Givens rotations
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = 1 implies failure in d_rot2_vec2gen
!                    INFO = -1 implies JOB is invalid
!                    INFO = -2 implies Q1 is invalid
!                    INFO = -3 implies Q2 is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_rot2_fuse(JOB,Q1,Q2,INFO)

  implicit none

  ! input variables
  character, intent(in) :: JOB
  integer, intent(inout) :: INFO
  real(8), intent(inout) :: Q1(2), Q2(2)

  ! compute variables
  real(8), parameter :: tol = epsilon(1d0)
  real(8) :: nrm

  ! initialize INFO
  INFO = 0

  ! check input in debug mode
  if (DEBUG) then 
  
    ! check JOB
    if ((JOB.NE.'L').AND.(JOB.NE.'R')) then
      INFO = -1
      call u_infocode_check(__FILE__,__LINE__,"JOB is invalid",INFO,INFO)
      return
    end if
    
    ! check Q1
    call d_1Darray_check(2,Q1,INFO)
    if (INFO.NE.0) then 
      call u_infocode_check(__FILE__,__LINE__,"Q1 is invalid",INFO,-2)
      return 
    end if 
    nrm = sqrt(Q1(1)**2 + Q1(2)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -2
      call u_infocode_check(__FILE__,__LINE__,"Q1 is not orthogonal to working precision",INFO,INFO)
      return
    end if
    
    ! check Q2
    call d_1Darray_check(2,Q2,INFO)
    if (INFO.NE.0) then 
      call u_infocode_check(__FILE__,__LINE__,"Q2 is invalid",INFO,-3)
      return 
    end if 
    nrm = sqrt(Q2(1)**2 + Q2(2)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -3
      call u_infocode_check(__FILE__,__LINE__,"Q2 is not orthogonal to working precision",INFO,INFO)
      return
    end if
  
  end if

  ! store in Q2
  if(JOB.EQ.'R')then
    ! compute new generators
    nrm = Q1(1)*Q2(1) - Q1(2)*Q2(2)
    Q1(2) = Q1(2)*Q2(1) + Q1(1)*Q2(2)
    Q1(1) = nrm
    
    ! compute new generators
    call d_rot2_vec2gen(Q1(1),Q1(2),Q2(1),Q2(2),nrm,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,1)
      if (INFO.NE.0) then 
        return 
      end if 
    end if

  ! store in Q2
  else
    ! compute new generators
    nrm = Q1(1)*Q2(1) - Q1(2)*Q2(2)
    Q2(2) = Q1(2)*Q2(1) + Q1(1)*Q2(2)
    Q2(1) = nrm
    
    ! compute new generators
    call d_rot2_vec2gen(Q2(1),Q2(2),Q1(1),Q1(2),nrm,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,1)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
  end if

end subroutine d_rot2_fuse
