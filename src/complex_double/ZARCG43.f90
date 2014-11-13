#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZARCG43 (Zomplex Auxiliary Routine Compute Givens generators 4->3)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generator for a Givens rotation represented
! by 3 real numbers: the real and imaginary parts of a complex cosine
! and a scrictly real sine. The first column is constructed to be 
! parallel with the complex vector [A,B]^T. The (4->3) refers to the 
! 4 double inputs and 3 double outputs (excluding the norm).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  AR, AI          REAL(8) 
!                    real and imaginary part of the first component 
!                    of the complex vector [A,B]^T
!
!  BR, BI          REAL(8) 
!                    real and imaginary part of the second component 
!                    of the complex vector [A,B]^T
!
! OUTPUT VARIABLES:
!
!  CR, CI          REAL(8)
!                    on exit contains the generators for the cosine
!                    component of the Givens rotation
!
!  S               REAL(8)
!                    on exit contains the generators for the sine
!                    component of the Givens rotation
!
!  NRM             REAL(8)
!                    on exit contains the norm of the vector [A,B]^T
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies AR is invalid
!                    INFO = -2 implies AI is invalid
!                    INFO = -3 implies BR is invalid
!                    INFO = -4 implies BI is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZARCG43(AR,AI,BR,BI,CR,CI,S,NRM,INFO)

  implicit none
  
  ! input variables
  integer, intent(inout) :: INFO
  real(8), intent(in) :: AR, AI, BR, BI
  real(8), intent(inout) :: CR, CI, S, NRM
  
  ! compute variables
  real(8) :: sr, si, temp
  
  ! initialize INFO
  INFO = 0
  
  ! print error in debug mode
  if (DEBUG) then
  
    ! check AR
    call DARNAN(AR,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"AR is invalid",INFO,-1)
      return
    end if
    call DARINF(AR,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"AR is invalid",INFO,-1)
      return
    end if   

    ! check AI
    call DARNAN(AI,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"AI is invalid",INFO,-2)
      return
    end if
    call DARINF(AI,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"AI is invalid",INFO,-2)
      return
    end if  

    ! check BR
    call DARNAN(BR,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"BR is invalid",INFO,-3)
      return
    end if
    call DARINF(BR,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"BR is invalid",INFO,-3)
      return
    end if
    
    ! check BI
    call DARNAN(BI,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"BI is invalid",INFO,-4)
      return
    end if
    call DARINF(BI,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"BI is invalid",INFO,-4)
      return
    end if

  end if 
  
  ! compute NRM
  CR = abs(AR)
  CI = abs(AI)
  sr = abs(BR)
  si = abs(BI)
  
  ! exit if nrm is zero
  if ((CR.EQ.0).AND.(CI.EQ.0).AND.(sr.EQ.0).AND.(si.EQ.0)) then
    nrm = 0d0
    CR = 1d0
    CI = 0d0
    S = 0d0
    return
  end if
  
  ! find maximum value
  NRM = CR
  if (CI > NRM) then
    NRM = CI
  end if
  if (sr > NRM) then
    NRM = sr
  end if
  if (si > NRM) then
    NRM = si
  end if
  
  ! scale
  CR = CR/NRM
  CI = CI/NRM
  sr = sr/NRM
  si = si/NRM
  
  ! NRM
  NRM = NRM*sqrt(CR*CR + CI*CI + sr*sr + si*si)
  
  ! update Generators
  CR = AR/NRM
  CI = AI/NRM
  sr = BR/NRM
  si = BI/NRM
  
  ! check si
  if (abs(si) .EQ. 0d0) then
    S = sr
    return
  end if
  
  ! compute phase
  if (abs(sr) > abs(si)) then
    S = abs(sr)*sqrt(1d0+abs(si/sr)**2)
  else
    S = abs(si)*sqrt(1d0+abs(sr/si)**2)
  end if
  sr = sr/S
  si = si/S
  
  ! update CR and CI
  temp = CR*sr + CI*si
  CI = -CR*si + CI*sr
  CR = temp  

end subroutine ZARCG43
