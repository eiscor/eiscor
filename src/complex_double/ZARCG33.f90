#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZARCG33 (Zomplex Auxiliary Routine Compute Givens generators 3->3)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generator for a Givens rotation represented
! by 3 real numbers: the real and imaginary parts of a complex cosine
! and a scrictly real sine. The first column is constructed to be 
! parallel with the complex vector [A,B]^T. The (3->3) refers to the 
! 3 double inputs and 3 double outputs (excluding the norm).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  AR, AI          REAL(8) 
!                    real and imaginary part of the first component 
!                    of the complex vector [A,B]^T
!
!  B               REAL(8) 
!                    the second component of the complex vector [A,B]^T
!
! OUTPUT VARIABLES:
!
!  CR, CI          REAL(8)
!                    on exit contAIns the generators for the cosine
!                    component of the Givens rotation
!
!  S               REAL(8)
!                    on exit contAIns the generators for the sine
!                    component of the Givens rotation
!
!  NRM             REAL(8)
!                    on exit contAIns the norm of the vector [A,B]^T
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies AR is invalid
!                    INFO = -2 implies AI is invalid
!                    INFO = -3 implies B is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZARCG33(AR,AI,B,CR,CI,S,NRM,INFO)

  implicit none
  
  ! input variables
  integer, intent(inout) :: INFO
  real(8), intent(in) :: AR, AI, B
  real(8), intent(inout) :: CR, CI, S, NRM
  
  ! compute variables
  real(8) :: tar, tai, tb
  real(8) :: nar, nai, nb
  real(8), parameter :: tol = epsilon(1d0)
  
  ! initialize INFO
  INFO = 0

  ! check error in debug mode
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

    ! check B
    call DARNAN(B,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"B is invalid",INFO,-3)
      return
    end if
    call DARINF(B,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"B is invalid",INFO,-3)
      return
    end if

  end if

  ! set local variables
  nar = abs(AR)
  nai = abs(AI)
  nb = abs(B)
  NRM = 1d0
  
  if(nar.EQ.0 .AND. nai.EQ.0 .AND. nb.EQ.0)then
    CR = 1d0
    CI = 0d0
    S = 0d0
    NRM = 0d0
  else if(nb.EQ.0 .AND. nar > nai)then
    tai = AI/AR
    NRM = sqrt(1d0 + tai*tai)
    if(AR < 0)then
      CR = -1d0/NRM
      CI = tai*CR
      S = 0d0
      NRM = -AR*NRM
    else
      CR = 1d0/NRM
      CI = tai*CR
      S = 0d0
      NRM = AR*NRM
    end if
  else if(nb.EQ.0)then
    tar = AR/AI
    NRM = sqrt(1d0 + tar*tar)
    if(AI < 0)then
      CI = -1d0/NRM
      CR = tar*CI
      S = 0d0
      NRM = -AI*NRM
    else
      CI = 1d0/NRM
      CR = tar*CI
      S = 0d0
      NRM = AI*NRM
    end if
  else if(nar >= nb .AND. nar >= nai)then
    tb = B/AR
    tai = AI/AR
    NRM = sqrt(1d0 + tb*tb + tai*tai)
    if(AR < 0)then
      CR = -1d0/NRM
      CI = tai*CR
      S = tb*CR
      NRM = -AR*NRM
    else
      CR = 1d0/NRM
      CI = tai*CR
      S = tb*CR
      NRM = AR*NRM
    end if
  else if(nai >= nb .AND. nai >= nar)then
    tb = B/AI
    tar = AR/AI
    NRM = sqrt(1d0 + tb*tb + tar*tar)
    if(AI < 0)then
      CI = -1d0/NRM
      CR = tar*CI
      S = tb*CI
      NRM = -AI*NRM
    else
      CI = 1d0/NRM
      CR = tar*CI
      S = tb*CI
      NRM = AI*NRM
    end if
  else
    tar = AR/B
    tai = AI/B
    NRM = sqrt(1d0 + tai*tai + tar*tar)
    if(B < 0)then
      S = -1d0/NRM
      CR = tar*S
      CI = tai*S
      NRM = -B*NRM
    else
      S = 1d0/NRM
      CR = tar*S
      CI = tai*S
      NRM = B*NRM
    end if
  end if

end subroutine ZARCG33
