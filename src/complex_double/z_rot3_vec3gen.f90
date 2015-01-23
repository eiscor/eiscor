#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot3_vec3gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generators for a complex rotation represented
! by 3 real numbers: the real and imaginary parts of a complex cosine,
! CR and CI, and a scrictly real sine, S. The CR, CI and S are 
! constructed to be parallel with the vector [AR+iAI,B]^T, where AR, 
! AI and B are real and i = sqrt(-1).
!
! If any part of the input vector [AR+iAI,B]^T contains a NAN then
! CR, CI, S and NRM are all set to NAN.
!
! If only one of AR, AI or B = +/-INF then the corresponding AR, AI or B is 
! first set to +/-1 and the remaining terms are set to 0. Then the CR, CI 
! and S are computed from the new vector containing +/-1 and 0. In 
! this case NRM is always set to INF.
!
! If more than one of AR, AI or B = +/- INF then CR = CI = S = NRM = NAN
!
! If AR = AI = B = 0 then CR = 1, CI = S = 0 and NRM = 0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! EXCEPTIONAL CASES
!
!    AR   |    AI   |    B    |    CR   |    CI   |    S    |   NRM 
! ------- | ------- | ------- | ------- | ------- | ------- | -------
!   0d0   |   0d0   |   0d0   |   1d0   |   0d0   |   0d0   |   0d0
! ------- | ------- | ------- | ------- | ------- | ------- | -------
! +-INF   |   XdX   |   XdX   | +-1d0   |   0d0   |   0d0   |   INF
!   XdX   | +-INF   |   XdX   |   0d0   | +-1d0   |   0d0   |   INF
!   XdX   |   XdX   | +-INF   |   0d0   |   0d0   | +-1d0   |   INF
! +-INF   | +-INF   |   XdX   |   NAN   |   NAN   |   NAN   |   NAN   
! +-INF   |   XdX   | +-INF   |   NAN   |   NAN   |   NAN   |   NAN   
!   XdX   | +-INF   | +-INF   |   NAN   |   NAN   |   NAN   |   NAN   
! +-INF   | +-INF   | +-INF   |   NAN   |   NAN   |   NAN   |   NAN   
! ------- | ------- | ------- | ------- | ------- | ------- | -------
!   NAN   |   XdX   |   XdX   |   NAN   |   NAN   |   NAN   |   NAN
!   XdX   |   NAN   |   XdX   |   NAN   |   NAN   |   NAN   |   NAN
!   XdX   |   XdX   |   NAN   |   NAN   |   NAN   |   NAN   |   NAN
!   NAN   |   NAN   |   XdX   |   NAN   |   NAN   |   NAN   |   NAN
!   XdX   |   NAN   |   NAN   |   NAN   |   NAN   |   NAN   |   NAN
!   NAN   |   XdX   |   NAN   |   NAN   |   NAN   |   NAN   |   NAN
!   NAN   |   NAN   |   NAN   |   NAN   |   NAN   |   NAN   |   NAN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  AR, AI          REAL(8) 
!                    real and imaginary part of the first component 
!                    of the complex vector [AR+iAI,B]^T
!
!  B               REAL(8) 
!                    the second component 
!                    of the complex vector [AR+iAI,B]^T
!
! OUTPUT VARIABLES:
!
!  CR, CI          REAL(8)
!                    on exit CR = AR/NRM, CI = AI/NRM
!
!  S               REAL(8)
!                    on exit S = B/NRM
!
!  NRM             REAL(8)
!                    on exit contains the norm 
!                    of the complex vector [AR+iAI,B]^T
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot3_vec3gen(AR,AI,B,CR,CI,S,NRM)

  implicit none
  
  ! input variables
  real(8), intent(in) :: AR, AI, B
  real(8), intent(inout) :: CR, CI, S, NRM
  
  ! compute variables
  real(8) :: tar, tai, tb
  real(8) :: nar, nai, nb

  ! set local variables
  nar = abs(AR)
  nai = abs(AI)
  nb = abs(B)
  NRM = 1d0
  
  ! 2 or more INFs
  if (((nar>EISCOR_DBL_INF).AND.((nai>EISCOR_DBL_INF).OR.(nb>EISCOR_DBL_INF))).OR.((nai>EISCOR_DBL_INF).AND.(nb>EISCOR_DBL_INF))) then
  
    CR = 0d0
    CR = 0d0/CR
    CI = CR
    S = CR
    NRM = CR
        
    return
  
  end if
  
  ! AR = AI = B = 0
  if(nar.EQ.0 .AND. nai.EQ.0 .AND. nb.EQ.0)then
  
    CR = 1d0
    CI = 0d0
    S = 0d0
    NRM = 0d0
  
  ! B = 0 and |AR| > |AI|
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
    
  ! B = 0 and |AR| <= |AI|    
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
    
  ! B /= 0, |AR| >= |B| and |AR| >= |AI| 
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
    
  ! B /= 0, |AI| >= |B| and |AI| >= |AR| 
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
    
  ! B /= 0, |B| >= |AR| and |B| >= |AI|     
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

end subroutine z_rot3_vec3gen
