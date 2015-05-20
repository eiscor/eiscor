#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_rot4_vec4gen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The description has to be written
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generators for a complex rotation represented
! by 4 real numbers: the real and imaginary parts of a complex cosine,
! CR and CI, and a the real and complex parts of a complex sine, SR and SI.
! The CR, CI, SR, and SI are constructed to be parallel with the 
! vector [AR+iAI,BR+iBI]^T, where AR, 
! AI, BR and BI are real and i = sqrt(-1).
!
! If any part of the input vector [AR+iAI,BR+iBI]^T contains a NAN then
! CR, CI, S and NRM are all set to NAN.
!
! If only one of AR, AI, BR or BI = +/-INF then the corresponding AR, AI, BR 
! or BI is first set to +/-1 and the remaining terms are set to 0. Then the 
! CR, CI and S are computed from the new vector containing +/-1 and 0. In 
! this case NRM is always set to INF.
!
! If more than one of AR, AI, BR or BI = +/- INF then CR = CI = S = NRM = NAN
!
! If AR = AI = BR = BI = 0 then CR = 1, CI = S = 0 and NRM = 0.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! EXCEPTIONAL CASES
!
!    AR   |    AI   |    BR   |    BI   |    CR   |    CI   |    SR   |    SI   |   NRM 
! ------- | ------- | ------- | ------- | ------- | ------- | ------- | ------- | -------
!   0d0   |   0d0   |   0d0   |   0d0   |   1d0   |   0d0   |   0d0   |   0d0
! ------- | ------- | ------- | ------- | ------- | ------- | ------- | ------- | -------
! +-INF   |   XdX   |   XdX   |   XdX   | +-1d0   |   0d0   |   0d0   |   INF
!   XdX   | +-INF   |   XdX   |   XdX   |   0d0   | +-1d0   |   0d0   |   INF
!   XdX   |   XdX   | +-INF   |   XdX   |   0d0   |   0d0   |   1d0   |   INF
!   XdX   |   XdX   |   XdX   | +-INF   |   0d0   |   0d0   |   1d0   |   INF
!          at least two +-INFs          |   NAN   |   NAN   |   NAN   |   NAN
! ------- | ------- | ------- | ------- | ------- | ------- | ------- | -------
!           at least one NAN            |   NAN   |   NAN   |   NAN   |   NAN
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  AR, AI          REAL(8) 
!                    real and imaginary part of the first component 
!                    of the complex vector [AR+iAI,BR+iBI]^T
!
!  BR, BI          REAL(8) 
!                    the second component 
!                    of the complex vector [AR+iAI,BR+iBI]^T
!
! OUTPUT VARIABLES:
!
!  CR, CI          REAL(8)
!                    on exit CR = AR/NRM, CI = AI/NRM
!
!  SR, SI          REAL(8)
!                    on exit SR = BR/NRM, SI = BI/NRM
!
!  NRM             REAL(8)
!                    on exit contains the norm 
!                    of the complex vector [AR+iAI,BR+iBI]^T
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_rot4_vec4gen(AR,AI,BR,BI,CR,CI,SR,SI,NRM)

  implicit none
  
  ! input variables
  real(8), intent(in) :: AR, AI, BR, BI
  real(8), intent(inout) :: CR, CI, SR, SI, NRM
  
  ! compute variables
  real(8) :: tar, tai, tbr, tbi
  real(8) :: nar, nai, nbr, nbi
  ! compute variables
  real(8), parameter :: inf = EISCOR_DBL_INF
  

  ! set local variables
  nar = abs(AR)
  nai = abs(AI)
  nbr = abs(BR)
  nbi = abs(BI)
  NRM = 1d0
  
  ! 2 or more INFs
  if (((nar>EISCOR_DBL_INF).AND.((nai>EISCOR_DBL_INF).OR.(nbr>EISCOR_DBL_INF).OR.(nbi>EISCOR_DBL_INF)))&
       &.OR.((nai>EISCOR_DBL_INF).AND.((nbr>EISCOR_DBL_INF).OR.(nbi>EISCOR_DBL_INF)))&
       &.OR.((nbr>EISCOR_DBL_INF).AND.(nbi>EISCOR_DBL_INF))) then
     
     CR = 0d0
     CR = 0d0/CR
     CI = CR
     SR = CR
     SI = CR
     NRM = CR
     
     return
     
  end if
  
  ! AR = AI = BR = BI = 0
  if(nar.EQ.0 .AND. nai.EQ.0 .AND. nbr.EQ.0 .AND. nbi.EQ.0)then
     
     CR = 1d0
     CI = 0d0
     SR = 0d0
     SI = 0d0
     NRM = 0d0
     
!!$  ! AR = AI = BR = 0, BI /=0
!!$  else
!!$     CR = 0d0
!!$     CI = 0d0
!!$     SR = 0d0
!!$     if (BI.LT.0d0) then
!!$        SI = -1d0
!!$     else
!!$        SI = 1d0
!!$     end if
!!$     NRM = SI*BI
     
  ! BR = BI = 0 
  else if(nbr.EQ.0 .AND. nbi.EQ.0)then
  
     call d_rot2_vec2gen(AR,AI,CR,CI,NRM)
     SR = 0d0
     SI = 0d0

   
  ! AR,AI,BR,BI /= 0, |AR| largest
  else if(nar >= nai .AND. nar >= nbr .AND. nar >= nbi)then
  
     tai = AI/AR
     tbr = BR/AR
     tbi = BI/AR
     
     NRM = sqrt(1d0 + tai*tai + tbr*tbr + tbi*tbi)
     if(AR < 0)then
        CR = -1d0/NRM
        CI = tai*CR
        SR = tbr*CR
        SI = tbi*CR
        NRM = -AR*NRM
     else
        CR = 1d0/NRM
        CI = tai*CR
        SR = tbr*CR
        SI = tbi*CR
        NRM = AR*NRM
     end if

     
     ! AR,AI,BR,BI /= 0, |AI| largest
  else if(nai >= nbr .AND. nai >= nbi)then
     
     tar = AR/AI
     tbr = BR/AI
     tbi = BI/AI
     
     NRM = sqrt(1d0 + tar*tar + tbr*tbr + tbi*tbi)
     if(AR < 0)then
        CI = -1d0/NRM
        CR = tar*CI
        SR = tbr*CI
        SI = tbi*CI
        NRM = -AI*NRM
     else
        CI = 1d0/NRM
        CR = tar*CI
        SR = tbr*CI
        SI = tbi*CI
        NRM = AI*NRM
    end if

    
     ! AR,AI,BR,BI /= 0, |BR| largest
  else if(nbr >= nbi)then
     
     tar = AR/BR
     tai = AI/BR
     tbi = BI/BR
     
     NRM = sqrt(1d0 + tar*tar + tai*tai + tbi*tbi)
     if(AR < 0)then
        SR = -1d0/NRM
        CR = tar*SR
        CI = tai*SR
        SI = tbi*SR
        NRM = -BR*NRM
     else
        SR = 1d0/NRM
        CR = tar*SR
        CI = tai*SR
        SI = tbi*SR
        NRM = BR*NRM
    end if


     ! AR,AI,BR,BI /= 0, |BI| largest
  else 
     
     tar = AR/BI
     tai = AI/BI
     tbr = BR/BI
     
     NRM = sqrt(1d0 + tar*tar + tai*tai + tbr*tbr)
     if(AR < 0)then
        SI = -1d0/NRM
        CR = tar*SI
        CI = tai*SI
        SR = tbi*SI
        NRM = -BI*NRM
     else
        SI = 1d0/NRM
        CR = tar*SI
        CI = tai*SI
        SR = tbr*SI
        NRM = BI*NRM
    end if

 end if

end subroutine z_rot4_vec4gen
