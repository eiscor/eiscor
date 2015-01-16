#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthfact_buildbulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the first Givens' rotations that start a real 
! Francis iteration. If JOB = 'S' a single Givens' rotation is 
! constructed. If JOB = 'D' two Givens' rotations are constructed for
! a real double shift step.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  JOB             CHARACTER
!                    'S': real singleshift
!                    'D': real double shift
!
!  N               INTEGER
!                    dimension of matrix
!
!  STR             INTEGER
!                    index of the top most givens rotation where 
!                    the iteration begins
!
!  Q               REAL(8) array of dimension (2*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  E               REAL(8) array of dimension (2)
!                    real and imaginary part of complex shift
!
! OUTPUT VARIABLES:
!
!  B1, B2          REAL(8) array of dimension (2)
!                    generators for givens rotations that represent
!                    the first transformation
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies JOB is invalid
!                    INFO = -2 implies N is invalid
!                    INFO = -3 implies STR is invalid
!                    INFO = -4 implies Q is invalid
!                    INFO = -5 implies D is invalid
!                    INFO = -6 implies E is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthfact_buildbulge(JOB,N,STR,Q,D,E,B1,B2,INFO)
  
  implicit none
  
  ! input variables
  character, intent(in) :: JOB
  integer, intent(in) :: N, STR
  integer, intent(inout) :: INFO
  real(8), intent(in) :: Q(2*(N-1)), D(2*N), E(2)
  real(8), intent(inout) :: B1(2), B2(2)
  
  ! compute variables
  real(8) :: rrho1, rrho2, irho1, irho2 
  real(8) :: nrm1, nrm2, COL(3), T(3,2)
  
  ! initialize INFO
  INFO = 0

  ! check input in debug mode
  if (DEBUG) then 
  
    ! check JOB
    if ((JOB.NE.'S').AND.(JOB.NE.'D')) then
      INFO = -1
      call u_infocode_check(__FILE__,__LINE__,"JOB is invalid",INFO,INFO)
      return
    end if  

    ! check N
    if (N < 2) then
      INFO = -2
      call u_infocode_check(__FILE__,__LINE__,"N is invalid",INFO,INFO)
      return
    end if         
  
    ! check STR
    if ((STR >= N).OR.(STR < 1)) then
      INFO = -3
      call u_infocode_check(__FILE__,__LINE__,"STR is invalid",INFO,INFO)
      return
    end if  

    ! check Q and D
    call d_orthfact_factorcheck(N,Q,D,INFO)
    if (INFO.EQ.-1) then
      call u_infocode_check(__FILE__,__LINE__,"N is invalid",INFO,-1)
      return
    end if
    if (INFO.EQ.-2) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,-4)
      return
    end if
    if (INFO.EQ.-3) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,-5)
      return
    end if
  
    ! check E
    call d_1Darray_check(2,E,INFO)
    call u_infocode_check(__FILE__,__LINE__,"E is invalid",INFO,-6)
    if (INFO.NE.0) then 
      return 
    end if 
    
  end if
  
  ! single shift
  if (JOB.EQ.'S') then
  
    ! check E in debug mode
    if (DEBUG) then
      if (abs(E(2)).NE.0d0) then
        INFO = -6
        call u_infocode_check(__FILE__,__LINE__,"E should be strictly real",INFO,INFO)
        return
      end if
    end if 
  
    ! compute first block of A
    call d_orthfact_2x2diagblock(N,STR,Q,D,T(1:2,1:2),INFO)  

    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_orthfact_2x2diagblock failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! compute first Givens' rotation
    call d_rot2_vec2gen(T(1,1)-E(1),T(2,1),B1(1),B1(2),nrm1,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if 
  
  ! double shift
  else
  
    ! set shifts
    rrho1 = E(1)
    irho1 = E(2)
    rrho2 = E(1)
    irho2 = -E(2) 
    
    ! compute first two columns of A
    call d_orthfact_2x2diagblock(N,STR+1,Q,D,T(1:2,1:2),INFO)  

    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_orthfact_2x2diagblock failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    T(3,1) = 0d0
    T(3,2) = T(2,1) 
    
    call d_orthfact_2x2diagblock(N,STR,Q,D,T(1:2,1:2),INFO)  

    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_orthfact_2x2diagblock failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if

    ! compute (A-E(1))(A-E(2))e_1
    COL(1) = T(1,1)*T(1,1) + T(1,2)*T(2,1) + rrho1*rrho2 - irho1*irho2 - T(1,1)*(rrho1+rrho2)
    COL(2) = T(2,1)*(T(1,1)+T(2,2)-(rrho1+rrho2))
    COL(3) = T(2,1)*T(3,2)
    
    ! compute first givens rotations
    call d_rot2_vec2gen(COL(2),COL(3),B1(1),B1(2),nrm1,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if

    ! compute second givens rotations
    call d_rot2_vec2gen(COL(1),nrm1,B2(1),B2(2),nrm2,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
  
  end if

end subroutine d_orthfact_buildbulge

