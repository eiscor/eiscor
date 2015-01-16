#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_rot3throughtri
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine passes a rotation through an upper-triangular matrix 
! stored as the product of a diagonal matrix and 2 sequences of 
! rotations.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  DIR             CHARACTER(3)
!                    'L2R': pass rotation from left to right
!                    'R2L': pass rotation from right to left
!
!  N               INTEGER
!                    dimension of matrix
!
!  K               INTEGER
!                    index where rotation is passed through
!
!  D               REAL(8) array of dimension (2*(N+1))
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C               REAL(8) array of dimension (3*N)
!                    array of generators for upper-triangular parts
!
!  B               REAL(8) array of dimension (3*N)
!                    array of generators for upper-triangular parts
!
!  G               REAL(8) array of dimension (3)
!                    generator for rotation
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies DIR is invalid
!                   INFO = -2 implies N is invalid
!                   INFO = -3 implies K is invalid
!                   INFO = -4 implies D is invalid
!                   INFO = -5 implies C is invalid
!                   INFO = -6 implies B is invalid
!                   INFO = -7 implies G is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_rot3throughtri(DIR,N,K,D,C,B,G,INFO)

  implicit none
  
  ! input variables
  character(3), intent(in) :: DIR
  integer, intent(in) :: N, K
  real(8), intent(inout) :: D(2*(N+1)), C(3*N), B(3*N), G(3)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ind1, ind2, ii
  
  ! initialize info
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
  
    ! check DIR
    if ((DIR.NE.'L2R').AND.(DIR.NE.'R2L')) then
      INFO = -1
      call u_infocode_check(__FILE__,__LINE__,"DIR must be 'L2R' or 'R2L'",INFO,INFO)
      return
    end if
    
    ! check N
    if (N < 2) then
      INFO = -2
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
      return
    end if 
  
    ! check K
    if ((K < 1).OR.(K > N-1)) then
      INFO = -3
      call u_infocode_check(__FILE__,__LINE__,"K must 1 <= K <= N-1",INFO,INFO)
      return
    end if 
    
    ! check D
    call z_1Darray_check(2*(N+1),D,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,-4)
      return
    end if
  
    ! check C
    call z_1Darray_check(3*N,C,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"C is invalid",INFO,-5)
      return
    end if
    
    ! check B
    call z_1Darray_check(3*N,B,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"B is invalid",INFO,-6)
      return
    end if
    
    ! check G
    call z_1Darray_check(3,G,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"G is invalid",INFO,-7)
      return
    end if

  end if
  
  ! L2R
  if (DIR.EQ.'L2R') then
  
    ! through D
    ind1 = 2*K-1
    ind2 = ind1+3
    call z_rot3_swapdiag('L',D(ind1:ind2),G,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_swapdiag failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! though C
    ind1 = 3*K-2
    ind2 = ind1+2
    call z_rot3_turnover(C(ind1:ind2),C((ind1+3):(ind2+3)),G,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_turnover",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! through B
    call z_rot3_turnover(B((ind1+3):(ind2+3)),B(ind1:ind2),G,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_turnover",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
  
  ! R2L
  else
  
    ! through B
    ind1 = 3*K-2
    ind2 = ind1+2
    call z_rot3_turnover(B(ind1:ind2),B((ind1+3):(ind2+3)),G,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_turnover",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! through C
    call z_rot3_turnover(C((ind1+3):(ind2+3)),C(ind1:ind2),G,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_turnover",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
    
    ! through D
    ind1 = 2*K-1
    ind2 = ind1+3
    call z_rot3_swapdiag('R',D(ind1:ind2),G,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_swapdiag failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
  
  end if

end subroutine z_upr1fact_rot3throughtri
