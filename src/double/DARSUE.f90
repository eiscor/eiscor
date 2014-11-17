#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DARSUE (Double Auxiliary Routine Sort Unimodular Eigenpairs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine sorts the eigenvalues and optionally eigenvectors of a
! real orthogonal matrix. Conjugate pairs C +/- Si are stored as 2x2
! diagonal blocks:
!
!  | C  -S |
!  | S   C |
!
! After sorting S will be positive and the blocks will be sorted in
! ascending order of their argument by the conjugate pair in the upper
! half plane C + Si.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  JOB            CHARACTER
!                    'E': just eigenvalues
!                    'V': eigenvectors also 
!
!  N              INTEGER
!                    dimension of matrix
!
!  E              REAL(8) array of dimension (2*N)
!                    array of eigenvalues
!
!  Z              REAL(8) array of dimension (N,N)
!                   if JOB = 'E' unused
!                   if JOB = 'V' sort Z as well 
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies JOB is invalid
!                    INFO = -2 implies N is invalid
!                    INFO = -3 implies E is invalid
!                    INFO = -4 implies Z is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DARSUE(JOB,N,E,Z,INFO)

  ! input variables
  character, intent(in) :: JOB
  integer, intent(in) :: N
  real(8), intent(inout) :: E(2*N), Z(N,N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii, jj, ind
  real(8), parameter :: PI = 3.14159265358979323846264338327950d0
  real(8) :: s1, s2
  real(8) :: tol, nrm, temp
  
  ! set tol
  tol = 10d0*epsilon(1d0)
  
  ! initialize INFO
  INFO = 0

  ! check COMPZ
  if ((JOB.NE.'E').AND.(JOB.NE.'V')) then
    INFO = -1
    ! print error in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"JOB must be 'E' or 'V'",INFO,INFO)
    end if   
    return
  end if
  
  ! check N
  if (N < 2) then
    INFO = -2
    ! print error in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
    end if   
    return
  end if  
  
  ! check E
  call DARACH1(2*N,E,INFO)
  if (INFO.NE.0) then
    ! print error in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"E is invalid",INFO,INFO)
    end if
    INFO = -3   
    return
  end if
  do ii=1,N
    nrm = sqrt(E(2*ii-1)**2 + E(2*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -3
      ! print error in debug mode
      if (DEBUG) then
        call UARERR(__FILE__,__LINE__,"E is not unimodular to machine precision",INFO,INFO)
      end if   
      return   
    end if
  end do
  
  ! check Z
  if (JOB.EQ.'V') then
    call DARACH2(N,N,Z,INFO)
    if (INFO.NE.0) then
      ! print error in debug mode
      if (DEBUG) then
        call UARERR(__FILE__,__LINE__,"Z is invalid",INFO,INFO)
      end if
      INFO = -4   
      return 
    end if
  end if 
  
  ! loop through E
  do ii=1,(N-1)
  
    ! find index of smallest arg
    ! initialize index
    ind = ii
  
    ! check rest of E
    do jj=(ii+1),N
      
      ! update if smaller
      if (E(2*jj-1) > E(2*ind-1)) then
        ind = jj
      else if ((E(2*jj-1).EQ.E(2*ind-1)).AND.(E(2*jj) > E(2*ind))) then 
        ind = jj
      end if
      
    end do
    
    ! swap eigenvalues
    if (ind.NE.ii) then
      temp = E(2*ii-1)
      E(2*ii-1) = E(2*ind-1)
      E(2*ind-1) = temp 
      temp = E(2*ii)
      E(2*ii) = E(2*ind)
      E(2*ind) = temp
    end if
    
    ! swap eigenvectors
    if ((ind.NE.ii).AND.(JOB.EQ.'V')) then
      do jj=1,N
        temp = Z(jj,ii)
        Z(jj,ii) = Z(jj,ind)
        Z(jj,ind) = temp
      end do
    end if
  end do

end subroutine DARSUE
