#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZARSUE (Zomplex Auxiliary Routine Sort Unimodular Eigenpairs)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine sorts the eigenvalues and optionally eigenvectors of a
! unitary matrix. The eigenpairs are sorted in ascending order of the
! argument of the eigenvalue. 
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
!  Z              COMPLEX(8) array of dimension (N,N)
!                   if JOB = 'E' unused
!                   if JOB = 'V' sort Z as well 
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies JOB is invalid
!                   INFO = -2 implies N is invalid
!                   INFO = -3 implies E is invalid
!                   INFO = -4 implies Z is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZARSUE(JOB,N,E,Z,INFO)

  ! input variables
  character, intent(in) :: JOB
  integer, intent(in) :: N
  real(8), intent(inout) :: E(2*N)
  complex(8), intent(inout) :: Z(N,N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii, jj, ind
  real(8), parameter :: PI = 3.14159265358979323846264338327950d0
  real(8) :: theta1, theta2
  real(8) :: tol, nrm, temp
  complex(8) :: ztmp
  
  ! set tol
  tol = 10d0*epsilon(1d0)
  
  ! initialize INFO
  INFO = 0

  ! check input in debug mode
  if (DEBUG) then
  
    ! check JOB
    if ((JOB.NE.'E').AND.(JOB.NE.'V')) then
      INFO = -1
      call UARERR(__FILE__,__LINE__,"JOB must be 'E' or 'V'",INFO,INFO)
      return
    end if
    
    ! check N
    call IARNAN(N,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"N is invalid",INFO,-2)
      return
    end if
    call IARINF(N,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"N is invalid",INFO,-2)
      return
    end if
    if (N < 2) then
      INFO = -2
      call UARERR(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
      return
    end if 
  
    ! check E
    call DARACH1(2*N,E,INFO)
    if (INFO.NE.0) then
      call UARERR(__FILE__,__LINE__,"E is invalid",INFO,-3)
      return
    end if
    do ii=1,N
      nrm = sqrt(E(2*ii-1)**2 + E(2*ii)**2)
      if (abs(nrm-1d0) > tol) then
        INFO = -3
        call UARERR(__FILE__,__LINE__,"E is not unimodular",INFO,INFO)
        return    
      end if
    end do
  
    ! check Z
    if (JOB.EQ.'V') then
      call ZARACH2(N,N,Z,INFO)
      if (INFO.NE.0) then
        call UARERR(__FILE__,__LINE__,"Z is invalid",INFO,-4)
        return
      end if
    end if 
  
  end if
  
  ! loop through E
  do ii=1,(N-1)
  
    ! find index of smallest arg
    ! initialize index
    ind = ii
    
    ! compute theta1
    call ZARNAR(E(2*ind-1),E(2*ind),theta1,INFO)
    
    ! check INFO in debug mode
    if (DEBUG) then
      call UARERR(__FILE__,__LINE__,"ZARNAR failed",INFO,INFO)
      if (INFO.NE.0) then 
        return 
      end if 
    end if
  
    ! check rest of E
    do jj=(ii+1),N
      ! compute theta2
      call ZARNAR(E(2*jj-1),E(2*jj),theta2,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call UARERR(__FILE__,__LINE__,"ZARNAR failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
      
      ! update if smaller
      if (theta2 < theta1) then
        ind = jj
        theta1 = theta2
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
        ztmp = Z(jj,ii)
        Z(jj,ii) = Z(jj,ind)
        Z(jj,ind) = ztmp
      end do
    end if
  end do

end subroutine ZARSUE
