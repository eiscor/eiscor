#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_factorcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks the factorization input into z_unifact_qr to make 
! sure it represents a unitary hessenberg matrix to machine precision. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                   INFO = 0 implies valid factorization
!                   INFO = -1 implies N is invalid
!                   INFO = -2 implies Q is invalid
!                   INFO = -3 implies D is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_factorcheck(N,Q,D,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(in) :: Q(3*N-3), D(2*N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii
  real(8),parameter :: tol = 10d0*epsilon(1d0) 
  real(8) :: nrm
  
  ! initialize INFO
  INFO = 0
  
  ! check N
  if (N < 2) then
    INFO = -1
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
    end if
    return
  end if
  
  ! check Q for NANs and INFs
  call d_1Darray_check(3*(N-1),Q,INFO)
  if (INFO.NE.0) then
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,INFO)
    end if
    INFO = -2
    return
  end if
  
  ! check Q for orthogonality
  do ii=1,(N-1)
    nrm = sqrt(Q(3*ii-2)**2 + Q(3*ii-1)**2 + Q(3*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -2
      ! print error in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"Q is not unitary",INFO,INFO)
      end if
      return
   end if
  end do
  
  ! check D
  call d_1Darray_check(2*N,D,INFO)
  if (INFO.NE.0) then
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,INFO)
    end if
    INFO = -3
    return
  end if
  
  ! check D for orthogonality
  do ii=1,N
    nrm = sqrt(D(2*ii-1)**2 + D(2*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -3
      ! print error in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"D is not unitary",INFO,INFO)
      end if
      return
   end if
  end do
  
end subroutine z_unifact_factorcheck
