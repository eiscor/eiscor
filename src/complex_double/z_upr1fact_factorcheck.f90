#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_factorcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks the factorization input into z_upr1fact_twistedqz
! to make sure it represents an extended hessenberg triangular pencil
! to machine precision. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  ALG             CHARACTER(2)
!                    'QR': second triangular factor is assumed to be identity
!                    'QZ': second triangular factor is assumed nonzero
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D               REAL(8) array of dimension (2,2*(N+1))
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  R               REAL(8) array of dimension (4,3*N)
!                    array of generators for upper-triangular parts
!                    of the pencil
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                   INFO = 0 implies valid factorization
!                   INFO = -1 implies ALG is invalid
!                   INFO = -2 implies N is invalid
!                   INFO = -3 implies Q is invalid
!                   INFO = -4 implies D is invalid
!                   INFO = -5 implies R is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_factorcheck(ALG,N,Q,D,R,INFO)

  implicit none
  
  ! input variables
  character(2), intent(in) :: ALG
  integer, intent(in) :: N
  real(8), intent(in) :: Q(3*(N-1)), D(2,2*(N+1)), R(4,3*N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii
  real(8),parameter :: tol = 10d0*epsilon(1d0) 
  real(8) :: nrm
  
  ! initialize INFO
  INFO = 0
  
  ! check ALG
  if ((ALG.NE.'QR').AND.(ALG.NE.'QZ')) then
    INFO = -1
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"ALG must be 'QR' or 'QZ'",INFO,INFO)
    end if
    return
  end if
  
  ! check N
  if (N < 2) then
    INFO = -2
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
    end if
    return
  end if
  
  ! check Q for NANs and INFs
  call d_1Darray_check(3*(N-1),Q,INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,INFO)
    end if
    INFO = -3
    return
  end if
  
  ! check Q for orthogonality
  do ii=1,(N-1)
    nrm = sqrt(Q(3*ii-2)**2 + Q(3*ii-1)**2 + Q(3*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -3
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"Q is not unitary",INFO,INFO)
      end if
      return
   end if
  end do
  
  ! check D
  call d_1Darray_check(2*(N+1),D(1,:),INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,INFO)
    end if
    INFO = -4
    return
  end if
  if (ALG.EQ.'QZ') then
    call d_1Darray_check(2*(N+1),D(2,:),INFO)
    if (INFO.NE.0) then
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,INFO)
      end if
      INFO = -4
      return
    end if
  end if
  
  ! check D for orthogonality
  do ii=1,(N+1)
    nrm = sqrt(D(1,2*ii-1)**2 + D(1,2*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -4
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"D is not unitary",INFO,INFO)
      end if
      return
   end if
  end do
  if (ALG.EQ.'QZ') then
    do ii=1,(N+1)
      nrm = sqrt(D(2,2*ii-1)**2 + D(2,2*ii)**2)
      if (abs(nrm-1d0) > tol) then
        INFO = -4
        ! print error message in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"D is not unitary",INFO,INFO)
        end if
        return
     end if
    end do
  end if
  
  ! check R for NANs and INFs
  call d_2Darray_check(2,3*N,R(1:2,:),INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"R is invalid",INFO,INFO)
    end if
    INFO = -5
    return
  end if
  if (ALG.EQ.'QZ') then
    call d_2Darray_check(2,3*N,R(3:4,:),INFO)
    if (INFO.NE.0) then
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"R is invalid",INFO,INFO)
      end if
      INFO = -5
      return
    end if
  end if
  
  ! check R for orthogonality
  do ii=1,N
    nrm = sqrt(R(1,3*ii-2)**2 + R(1,3*ii-1)**2 + R(1,3*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -5
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"R is not unitary",INFO,INFO)
      end if
      return
   end if
  end do
  do ii=1,N
    nrm = sqrt(R(2,3*ii-2)**2 + R(2,3*ii-1)**2 + R(2,3*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -5
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"R is not unitary",INFO,INFO)
      end if
      return
   end if
  end do
  if (ALG.EQ.'QZ') then
    do ii=1,N
      nrm = sqrt(R(3,3*ii-2)**2 + R(3,3*ii-1)**2 + R(3,3*ii)**2)
      if (abs(nrm-1d0) > tol) then
        INFO = -5
        ! print error message in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"R is not unitary",INFO,INFO)
        end if
        return
     end if
    end do
    do ii=1,N
      nrm = sqrt(R(4,3*ii-2)**2 + R(4,3*ii-1)**2 + R(4,3*ii)**2)
      if (abs(nrm-1d0) > tol) then
        INFO = -5
        ! print error message in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"R is not unitary",INFO,INFO)
        end if
        return
     end if
    end do  
  end if
  
  ! check R is upper-triangular
  do ii=1,N
    if (abs(R(2,3*ii)) <= tol) then
      INFO = -5
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"R is not upper-triangular",INFO,INFO)
      end if
      return
   end if
  end do
  if (ALG.EQ.'QZ') then  
    do ii=1,N
      if (abs(R(4,3*ii)) <= tol) then
        INFO = -5
        ! print error message in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"R is not upper-triangular",INFO,INFO)
        end if
        return
     end if
    end do
  end if
  
end subroutine z_upr1fact_factorcheck
