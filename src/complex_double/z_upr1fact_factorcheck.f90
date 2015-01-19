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
!  D1, D2          REAL(8) array of dimension (2*(N+1))
!                    arrays of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1, B1, C2, B2  REAL(8) array of dimension (4,3*N)
!                    arrays of generators for upper-triangular parts
!                    of the pencil
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                   INFO = 1 implies subroutine failed
!                   INFO = 0 implies valid factorization
!                   INFO = -1 implies ALG is invalid
!                   INFO = -2 implies N is invalid
!                   INFO = -3 implies Q is invalid
!                   INFO = -4 implies D1 is invalid
!                   INFO = -5 implies C1 is invalid
!                   INFO = -6 implies B1 is invalid
!                   INFO = -7 implies D2 is invalid
!                   INFO = -8 implies C2 is invalid
!                   INFO = -9 implies B2 is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_factorcheck(ALG,N,Q,D1,C1,B1,D2,C2,B2,INFO)

  implicit none
  
  ! input variables
  character(2), intent(in) :: ALG
  integer, intent(in) :: N
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*(N+1)), D2(2*(N+1))
  real(8), intent(inout) :: C1(3*N), B1(3*N), C2(3*N) ,B2(3*N)
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
  
  ! check D1 for INFs or NANs
  call d_1Darray_check(2*(N+1),D1,INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"D1 is invalid",INFO,INFO)
    end if
    INFO = -4
    return
  end if

  ! check D1 for orthogonality
  do ii=1,(N+1)
    nrm = sqrt(D1(2*ii-1)**2 + D1(2*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -4
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"D1 is not unitary",INFO,INFO)
      end if
      return
    end if
  end do  
  
  ! check C1 for NANs and INFs
  call d_1Darray_check(3*N,C1,INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"C1 is invalid",INFO,INFO)
    end if
    INFO = -5
    return
  end if
  
  ! check C1 for orthogonality
  do ii=1,N
    nrm = sqrt(C1(3*ii-2)**2 + C1(3*ii-1)**2 + C1(3*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -5
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"C1 is not unitary",INFO,INFO)
      end if
      return
    end if
  end do
  
  ! check C1 for inf diagonal elements
  do ii=1,N
    if (abs(C1(3*ii)) <= tol) then
      INFO = -5
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"C1 has an infinite diagonal element",INFO,INFO)
      end if
      return
    end if
  end do
  
  ! check B1 for NANs and INFs
  call d_1Darray_check(3*N,B1,INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"B1 is invalid",INFO,INFO)
    end if
    INFO = -6
    return
  end if
  
  ! check B1 for orthogonality
  do ii=1,N
    nrm = sqrt(B1(3*ii-2)**2 + B1(3*ii-1)**2 + B1(3*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -6
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"B1 is not unitary",INFO,INFO)
      end if
      return
    end if
  end do
  
  ! check B1 for zero diagonal elements
  do ii=1,N
    if (abs(B1(3*ii)) <= tol) then
      INFO = -6
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"B1 has a zero diagonal element",INFO,INFO)
      end if
      return
    end if
  end do
  
  ! check D2, C2, B2 if ALG == QZ
  if (ALG.EQ.'QZ') then
  
    ! check D2 for INFs or NANs
    call d_1Darray_check(2*(N+1),D2,INFO)
    if (INFO.NE.0) then
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"D2 is invalid",INFO,INFO)
      end if
      INFO = -7
      return
    end if
  
    ! check D2 for orthogonality
    do ii=1,(N+1)
      nrm = sqrt(D2(2*ii-1)**2 + D2(2*ii)**2)
      if (abs(nrm-1d0) > tol) then
        INFO = -7
        ! print error message in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"D2 is not unitary",INFO,INFO)
        end if
        return
      end if
    end do
  
   ! check C2 for NANs and INFs  
    call d_1Darray_check(3*N,C2,INFO)
    if (INFO.NE.0) then
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"C2 is invalid",INFO,INFO)
      end if
      INFO = -8
      return
    end if
    
    ! check C2 for orthogonality
    do ii=1,N
      nrm = sqrt(C2(3*ii-2)**2 + C2(3*ii-1)**2 + C2(3*ii)**2)
      if (abs(nrm-1d0) > tol) then
        INFO = -8
        ! print error message in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"C2 is not unitary",INFO,INFO)
        end if
        return
      end if
    end do
     
    ! check C2 for inf diagonal elements
    do ii=1,N
      if (abs(C2(3*ii)) <= tol) then
        INFO = -8
        ! print error message in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"C2 has an infinite diagonal element",INFO,INFO)
        end if
        return
     end if
    end do
  
   ! check B2 for NANs and INFs  
    call d_1Darray_check(3*N,B2,INFO)
    if (INFO.NE.0) then
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"B2 is invalid",INFO,INFO)
      end if
      INFO = -9
      return
    end if
 
    ! check B2 for orthogonality
    do ii=1,N
      nrm = sqrt(B2(3*ii-2)**2 + B2(3*ii-1)**2 + B2(3*ii)**2)
      if (abs(nrm-1d0) > tol) then
        INFO = -9
        ! print error message in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"B2 is not unitary",INFO,INFO)
        end if
        return
     end if
    end do  

    ! check B2 for zero diagonal elements
    do ii=1,N
      if (abs(B2(3*ii)) <= tol) then
        INFO = -9
        ! print error message in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"B2 has a zero diagonal element",INFO,INFO)
        end if
        return
     end if
    end do
  
  end if 

end subroutine z_upr1fact_factorcheck
