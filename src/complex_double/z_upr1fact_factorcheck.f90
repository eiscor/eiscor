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
!  FLAG            LOGICAL
!                    .TRUE. second triangular factor is assumed nonzero
!                    .FALSE. second triangular factor is assumed to be identity
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D1,D2           REAL(8) arrays of dimension (2*(N+1))
!                    array of generators for complex diagonal matrices
!                    in the upper-triangular factors
!                    If FLAG = .FALSE., D2 is unused.
!
!  C1,B1,C2,B2     REAL(8) arrays of dimension (3*N)
!                    array of generators for upper-triangular parts of the pencil
!                    If FLAG = .FALSE., C2 and B2 are unused.
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                   INFO = 0 implies valid factorization
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
subroutine z_upr1fact_factorcheck(FLAG,N,Q,D1,C1,B1,D2,C2,B2,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: FLAG
  integer, intent(in) :: N
  real(8), intent(in) :: Q(3*(N-1)), D1(2*(N+1)), C1(3*N), B1(3*N)
  real(8), intent(in) :: D2(2*(N+1)), C2(3*N), B2(3*N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  logical :: flg
  integer :: ii
  real(8),parameter :: tol = 10d0*EISCOR_DBL_EPS 
  
  ! initialize INFO
  INFO = 0
  
  ! check N
  if (N < 2) then
    INFO = -2
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
    end if
    return
  end if
  
  ! check Q 
  call z_rot3array_check(N-1,Q,flg)
  if (.NOT.flg) then
    INFO = -3
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,INFO)
    end if
    return
  end if
  
  ! check D1 
  call d_rot2array_check(N+1,D1,flg)
  if (.NOT.flg) then
    INFO = -4
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"D1 is invalid",INFO,INFO)
    end if
    return
  end if

  ! check C1
  call z_rot3array_check(N,C1,flg)
  if (.NOT.flg) then
    INFO = -5
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"C1 is invalid",INFO,INFO)
    end if
    return
  end if

  ! check C1 for diagonal rotations
  do ii=1,N
    if (C1(3*ii).EQ.0) then
      INFO = -5
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"C1 contains a diagonal rotation",INFO,INFO)
      end if
      return
   end if
  end do

  ! check B1
  call z_rot3array_check(N,B1,flg)
  if (.NOT.flg) then
    INFO = -6
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"B1 is invalid",INFO,INFO)
    end if
    return
  end if

  ! check B1 for diagonal rotations
  do ii=1,N
    if (B1(3*ii).EQ.0) then
      INFO = -6
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"B1 contains a diagonal rotation",INFO,INFO)
      end if
      return
   end if
  end do

  ! check second triangular matrix
  if (FLAG) then

    ! check D2 
    call d_rot2array_check(N+1,D2,flg)
    if (.NOT.flg) then
      INFO = -7
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"D2 is invalid",INFO,INFO)
      end if
      return
    end if

    ! check C2
    call z_rot3array_check(N,C2,flg)
    if (.NOT.flg) then
      INFO = -8
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"C2 is invalid",INFO,INFO)
      end if
      return
    end if

    ! check C2 for diagonal rotations
    do ii=1,N
      if (C2(3*ii).EQ.0) then
        INFO = -8
        ! print error message in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"C2 contains a diagonal rotation",INFO,INFO)
        end if
        return
      end if
    end do

    ! check B2
    call z_rot3array_check(N,B2,flg)
    if (.NOT.flg) then
      INFO = -9
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"B2 is invalid",INFO,INFO)
      end if
      return
    end if

    ! check B2 for diagonal rotations
    do ii=1,N
      if (B2(3*ii).EQ.0) then
        INFO = -9
        ! print error message in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"B2 contains a diagonal rotation",INFO,INFO)
        end if
        return
      end if
    end do

  end if

end subroutine z_upr1fact_factorcheck
