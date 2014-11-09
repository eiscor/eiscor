!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZPFCHF (Zomplex unitary Plus rank 1 hessenberg Factored CHeck Factorization)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks the factorization input into ZUFFQR to make sure
! is represents a unitary hessenberg matrix to machine precision. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (3*N)
!                    array of generators for first sequence of Givens' 
!                    rotations
!
!  D               REAL(8) array of dimension (2*(N+2))
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
!  C               REAL(8) array of dimension (3*N)
!                    array of generators for second sequence of Givens' 
!                    rotations
!
!  B               REAL(8) array of dimension (3*N)
!                    array of generators for third sequence of Givens' 
!                    rotations
!
! OUTPUT VARIABLES:
!
!  INFO           INTEGER
!                   INFO equal to 0 implies the factorization is acceptable.
!                   INFO negative implies that an input variable has
!                   an improper value, i.e. INFO=-2 => Q is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZPFCHF(N,Q,D,C,B,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  real(8), intent(inout) :: Q(3*N), D(2*(N+1)), C(3*N), B(3*N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii
  real(8) :: tol, nrm
  
  ! set tol
  tol = 10d0*epsilon(1d0)
  
  ! initialize INFO
  INFO = 0
  
  ! check N
  call IARNAN(N,INFO)
  if (INFO.NE.0) then
    INFO = -1
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "N contains a NAN."
    write(*,*) ""
    return
  end if
  call IARINF(N,INFO)
  if (INFO.NE.0) then
    INFO = -1
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "N contains an INF."
    write(*,*) ""
    return
  end if
  if (N < 2) then
    INFO = -1
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "N must be at least 2."
    write(*,*) ""
    return
  end if
  
  ! check Q for NANs and INFs
  call DARACH1(3*N,Q,INFO)
  if (INFO.NE.0) then
    INFO = -2
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "Q contains an INF or NAN"
    write(*,*) ""
    return
  end if
  
  ! check Q for orthogonality
  do ii=1,N
    nrm = sqrt(Q(3*ii-2)**2 + Q(3*ii-1)**2 + Q(3*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -2
      write(*,*) "Error in "//__FILE__//" line:",__LINE__
      write(*,*) "Q is not unitary to working precision"
      write(*,*) ""
      return
   end if
  end do
  
  ! check D
  call DARACH1(2*(N+1),D,INFO)
  if (INFO.NE.0) then
    INFO = -3
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "D contains an INF or NAN"
    write(*,*) ""
    return
  end if
  
  ! check D for orthogonality
  do ii=1,(N+1)
    nrm = sqrt(D(2*ii-1)**2 + D(2*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -3
      write(*,*) "Error in "//__FILE__//" line:",__LINE__
      write(*,*) "D is not unitary to working precision"
      write(*,*) ""
      return
   end if
  end do
  
  ! check C for NANs and INFs
  call DARACH1(3*N,C,INFO)
  if (INFO.NE.0) then
    INFO = -4
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "C contains an INF or NAN"
    write(*,*) ""
    return
  end if
  
  ! check C for orthogonality
  do ii=1,N
    nrm = sqrt(C(3*ii-2)**2 + C(3*ii-1)**2 + C(3*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -4
      write(*,*) "Error in "//__FILE__//" line:",__LINE__
      write(*,*) "C is not unitary to working precision"
      write(*,*) ""
      return
   end if
  end do
  
  ! check C for zeros off-diagonals
  do ii=1,N
    if (abs(C(3*ii)).EQ.0d0) then
      INFO = -4
      write(*,*) "Error in "//__FILE__//" line:",__LINE__
      write(*,*) "C must not be diagonal."
      write(*,*) ""
      return
   end if
  end do
  
  ! check B for NANs and INFs
  call DARACH1(3*N,B,INFO)
  if (INFO.NE.0) then
    INFO = -5
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "B contains an INF or NAN"
    write(*,*) ""
    return
  end if
  
  ! check B for orthogonality
  do ii=1,N
    nrm = sqrt(B(3*ii-2)**2 + B(3*ii-1)**2 + B(3*ii)**2)
    if (abs(nrm-1d0) > tol) then
      INFO = -5
      write(*,*) "Error in "//__FILE__//" line:",__LINE__
      write(*,*) "B is not unitary to working precision"
      write(*,*) ""
      return
   end if
  end do
  
end subroutine ZPFCHF
