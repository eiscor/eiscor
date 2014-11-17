!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZPFCFT (Zomplex unitary Plus rank 1 hessenberg Factored Compute First Transformation)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a two by two diagonal block of an upper 
! hessenberg matrix that is the sum of a unitary matrix and a rank one 
! matrix. This sum is stored as a product of three sequences of Givens' 
! rotations and a complex unimodular diagonal matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  K               INTEGER
!                    index of the diagonal block
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
!  SHFT            COMPLEX(8) 
!                    contains the shift need for the first transformation
!
! OUTPUT VARIABLES:
!
!  B1              REAL(8) array of dimension (3)
!                    on exit contains the generators for the first
!                    transformation
!
!  INFO            INTEGER
!                    INFO equal to 0 implies successful computation.
!                    INFO negative implies that an input variable has
!                    an improper value, i.e. INFO=2 => N is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZPFCFT(N,K,Q,D,C,B,SHFT,B1,INFO)
  
  implicit none
  
  ! input variables
  integer, intent(in) :: N, K
  integer, intent(inout) :: INFO
  real(8), intent(in) :: Q(3*N), D(2*(N+1)), C(3*N), B(3*N)
  complex(8), intent(in) :: SHFT
  real(8), intent(inout) :: B1(3)
  
  ! compute variables
  real(8) :: nrm
  complex(8) :: block(2,2)
  
  ! initialize info
  INFO = 0
  
  ! get top block
  call ZPFTDB(N,K,Q,D,C,B,block,INFO) 
      
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "INFO =",INFO
    write(*,*) ""
    return
  end if
  
  ! shift first entry
  block(1,1) = block(1,1) - SHFT
  
  ! bulge
<<<<<<< HEAD
  call ZARCG43(dble(block(1,1)),aimag(block(1,1)),dble(block(2,1)), &
    aimag(block(2,1)),B1(1),B1(2),B1(3),nrm,INFO)
=======
  call ZARCG43(dble(block(1,1)),dimag(block(1,1)),dble(block(2,1)), &
    dimag(block(2,1)),B1(1),B1(2),B1(3),nrm,INFO)
>>>>>>> added files for complex unitary plus rank one from svn repo
      
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "INFO =",INFO
    write(*,*) ""
    return
  end if

end subroutine ZPFCFT
