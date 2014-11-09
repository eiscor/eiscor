!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZPFDTB (Zomplex unitary Plus rank 1 hessenberg Factored Deflate Two by two Block)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This an isolated 2x2 block of an upper hessenberg matrix that is the 
! sum of a unitary matrix and a rank one matrix. This sum is stored as 
! a product of three sequences of Givens' rotations and a complex 
! unimodular diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  COMPZ           CHARACTER
!                    'N': no eigenvectors
!                    'I': eigenvectors, initializing Z to the identity
!                    'V': eigenvectors, assume Z already initialized
!
!  N               INTEGER
!                    dimension of matrix
!
!  K               INTEGER
!                    index of the block to be deflated
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
!  Z               COMPLEX(8) array of dimension (N,N)
!                    if COMPZ = 'N' unused
!                    if COMPZ = 'I' stores eigenvectors in Z 
!                    if COMPZ = 'V' update Z to store eigenvectors 
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                   INFO equal to 0 implies successful computation.
!                   INFO negative implies that an input variable has
!                   an improper value, i.e. INFO=-2 => N is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZPFDTB(COMPZ,N,K,Q,D,C,B,Z,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: COMPZ
  integer, intent(in) :: N, K
  integer, intent(inout) :: INFO
  real(8), intent(inout) :: Q(3*N), D(2*(N+1)), C(3*N), B(3*N)
  complex(8), intent(inout) :: Z(N,N)
  
  ! compute variables
  integer :: ind1, ind2
  real(8) :: nrm, b1(3), binv(3)
  complex(8) :: block(2,2), eigs(2), temp(2,2)
  
  ! initialize INFO
  INFO = 0
  
  ! check COMPZ
  if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
    INFO = -1
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "COMPZ must be 'N', 'I' or 'V'"
    write(*,*) ""
    return
  end if
  
  ! check N
  if (N < 2) then
    INFO = -2
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "N must be at least 2"
    write(*,*) ""
    return
  end if
  
  ! check K
  if ((K < 1).OR.(K > N-1)) then
    INFO = -3
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "Must hold: 1 <= K <= N-1"
    write(*,*) ""
    return
  end if  
  
  ! get 2x2 block
  call ZPFTDB(N,K,Q,D,C,B,block,INFO) 
      
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "INFO =",INFO
    write(*,*) ""
    return
  end if
        
  ! compute eigenvalues and eigenvectors
  call ZTTEEV(block,eigs,temp,INFO)
      
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "INFO =",INFO
    write(*,*) ""
    return
  end if
      
  ! adjust eigenvectors to be a Givens' rotation
  call ZARCG43(dble(temp(1,1)),dimag(temp(1,1)),dble(temp(2,1)),dimag(temp(2,1)), &
        b1(1),b1(2),b1(3),nrm,INFO)    
      
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "INFO =",INFO
    write(*,*) ""
    return
  end if  
      
  ! update eigenvectors
  temp(1,1) = complex(b1(1),b1(2))
  temp(2,1) = complex(b1(3),0d0)
  temp(1,2) = -temp(2,1)
  temp(2,2) = conjg(temp(1,1))
  if ((COMPZ.EQ.'I').OR.(COMPZ.EQ.'V')) then
    Z(:,K:(K+1)) = matmul(Z(:,K:(K+1)),temp)
  end if
      
  ! set binv
  binv(1) = b1(1)
  binv(2) = -b1(2)
  binv(3) = -b1(3)
      
  ! fusion at top
  call ZUFFGR('T',N,K,K,Q,D,binv,INFO)
      
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "INFO =",INFO
    write(*,*) ""
    return
  end if
      
  ! pass through B
  ind1 = 3*(K-1) + 1
  ind2 = ind1+2
  call ZARGTO(B(ind1:ind2),B((ind1+3):(ind2+3)),b1,INFO)

  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "INFO =",INFO
    write(*,*) ""
    return
  end if
      
  ! pass through C
  call ZARGTO(C((ind1+3):(ind2+3)),C(ind1:ind2),b1,INFO)

  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "INFO =",INFO
    write(*,*) ""
    return
  end if
      
  ! fusion at bottom
  call ZUFFGR('B',N,K,K,Q,D,b1,INFO)
      
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "INFO =",INFO
    write(*,*) ""
    return
  end if
  
end subroutine ZPFDTB
