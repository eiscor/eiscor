!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZPFFQR (Zomplex unitary Plus rank 1 hessenberg Factored Fast QR)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the schur decomposition of an upper hessenberg 
! matrix that is the sum of a unitary matrix and a rank one matrix.
! This sum is stored as a product of three sequences of Givens' 
! rotations and a complex unimodular diagonal matrix. 
! 
! This factorization is described in:
!
!
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
!  Q               REAL(8) array of dimension (3*N)
!                    array of generators for first sequence of Givens' 
!                    rotations
!
!  D               REAL(8) array of dimension (2*(N+2))
!                    array of generators for complex diagonal matrix
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
!  Z              COMPLEX(8) array of dimension (N,N)
!                   if COMPZ = 'N' unused
!                   if COMPZ = 'I' stores eigenvectors in Z 
!                   if COMPZ = 'V' update Z to store eigenvectors 
!
!  ITS            INTEGER array of dimension (N-1)
!                   Contains the number of iterations per deflation
!
!  INFO           INTEGER
!                   INFO = 0 implies successful computation.
!                   INFO = 1 implies no convergence in maximum 
!                   allowed iterations.
!                   INFO negative implies that an input variable has
!                   an improper value, i.e. INFO=-2 => N is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZPFFQR(COMPZ,N,Q,D,C,B,Z,ITS,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: COMPZ
  integer, intent(in) :: N
  real(8), intent(inout) :: Q(3*N), D(2*(N+1)), C(3*N), B(3*N)
  integer, intent(inout) :: INFO, ITS(N-1)
  complex(8), intent(inout) :: Z(N,N)
  
  ! compute variables
  integer :: ii, jj, kk, ind1, ind2, ll, strt, k
  integer :: start_index, stop_index, zero_index, bc_count, it_max, it_count
  real(8) :: nrm, b1(3), binv(3)
  complex(8) :: block(2,2), eigs(2), temp(2,2)
  
  ! initialize info
  INFO = 0
  
  ! check COMPZ
  if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
    INFO = -1
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "COMPZ must be 'N', 'I' or 'V'"
    write(*,*) ""
    return
  end if
  
  ! check factorization
  call ZPFCHF(N,Q,D,C,B,INFO)
  if (INFO.NE.0) then
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "Invalid factorization."
    write(*,*) ""
    return
  end if
  
  ! check Z
  if (COMPZ.EQ.'V') then
    call ZARACH2(N,N,Z,INFO)
    if (INFO.NE.0) then
      INFO = -5
      write(*,*) "Error in "//__FILE__//" line:",__LINE__
      write(*,*) "Z contains an INF or NAN"
      write(*,*) ""
      return
    end if
  end if   

  ! initialize storage
  ITS = 0
  
  if (COMPZ.EQ.'I') then
<<<<<<< HEAD
    Z = cmplx(0d0,0d0,kind=8)
    do ii=1,n
      Z(ii,ii) = cmplx(1d0,0d0,kind=8)
=======
    Z = complex(0d0,0d0)
    do ii=1,n
      Z(ii,ii) = complex(1d0,0d0)
>>>>>>> added files for complex unitary plus rank one from svn repo
    end do
  end if
  
  ! initialize indices
  start_index = 1
  stop_index = N-1
  zero_index = 0
  it_max = 20*N
  it_count = 0
  bc_count = N-1
  
  ! loop for bulgechasing
  do kk=1,it_max
  
    ! check for completion
    if(stop_index <= 0)then
      exit
    end if

    ! check for deflation
    !call ZPFGRD(N,start_index,stop_index,zero_index,Q,D,C,B,it_count,ITS,INFO)
    call ZUFGRD(N,start_index,stop_index,zero_index,Q,D,it_count,ITS,INFO)
    
    ! check INFO
    if (INFO .NE. 0) then
      write(*,*) "Error in "//__FILE__//" line:",__LINE__
      write(*,*) "INFO =",INFO
      write(*,*) ""
      return
    end if
       
    ! if 1x1 block remove and check again 
    if(stop_index == zero_index)then
                
      ! update indices
      stop_index = stop_index - 1
      zero_index = 0
      start_index = 1

    ! if 2x2 block remove and check again
    else if(stop_index-1 == zero_index)then
     
      ! deflate 2x2 block
      call ZPFDTB(COMPZ,N,stop_index,Q,D,C,B,Z,INFO)
      
      ! check INFO
      if (INFO .NE. 0) then
        write(*,*) "Error in "//__FILE__//" line:",__LINE__
        write(*,*) "INFO =",INFO
        write(*,*) ""
        return
      end if
           
    ! if greater than 2x2 chase a bulge and check again
    else
        
      ! chase bulge
      call ZPFSFI(COMPZ,N,start_index,stop_index,Q,D,C,B,bc_count,Z,it_count,INFO)

      ! check INFO
      if (INFO .NE. 0) then
        write(*,*) "Error in "//__FILE__//" line:",__LINE__
        write(*,*) "INFO =",INFO
        write(*,*) ""
        return
      end if
                
    end if
    
    ! if it_max hit
    if (kk == it_max) then
      INFO = 1
      its(stop_index) = it_count
    end if
    
  end do

end subroutine ZPFFQR
