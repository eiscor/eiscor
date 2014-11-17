!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZPFSFI (Zomplex unitary Plus rank 1 hessenberg Factored Singleshift Francis Iteration)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift 
! algorithm on an upper hessenberg matrix that is the sum of a unitary 
! matrix and a rank one matrix. This sum is stored as a product of 
! three sequences of Givens' rotations and a complex unimodular 
! diagonal matrix. 
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
!  STR             INTEGER
!                    index of the top most givens rotation where 
!                    the iteration begins
!
!  STP             INTEGER
!                    index of the bottom most givens rotation where 
!                    the iteration ends
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
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                   INFO equal to 0 implies successful computation.
!                   INFO negative implies that an input variable has
!                   an improper value, i.e. INFO=-2 => N is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZPFSFI(COMPZ,N,STR,STP,Q,D,C,B,BCCNT,Z,ITCNT,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: COMPZ
  integer, intent(in) :: N, STR, STP
  integer, intent(inout) :: BCCNT, ITCNT, INFO
  real(8), intent(inout) :: Q(3*N), D(2*(N+1)), C(3*N), B(3*N)
  complex(8), intent(inout) :: Z(N,N)
  
  ! compute variables
  integer :: ii, ind1, ind2
  real(8) :: s1, s2
  real(8) :: bulge(3),binv(3)
  complex(8) :: shift
  complex(8) :: block(2,2), temp(2,2), eigs(2)
  
  ! initialize INFO
  INFO = 0
  
  ! check COMPZ
!  if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
!    INFO = -1
!    write(*,*) "Error in "//__FILE__//" line:",__LINE__
!    write(*,*) "COMPZ must be 'N', 'I' or 'V'"
!    write(*,*) ""
!    return
!  end if
  
  ! check STR
!  if ((STR < 1).OR.(STR > N-1)) then
!    INFO = -3
!    write(*,*) "Error in "//__FILE__//" line:",__LINE__
!    write(*,*) "Must hold: 1 <= STR <= N-1"
!    write(*,*) ""
!    return
!  end if
  
  ! check STP
!  if ((STP < STR).OR.(STP > N-1)) then
!    INFO = -4
!    write(*,*) "Error in "//__FILE__//" line:",__LINE__
!    write(*,*) "Must hold: STR <= STP <= N-1"
!    write(*,*) ""
!    return
!  end if
  
  ! compute a nonzero shift
  ! random shift
  if(mod(ITCNT+1,11) == 0)then
    call random_number(s1)
    call random_number(s2)
<<<<<<< HEAD
    shift = cmplx(s1,s2,kind=8)
=======
    shift = complex(s1,s2)
>>>>>>> added files for complex unitary plus rank one from svn repo
          
  ! wilkinson shift
  else
    ! get 2x2 block
    call ZPFTDB(N,STP,Q,D,C,B,block,INFO) 
      
    ! check INFO
!    if (INFO .NE. 0) then
!      write(*,*) "Error in "//__FILE__//" line:",__LINE__
!      write(*,*) "INFO =",INFO
!      write(*,*) ""
!      return
!    end if
        
    ! compute eigenvalues and eigenvectors
    call ZTTEEV(block,eigs,temp,INFO)
      
    ! check INFO
!    if (INFO .NE. 0) then
!      write(*,*) "Error in "//__FILE__//" line:",__LINE__
!      write(*,*) "INFO =",INFO
!      write(*,*) ""
!      return
!    end if
          
    ! choose wikinson shift
    if(abs(block(2,2)-eigs(1)) < abs(block(2,2)-eigs(2)))then
      shift = eigs(1)
    else
      shift = eigs(2)
    end if

  end if

  ! build bulge
  call ZPFCFT(N,STR,Q,D,C,B,shift,bulge,INFO)
        
  ! check INFO
!  if (INFO .NE. 0) then
!    write(*,*) "Error in "//__FILE__//" line:",__LINE__
!    write(*,*) "INFO =",INFO
!    write(*,*) ""
!    return
!  end if
<<<<<<< HEAD
  
=======
	
>>>>>>> added files for complex unitary plus rank one from svn repo
  ! bulge inverse
  binv(1) = bulge(1)
  binv(2) = -bulge(2)
  binv(3) = -bulge(3)
  
  ! fusion at top
  call ZUFFGR('T',N,STR,STP,Q,D,binv,INFO)
  
  ! check INFO
!  if (INFO .NE. 0) then
!    write(*,*) "Error in "//__FILE__//" line:",__LINE__
!    write(*,*) "INFO =",INFO
!    write(*,*) ""
!    return
!  end if
  
  ! main chasing loop
  do ii=STR,(STP-1)
  
<<<<<<< HEAD
    ! update eigenvectors
    if (COMPZ .NE. 'N')then
      temp(1,1) = cmplx(bulge(1),bulge(2),kind=8)
      temp(2,1) = cmplx(bulge(3),0d0,kind=8)
      temp(1,2) = -temp(2,1)
      temp(2,2) = conjg(temp(1,1))
      Z(:,ii:(ii+1)) = matmul(Z(:,ii:(ii+1)),temp)
    end if
    
    ! check BCCNT
    if (ii >= BCCNT) then
    
=======
  	! update eigenvectors
  	if (COMPZ .NE. 'N')then
  	  temp(1,1) = complex(bulge(1),bulge(2))
  	  temp(2,1) = complex(bulge(3),0d0)
  	  temp(1,2) = -temp(2,1)
  	  temp(2,2) = conjg(temp(1,1))
  	  Z(:,ii:(ii+1)) = matmul(Z(:,ii:(ii+1)),temp)
  	end if
  	
  	! check BCCNT
  	if (ii >= BCCNT) then
  	
>>>>>>> added files for complex unitary plus rank one from svn repo
      ! pass through B
      ind1 = 3*(ii-1) + 1
      ind2 = ind1+2
      call ZARGTO(B(ind1:ind2),B((ind1+3):(ind2+3)),bulge,INFO)

      ! check INFO
  !    if (INFO .NE. 0) then
  !      write(*,*) "Error in "//__FILE__//" line:",__LINE__
  !      write(*,*) "INFO =",INFO
  !      write(*,*) ""
  !      return
  !    end if
        
      ! pass through C
      call ZARGTO(C((ind1+3):(ind2+3)),C(ind1:ind2),bulge,INFO)

      ! check INFO
  !    if (INFO .NE. 0) then
  !      write(*,*) "Error in "//__FILE__//" line:",__LINE__
  !      write(*,*) "INFO =",INFO
  !      write(*,*) ""
  !      return
  !    end if
  
      ! update BCCNT
      BCCNT = BCCNT - 1
  
    endif
     
    ! set indices
    ind1 = 2*(ii-1) + 1
    ind2 = ind1+3
     
    ! through diag
    call ZARGTD('R',D(ind1:ind2),bulge,INFO)

    ! check INFO
!    if (INFO .NE. 0) then
!      write(*,*) "Error in "//__FILE__//" line:",__LINE__
!      write(*,*) "INFO =",INFO
!      write(*,*) ""
!      return
!    end if
    
    ! set indices
    ind1 = 3*(ii-1) + 1
    ind2 = ind1+2
     
    ! through Q
    call ZARGTO(Q(ind1:ind2),Q((ind1+3):(ind2+3)),bulge,INFO)

    ! check INFO
!    if (INFO .NE. 0) then
!      write(*,*) "Error in "//__FILE__//" line:",__LINE__
!      write(*,*) "INFO =",INFO
!      write(*,*) ""
!      return
!    end if

  end do
  
  ! update eigenvectors
  if (COMPZ .NE. 'N')then
<<<<<<< HEAD
    temp(1,1) = cmplx(bulge(1),bulge(2),kind=8)
    temp(2,1) = cmplx(bulge(3),0d0,kind=8)
    temp(1,2) = -temp(2,1)
    temp(2,2) = conjg(temp(1,1))
    Z(:,STP:(STP+1)) = matmul(Z(:,STP:(STP+1)),temp)
=======
    temp(1,1) = complex(bulge(1),bulge(2))
  	temp(2,1) = complex(bulge(3),0d0)
  	temp(1,2) = -temp(2,1)
  	temp(2,2) = conjg(temp(1,1))
  	Z(:,STP:(STP+1)) = matmul(Z(:,STP:(STP+1)),temp)
>>>>>>> added files for complex unitary plus rank one from svn repo
  end if
  
  ! pass through B
  ind1 = 3*(STP-1) + 1
  ind2 = ind1+2
  call ZARGTO(B(ind1:ind2),B((ind1+3):(ind2+3)),bulge,INFO)

  ! check INFO
!  if (INFO .NE. 0) then
!    write(*,*) "Error in "//__FILE__//" line:",__LINE__
!    write(*,*) "INFO =",INFO
!    write(*,*) ""
!    return
!  end if
      
  ! pass through C
  call ZARGTO(C((ind1+3):(ind2+3)),C(ind1:ind2),bulge,INFO)

  ! check INFO
!  if (INFO .NE. 0) then
!    write(*,*) "Error in "//__FILE__//" line:",__LINE__
!    write(*,*) "INFO =",INFO
!    write(*,*) ""
!    return
!  end if
  
  ! fusion at bottom
  call ZUFFGR('B',N,STR,STP,Q,D,bulge,INFO)
  
  ! check INFO
!  if (INFO .NE. 0) then
!    write(*,*) "Error in "//__FILE__//" line:",__LINE__
!    write(*,*) "INFO =",INFO
!    write(*,*) ""
!    return
!  end if
  
  ! update ITCNT
  ITCNT = ITCNT + 1

end subroutine ZPFSFI
