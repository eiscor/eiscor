!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZPFTDB (Zomplex unitary Plus rank 1 hessenberg Factored Two by two Diagonal Block)
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
! OUTPUT VARIABLES:
!
<<<<<<< HEAD
!  H              complex(8) array of dimension (2,2)
=======
!  H              COMPLEX(8) array of dimension (2,2)
>>>>>>> added files for complex unitary plus rank one from svn repo
!                   on exit contains the desired 2x2 block
!
!  INFO           INTEGER
!                   INFO equal to 0 implies successful computation.
!                   INFO negative implies that an input variable has
!                   an improper value, i.e. INFO=-2 => K is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine ZPFTDB(N,K,Q,D,C,B,H,INFO)
  
  implicit none
  
  ! input variables
  integer, intent(in) :: N, K
  integer, intent(inout) :: INFO
  real(8), intent(in) :: Q(3*N), D(2*(N+1)), C(3*N), B(3*N)
  complex(8), intent(inout) :: H(2,2)
  
  ! compute variables
  integer :: ind
  complex(8) :: R(3,2)
  
  ! initialize INFO
  INFO = 0
  
  ! check N
  if (N < 2) then
    INFO = -1
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "N must be at least 2"
    write(*,*) ""
    return
  end if
  
  ! check K
  if ((K < 1).OR.(K > N-1)) then
    INFO = -2
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "Must hold: 1 <= K <= N-1"
    write(*,*) ""
    return
  end if
  
  ! check C
  if (abs(C(3*(K-1)+3)).EQ.0d0) then
    INFO = -5
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "C is not properly upper hessenberg."
    write(*,*) ""
    return 
  end if  
  if (abs(C(3*(K)+3)).EQ.0d0) then
    INFO = -5
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "C is not properly upper hessenberg."
    write(*,*) ""
    return 
  end if 
  if (K > 1) then
    if (abs(C(3*(K-2)+3)).EQ.0d0) then
      INFO = -5
      write(*,*) "Error in "//__FILE__//" line:",__LINE__
      write(*,*) "C is not properly upper hessenberg."
      write(*,*) ""
      return 
    end if  
  end if
  
  ! initialize 
<<<<<<< HEAD
  R = cmplx(0d0,0d0,kind=8)

  ! first column of R
  ind = 3*(K-1)
  R(2,1) = cmplx(-B(ind+3)/C(ind+3),0d0,kind=8)

  ! if not at top  
  if (K > 1) then
    R(1,1) = (cmplx(-B(ind-2),B(ind-1),kind=8)*cmplx(B(ind+1),B(ind+2),kind=8) &
      + R(2,1)*cmplx(C(ind-2),C(ind-1),kind=8)*cmplx(C(ind+1),-C(ind+2),kind=8))/cmplx(C(ind),0d0,kind=8)
  end if
      
  ! second column of R
  ind = 3*K
  R(3,2) = cmplx(-B(ind+3)/C(ind+3),0d0)
  R(2,2) = (cmplx(-B(ind-2),B(ind-1),kind=8)*cmplx(B(ind+1),B(ind+2),kind=8) &
      + R(3,2)*cmplx(C(ind-2),C(ind-1),kind=8)*cmplx(C(ind+1),-C(ind+2),kind=8))/cmplx(C(ind),0d0,kind=8)
  
  ! if not at top
  if (K > 1) then    
    R(1,2) = (cmplx(B(ind-5),-B(ind-4),kind=8)*cmplx(B(ind),0d0,kind=8)*cmplx(B(ind+1),B(ind+2),kind=8) - &
      cmplx(C(ind-5),C(ind-4),kind=8)/cmplx(C(ind),0d0,kind=8)* &
      (cmplx(C(ind-2),-C(ind-1),kind=8)*cmplx(B(ind-2),-B(ind-1),kind=8)*cmplx(B(ind+1),B(ind+2),kind=8) - &
      cmplx(C(ind+1),-C(ind+2),kind=8)*R(3,2)))/cmplx(C(ind-3),0d0,kind=8)
=======
  R = complex(0d0,0d0)
	
  ! first column of R
  ind = 3*(K-1)
  R(2,1) = complex(-B(ind+3)/C(ind+3),0d0)
	
	! if not at top	
  if (K > 1) then
	  R(1,1) = (complex(-B(ind-2),B(ind-1))*complex(B(ind+1),B(ind+2)) &
			+ R(2,1)*complex(C(ind-2),C(ind-1))*complex(C(ind+1),-C(ind+2)))/complex(C(ind),0d0)
  end if
			
  ! second column of R
  ind = 3*K
  R(3,2) = complex(-B(ind+3)/C(ind+3),0d0)
  R(2,2) = (complex(-B(ind-2),B(ind-1))*complex(B(ind+1),B(ind+2)) &
			+ R(3,2)*complex(C(ind-2),C(ind-1))*complex(C(ind+1),-C(ind+2)))/complex(C(ind),0d0)
	
	! if not at top
	if (K > 1) then		
		R(1,2) = (complex(B(ind-5),-B(ind-4))*complex(B(ind),0d0)*complex(B(ind+1),B(ind+2)) - &
			complex(C(ind-5),C(ind-4))/complex(C(ind),0d0)* &
			(complex(C(ind-2),-C(ind-1))*complex(B(ind-2),-B(ind-1))*complex(B(ind+1),B(ind+2)) - &
			complex(C(ind+1),-C(ind+2))*R(3,2)))/complex(C(ind-3),0d0)
>>>>>>> added files for complex unitary plus rank one from svn repo
  end if
  
  ! apply diagonal
  ind = 2*(K-1)
<<<<<<< HEAD
  R(2,:) = cmplx(D(ind+1),D(ind+2),kind=8)*R(2,:)
  R(3,:) = cmplx(D(ind+3),D(ind+4),kind=8)*R(3,:)
  
  ! if not at top
  if (K > 1) then
    R(1,:) = cmplx(D(ind-1),D(ind),kind=8)*R(1,:)
  end if

  ! apply Q
  ind = 3*(K-1)  
  R(3,2) = cmplx(Q(ind+4),Q(ind+5),kind=8)*R(3,2)

  H(1,1) = cmplx(Q(ind+1),Q(ind+2),kind=8)
  H(2,1) = cmplx(Q(ind+3),0d0,kind=8)
  H(1,2) = cmplx(-Q(ind+3),0d0,kind=8)
  H(2,2) = cmplx(Q(ind+1),-Q(ind+2),kind=8)
  
=======
  R(2,:) = complex(D(ind+1),D(ind+2))*R(2,:)
  R(3,:) = complex(D(ind+3),D(ind+4))*R(3,:)
  
  ! if not at top
  if (K > 1) then
    R(1,:) = complex(D(ind-1),D(ind))*R(1,:)
  end if

  ! apply Q
  ind = 3*(K-1)	
  R(3,2) = complex(Q(ind+4),Q(ind+5))*R(3,2)

  H(1,1) = complex(Q(ind+1),Q(ind+2))
  H(2,1) = complex(Q(ind+3),0d0)
  H(1,2) = complex(-Q(ind+3),0d0)
  H(2,2) = complex(Q(ind+1),-Q(ind+2))
	
>>>>>>> added files for complex unitary plus rank one from svn repo
  R(2:3,1:2) = matmul(H,R(2:3,1:2))

  ! if not at top
  if (K > 1) then
<<<<<<< HEAD
    H(1,1) = cmplx(Q(ind-2),Q(ind-1),kind=8)
    H(2,1) = cmplx(Q(ind),0d0,kind=8)
    H(1,2) = cmplx(-Q(ind),0d0,kind=8)
    H(2,2) = cmplx(Q(ind-2),-Q(ind-1),kind=8)
      
    R(1:2,1:2) = matmul(H,R(1:2,1:2))
  end if
  
  ! set output
  H = R(2:3,1:2)
=======
    H(1,1) = complex(Q(ind-2),Q(ind-1))
		H(2,1) = complex(Q(ind),0d0)
		H(1,2) = complex(-Q(ind),0d0)
		H(2,2) = complex(Q(ind-2),-Q(ind-1))
			
		R(1:2,1:2) = matmul(H,R(1:2,1:2))
	end if
	
	! set output
	H = R(2:3,1:2)
>>>>>>> added files for complex unitary plus rank one from svn repo

end subroutine ZPFTDB
