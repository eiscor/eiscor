!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_extracttri
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine constructs the triangular matrix in the schur 
! decomposition of an upper triangular matrix that is the sum of a 
! unitary matrix and a rank one matrix. This sum is stored as a product 
! of two sequences of Givens' rotations and a complex unimodular 
! diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  JOB             CHARACTER
!                    'E': eigenvalues only
!                    'T': eigenvalues and triangular part
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (3*N)
!                    array of generators for second sequence of Givens' 
!                    rotations, assumed to be the identity matrix
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
!  T              COMPLEX(8) array of dimension (N,N)
!                   if JOB = 'E' unused
!                   if JOB = 'T' stores the upper triangular matrix
!
!  INFO           INTEGER
!                   INFO = 0 implies successful computation.
!                   INFO = 1 implies no convergence in maximum 
!                   allowed iterations.
!                   INFO negative implies that an input variable has
!                   an improper value, i.e. INFO=-2 => N is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_extracttri(JOB,N,Q,D,C,B,T,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: JOB
  integer, intent(in) :: N
  real(8), intent(inout) :: Q(3*N), D(2*(N+1)), C(3*N), B(3*N)
  integer, intent(inout) :: INFO
  complex(8), intent(inout) :: T(N,N)
  
  ! compute variables
  integer :: ii,ind
  complex(8) :: g, p, temp(2,2)
  
  ! initialize info
  INFO = 0
  
  ! if JOB == T
  if (JOB.EQ.'T') then
  
    ! initialize T
    T = cmplx(0d0,0d0,kind=8)
    do ii=1,N
      T(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  
    ! fill with B
    do ii=1,(N-1)
      ind = 3*(ii-1)
      temp(1,1) = cmplx(B(ind+1),B(ind+2),kind=8)
      temp(2,1) = cmplx(B(ind+3),0d0,kind=8)
      temp(1,2) = -temp(2,1)
      temp(2,2) = conjg(temp(1,1))
      T(1:(ii+1),ii:(ii+1)) = matmul(T(1:(ii+1),ii:(ii+1)),temp)
    end do
    T(:,N) = T(:,N)*cmplx(B(3*(N-1)+1),B(3*(N-1)+2),kind=8)
    
    ! update first row
    g = cmplx(1d0,0d0,kind=8)
    p = g
    do ii=1,(N-1)
      ind = 3*(ii-1)
      T(1,ii) = T(1,ii)-(g*cmplx(B(ind+1),B(ind+2),kind=8)+cmplx(C(ind+1),-C(ind+2),kind=8)*B(ind+3)/C(ind+3))/p
      p = p*C(ind+3)
      g = cmplx(B(ind+1),-B(ind+2),kind=8)*cmplx(C(ind+1),-C(ind+2),kind=8)-g*B(ind+3)*C(ind+3)
    end do
    ind = 3*(N-1)
    T(1,N) = T(1,N)-(g*cmplx(B(ind+1),B(ind+2),kind=8)+cmplx(C(ind+1),-C(ind+2),kind=8)*B(ind+3)/C(ind+3))/p
    
    ! apply C
    do ii=1,(N-1)
      ind = 3*(ii-1)
      temp(1,1) = cmplx(C(ind+1),C(ind+2),kind=8)
      temp(2,1) = cmplx(C(ind+3),0d0,kind=8)
      temp(1,2) = -temp(2,1)
      temp(2,2) = conjg(temp(1,1))
      T(ii:(ii+1),(ii+1):N) = matmul(temp,T(ii:(ii+1),(ii+1):N))
      T(ii,ii) = temp(1,1)*T(ii,ii) + temp(1,2)*T(ii+1,ii)
      T(ii+1,ii) = cmplx(0d0,0d0,kind=8)
    end do
    T(N,N) = T(N,N)*cmplx(C(3*(N-1)+1),C(3*(N-1)+2),kind=8) + cmplx(B(3*(N-1)+3),0d0,kind=8)*cmplx(-C(3*(N-1)+3),0d0,kind=8)
    
    ! apply D
    do ii=1,N
      ind = 2*(ii-1)
      T(ii,:) = T(ii,:)*cmplx(D(ind+1),D(ind+2),kind=8)
    end do
  
  end if  
  
  ! update D with eigenvalues
  do ii=1,N
    ind = 3*(ii-1)
    p = cmplx(B(ind+3)/C(ind+3),0d0,kind=8)
    ind = 2*(ii-1)
    p = p*cmplx(D(ind+1),D(ind+2),kind=8)
    D(ind+1) = dble(p)
    D(ind+2) = aimag(p)
  end do
  
end subroutine z_upr1fact_extracttri
