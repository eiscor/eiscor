!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_factorcompmat
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a compressed QZ factorization for the 
! companion matrix of a polynomial expressed in the monomial basis. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  COEFFS          COMPLEX(8) array of dimension (N)
!                   coefficients of polynomial, assumed to be of degree
!                   exactly N, have leading coefficient 1 and have no 
!                   zero roots
!
!  Q               REAL(8) array of dimension (3*N)
!                    array of generators for first sequence of Givens' 
!                    rotations
!
!  D               REAL(8) array of dimension (2*(N+1))
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
!  INFO           INTEGER
!                   INFO = 0 implies successful computation.
!                   INFO negative implies that an input variable has
!                   an improper value, i.e. INFO=-2 => N is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_factorcompmat(COMPZ,N,COEFFS,Q,D,C,B,Z,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: COMPZ
  integer, intent(in) :: N
  complex(8), intent(in) :: COEFFS(N)
  real(8), intent(inout) :: Q(3*N), D(2*(N+1)), C(3*N), B(3*N)
  complex(8), intent(inout) :: Z(N,N)
  integer, intent(inout) :: INFO
  
  ! compute variables
  integer :: ii, ind
  real(8) :: nrm
  complex(8) :: t1, t2, phase
  
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
  call IARNAN(N,INFO)
  if (INFO.NE.0) then
    INFO = -2
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "N contains a NAN."
    write(*,*) ""
    return
  end if
  call IARINF(N,INFO)
  if (INFO.NE.0) then
    INFO = -2
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "N contains an INF."
    write(*,*) ""
    return
  end if
  if (N < 2) then
    INFO = -2
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "N must be at least 2."
    write(*,*) ""
    return
  end if
  
  ! check COEFFS
  call ZARACH1(N,COEFFS,INFO)
  if (INFO.NE.0) then
    INFO = -3
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "COEFFS contains an INF or NAN"
    write(*,*) ""
    return
  end if
  if (abs(COEFFS(N)).EQ.0d0) then
    INFO = -3
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "COEFFS has a zero root."
    write(*,*) ""
    return
  end if

  ! check Z
  if (COMPZ.EQ.'V') then
    call ZARACH2(N,N,Z,INFO)
    if (INFO.NE.0) then
      INFO = -8
      write(*,*) "Error in "//__FILE__//" line:",__LINE__
      write(*,*) "Z contains an INF or NAN"
      write(*,*) ""
      return
    end if
  end if      
  
  ! compute the phase of COEFFS(N)
  phase = COEFFS(N)
  nrm = abs(phase)
  phase = phase/nrm

  ! initialize Z
  if (COMPZ.EQ.'I') then
    Z = cmplx(0d0,0d0,kind=8)
    do ii=1,N
      Z(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if
  
  ! update Z
  if (COMPZ.NE.'N') then
    Z(:,N) = conjg(phase)*Z(:,N)
  end if
  
  ! set Q
  do ii=1,(N-1)
    ind = 3*(ii-1)
    Q(ind+1) = 0d0
    Q(ind+2) = 0d0
    Q(ind+3) = 1d0
  end do
  ind = 3*(N-1)
  Q(ind+1) = 1d0
  Q(ind+2) = 0d0
  Q(ind+3) = 0d0
  
  ! set D
  do ii=1,N+1
    ind = 2*(ii-1)
    D(ind+1) = 1d0
    D(ind+2) = 0d0
  end do
  ind = 2*(N-2)
  D(ind+1) = dble(phase)
  D(ind+2) = aimag(phase)
  
  ! initialize B and C
  t1 = cmplx(nrm*(-1d0)**(N),0d0,kind=8)
  t2 = cmplx((-1d0)**(N-1),0d0,kind=8)
  ind = 3*(N-1)
  call DARCG22(nrm*(-1d0)**(N),(-1d0)**(N-1),C(ind+1),C(ind+3),nrm,INFO)
  C(ind+2) = 0d0
  B(ind+1) = C(ind+3)*(-1d0)**(N)
  B(ind+2) = 0d0
  B(ind+3) = C(ind+1)*(-1d0)**(N)
  
  do ii=2,N
    ind = 3*(N-ii+1)
    t2 = cmplx(C(ind+1),-C(ind+2),kind=8)*t1 + cmplx(C(ind+3),0d0,kind=8)*t2
    t1 = -COEFFS(ii-1)*conjg(phase)
    ind = 3*(N-ii)
    call ZARCG43(dble(t1),aimag(t1),dble(t2),aimag(t2),C(ind+1),C(ind+2),C(ind+3),nrm,INFO)
    B(ind+1) = C(ind+1)
    B(ind+2) = -C(ind+2)
    B(ind+3) = -C(ind+3)
  end do

end subroutine z_upr1fact_factorcompmat
