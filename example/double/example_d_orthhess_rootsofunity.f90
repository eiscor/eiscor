!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_orthhess_rootsofunity
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Nth order roots of unity by solving a
! corresponding unitary eigenvalue problem two different ways.
!
! 1) Form upper hessenberg permutation matrix and compute its 
!    eigenvalues using d_orthhess_qr
!
! 2) Construct the factorization directly and compute its 
!    eigenvalues using d_orthfact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_orthhess_rootsofunity

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO, cpair
  real(8) :: WORK(3*N), Q(2*(N-1)), D(N)
  real(8) :: H(N,N), Z
  complex(8) :: E(N), V
  integer :: ITS(N-1)
  
  ! print banner
  print*,""
  print*,"example_d_orthhess_rootsofunity:"
  print*,""
  
  ! initialize H to be an upper hessenberg permutation matrix
  H = 0d0
  do ii=1,N-1
    H(ii+1,ii) = 1d0
  end do
  H(1,N) = 1d0
  
  ! call d_orthhess_qr
  call d_orthhess_qr(.FALSE.,.FALSE.,N,H,WORK,1,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"d_orthhess_qr failed."
    print*,"INFO:",INFO
  end if
  
  ! call d_orthhess_real2complex
  call d_orthhess_real2complex(.FALSE.,N,H,1,Z,E,V,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"d_orthhess_real2complex failed."
    print*,"INFO:",INFO
  end if
  
  ! print E
  print*,"Roots computed using d_orthhess_qr:"
  do ii=1,N
    print*,E(ii)
  end do
  print*,""
  
  ! initialize Q and D to be an upper hessenberg permutation matrix
  do ii=1,N-1
    Q(2*ii-1) = 0d0
    Q(2*ii) = 1d0
  end do
  D = 1d0
  D(N) = (-1d0)**(N-1)
  
  ! call d_orthfact_qr
  call d_orthfact_qr(.FALSE.,.FALSE.,N,Q,D,1,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"d_orthfact_qr failed."
    print*,"INFO:",INFO
  end if
  
  ! call d_orthfact_real2complex
  call d_orthfact_real2complex(.FALSE.,N,Q,D,1,Z,E,V,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"d_orthfact_real2complex failed."
    print*,"INFO:",INFO
  end if
  
  ! print E
  print*,"Roots computed using d_orthfact_qr:"
  do ii=1,N
    print*,E(ii)
  end do
  print*,""
    
end program example_d_orthhess_rootsofunity
