!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_orthfact_rootsofunity
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Nth order roots of unity by solving a
! corresponding unitary eigenvalue problem using d_orthfact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_orthfact_rootsofunity

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO
  real(8) :: Q(2*(N-1)), D(N), Z
  complex(8) :: E(N), V
  integer :: ITS(N-1)
  
  ! print banner
  print*,""
  print*,"example_d_orthfact_rootsofunity:"
  print*,""
  
  ! initialize Q and D to be an upper-Hessenberg permutation matrix
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
    
end program example_d_orthfact_rootsofunity
