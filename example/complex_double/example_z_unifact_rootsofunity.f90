!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_z_unifact_rootsofunity
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Nth order roots of unity by solving a
! corresponding unitary eigenvalue problem using z_unifact_qr.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_z_unifact_rootsofunity

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO
  real(8) :: Q(3*(N-1)), D(2*N)
  complex(8) :: Z
  integer :: ITS(N-1)
  
  ! print banner
  print*,""
  print*,"example_z_unifact_rootsofunity:"
  print*,""
  
  ! initialize Q and D to be an upper-Hessenberg permutation matrix
  do ii=1,N-1
    Q(3*ii-2) = 0d0
    Q(3*ii-1) = 0d0
    Q(3*ii) = 1d0
  end do
  do ii=1,N-1
    D(2*ii-1) = 1d0
    D(2*ii) = 0d0
  end do
  D(2*N-1) = (-1d0)**(N-1)
  D(2*N) = 0d0
  
  ! call z_unifact_qr
  call z_unifact_qr(.FALSE.,.FALSE.,N,Q,D,N,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"z_unifact_qr failed."
    print*,"INFO:",INFO
  end if
  
  ! print D
  print*,"Roots computed using z_unifact_qr:"
  do ii=1,N
    print*,D(2*ii-1),D(2*ii)
  end do
  print*,""
    
end program example_z_unifact_rootsofunity
