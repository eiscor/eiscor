!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_orthfact_gaussquad
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the N nodes and weights for Gauss-Szego rule
! corresponding to the unform measure on the unit circle.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_orthfact_gaussquad

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO
  real(8) :: Q(2*(N-1)), D(N), Z(1,N)
  complex(8) :: E(N), V(1,N)
  integer :: ITS(N-1)
  
  ! print banner
  print*,""
  print*,"example_d_orthfact_gaussquad:"
  print*,""
  
  ! initialize Q and D to be an upper-Hessenberg permutation matrix
  do ii=1,N-1
    Q(2*ii-1) = 0d0
    Q(2*ii) = 1d0
  end do
  D = 1d0
  D(N) = sign((-1d0)**(N-1),1d0)
  
  ! call d_orthfact_qr
  call d_orthfact_qr(.TRUE.,.TRUE.,N,Q,D,1,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"d_orthfact_qr failed."
    print*,"INFO:",INFO
  end if

  ! call d_orthfact_real2complex
  call d_orthfact_real2complex(.TRUE.,N,Q,D,1,Z,E,V,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"d_orthfact_real2complex failed."
    print*,"INFO:",INFO
  end if
  
  ! print nodes
  print*,"Nodes computed using d_orthfact_qr:"
  do ii=1,N
    print*,E(ii)
  end do
  print*,""
    
  ! print weights
  print*,"Weights computed using d_orthfact_qr:"
  do ii=1,N
    print*,abs(V(1,ii))**2
  end do
  print*,""
    
end program example_d_orthfact_gaussquad
