!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_z_unifact_gaussquad
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the N nodes and weights for a Gauss-Szego
! quadrature rule corresponding to the uniform measure on the
! unit circle (trapezoid rule).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_z_unifact_gaussquad

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO
  real(8) :: Q(3*(N-1)), D(2*N)
  complex(8) :: Z(1,N)
  integer :: ITS(N-1)
  
  ! print banner
  print*,""
  print*,"example_z_unifact_gaussquad:"
  print*,""
  
  ! initialize Q and D to be an upper-hessenberg permutation matrix
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
  call z_unifact_qr(.TRUE.,.TRUE.,N,Q,D,1,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"z_unifact_qr failed."
    print*,"INFO:",INFO
  end if
  
  ! print D
  print*,"Nodes computed using z_unifact_qr:"
  do ii=1,N
    print*,D(2*ii-1),D(2*ii)
  end do
  print*,""

  ! print Z
  print*,"Weights computed using z_unifact_qr:"
  do ii=1,N
    print*,abs(Z(1,ii))**2
  end do
  print*,""
    
end program example_z_unifact_gaussquad
