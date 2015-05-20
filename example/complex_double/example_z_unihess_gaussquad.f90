!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_z_unihess_gaussquad
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the N nodes and weights for a Gauss-Szego
! quadrature rule corresponding to the uniform measure on the
! unit circle (trapezoid rule).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_z_unihess_gaussquad

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO
  real(8) :: WORK(5*N)
  complex(8) :: H(N,N), Z(1,N)
  integer :: ITS(N-1)
  
  ! print banner
  print*,""
  print*,"example_z_unihess_gaussquad:"
  print*,""
  
  ! initialize H to be an upper-Hessenberg permutation matrix
  H = cmplx(0d0,0d0,kind=8)
  do ii=1,N-1
    H(ii+1,ii) = cmplx(1d0,0d0,kind=8)
  end do
  H(1,N) = cmplx(1d0,0d0,kind=8)
  
  ! call z_unihess_qr
  call z_unihess_qr(.TRUE.,.TRUE.,N,H,WORK,1,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"z_unihess_qr failed."
    print*,"INFO:",INFO
  end if
  
  ! print diag of H
  print*,"Nodes computed using z_unihess_qr:"
  do ii=1,N
    print*,dble(H(ii,ii)),aimag(H(ii,ii))
  end do
  print*,""
  
  ! print Z
  print*,"Weights computed using z_unihess_qr:"
  do ii=1,N
    print*,abs(Z(1,ii))**2
  end do
  print*,""
    
end program example_z_unihess_gaussquad
