!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_orthhess_gaussquad
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the N nodes and weights of the Gauss-Szego 
! rule corresponding to the uniform mesaure on the unit circle. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_orthhess_gaussquad

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO
  real(8) :: WORK(3*N), H(N,N), Z(1,N)
  complex(8) :: E(N), V(1,N)
  integer :: ITS(N-1)
  
  ! print banner
  print*,""
  print*,"example_d_orthhess_gaussquad:"
  print*,""
  
  ! initialize H to be an upper-Hessenberg permutation matrix
  H = 0d0
  do ii=1,N-1
    H(ii+1,ii) = 1d0
  end do
  H(1,N) = 1d0
  
  ! call d_orthhess_qr
  call d_orthhess_qr(.TRUE.,.TRUE.,N,H,WORK,1,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"d_orthhess_qr failed."
    print*,"INFO:",INFO
  end if
  
  ! call d_orthhess_real2complex
  call d_orthhess_real2complex(.TRUE.,N,H,1,Z,E,V,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"d_orthhess_real2complex failed."
    print*,"INFO:",INFO
  end if
  
  ! print nodes
  print*,"Nodes computed using d_orthhess_qr:"
  do ii=1,N
    print*,E(ii)
  end do
  print*,""
  
  ! print weights
  print*,"Nodes computed using d_orthhess_qr:"
  do ii=1,N
    print*,abs(V(1,ii))**2
  end do
  print*,""
  
end program example_d_orthhess_gaussquad
