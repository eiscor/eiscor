!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_z_urffact_rootsofunity
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Nth order roots of unity by solving a
! corresponding unitary eigenvalue problem using z_urffact_qr.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_z_urffact_rootsofunity

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO
  real(8) :: VV(N)
  complex(8) :: U(N)
  integer :: ITS(N-1)
  
  ! print banner
  print*,""
  print*,"example_z_urffact_rootsofunity:"
  print*,""
  
  ! initialize U and VV to be an upper-Hessenberg permutation matrix
  U = cmplx(0d0,0d0,kind=8)
  U(N) = cmplx(dble((-1)**(N-1)),0d0,kind=8)
  VV = 1d0
  VV(N) = 0d0
  
  ! call z_urffact_qr
  call z_urffact_qr(N,U,VV,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"z_urffact_qr failed."
    print*,"INFO:",INFO
  end if
  
  ! print D
  print*,"Roots computed using z_urffact_qr:"
  do ii=1,N
    print*,U(ii)
  end do
  print*,""
    
end program example_z_urffact_rootsofunity
