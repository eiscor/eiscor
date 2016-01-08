!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_z_unifact_backwarderror
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_z_unifact_backwarderror

  implicit none
  
  ! compute variables
  integer, parameter :: N = 512
  integer :: ii, INFO
  real(8) :: Q(3*(N-1)), D(2*N)
  complex(8) :: Z(N,N), he, temp(2,2)
  integer :: ITS(N-1), ind, jj
  double precision :: rm,ro,rp,pi = 3.141592653589793239d0
  double precision :: diff, her
  double precision :: nrm,b1(3),b2(3),b3(3),tb(3),tol=1e-16

  ! real and imag part of eigenvalues
  double precision :: rev(N), iev(N)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print banner
  print*,""
  print*,"example_z_unifact_backwarderror:"
  print*,""
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initialize Q and D to be an identity matrix

  open (unit=7, file="QD 512.txt", status='unknown', position="rewind")
  do ii=1,N-1
     read (7,*) Q(3*ii-2), Q(3*ii-1), Q(3*ii)
!     print*, ii, Q(3*ii-2), Q(3*ii-1), Q(3*ii)  
  end do
  do ii=1,N
     read (7,*) D(2*ii-1), D(2*ii)
!     print*, ii, D(2*ii-1), D(2*ii)  
  end do
  close(7)


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! call z_unifact_qr
  call z_unifact_qr(.TRUE.,.TRUE.,N,Q,D,N,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"z_unifact_qr failed."
    print*,"INFO:",INFO
  end if
  

end program example_z_unifact_backwarderror
