!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_orthhess_rootsofunity
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Nth order roots of unity by solving a
! corresponding unitary eigenvalue problem using d_orthhess_qr.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_orthhess_rootsofunity

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO, cpair
  real(8) :: WORK(3*N), H(N,N), Z
  complex(8) :: E(N), V
  integer :: ITS(N-1)
  
  ! print banner
  print*,""
  print*,"example_d_orthhess_rootsofunity:"
  print*,""
  
  ! initialize H to be an upper-Hessenberg permutation matrix
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
  
end program example_d_orthhess_rootsofunity
