!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_z_unifact_knowneigs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine chooses N random eigenvalues on the unit 
! circle, constructs a unitary matrix with prescribed 
! eigenvalues, and computes the N eigenvalues by solving
! a corresponding unitary eigenvalue problem using z_unifact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_z_unifact_knowneigs

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO
  real(8) :: Q(3*(N-1)), D(2*N)
  complex(8) :: Z, he, temp(2,2)
  integer :: ITS(N-1), ind, jj
  double precision :: rm,ro,rp,pi = 3.141592653589793239d0
  double precision :: diff, her
  double precision :: nrm,b1(3),b2(3),b3(3),tb(3),tol=1e-16

  ! real and imag part of eigenvalues
  double precision :: rev(N), iev(N)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print banner
  print*,""
  print*,"example_z_unifact_knowneigs:"
  print*,""
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initialize Q and D to be an identity matrix
  Q = 0d0
  D = 0d0
  do ii=1,n-1
     Q(3*ii-2) = 1d0
     D(2*ii-1) = 1d0
  end do
  D(2*n-1) = 1d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! random eigenvalues, on the unit circle
  do ii = 1,n
     call random_number(rm)
     rev(ii) = sin(2*pi*rm)
     iev(ii) = cos(2*pi*rm)
     ! project on the unit circle
     call d_rot2_vec2gen(rev(ii),iev(ii),rev(ii),iev(ii),nrm)
     D(2*ii-1) = rev(ii)
     D(2*ii)   = iev(ii)
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! inverse eigenvalue problem 
  do ii=n-1,1,-1
     ! generate random bulge
     call random_number(rp)     
     call random_number(ro)     
     call random_number(rm)     
     call z_rot3_vec3gen(rm,ro,rp,b1(1),b1(2),b1(3),NRM)
     ! fuse on the top of the current Q sequence    
     ind = 3*ii-3
     Q(ind+1) = b1(1)
     Q(ind+2) = -b1(2)
     Q(ind+3) = -b1(3)

     ! main chasing loop
     do jj=ii,(n-2)
        ! set indices
        ind = 3*(jj-1)
        ! through diag
        call z_rot3_swapdiag(D(2*jj-1:2*jj+2),b1)
        call z_rot3_turnover(Q((ind+1):(ind+3)),Q((ind+4):(ind+6)),b1)
     end do
     ind = 3*(n-1)
     ! fusion at bottom
     call z_rot3_swapdiag(D((2*n-3):(2*n)),b1)
     call z_unifact_mergebulge(.FALSE.,Q((ind-2):(ind)),D((2*n-3):(2*n)),b1)
  end do
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print*,"Eigenvalues:"
  do ii=1,N
     print*, rev(ii), iev(ii)
  end do
  print*,""
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! call z_unifact_qr
  call z_unifact_qr(.FALSE.,.FALSE.,N,Q,D,1,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"z_unifact_qr failed."
    print*,"INFO:",INFO
  end if
  
  ! print eigenvalues
  print*,"Eigenvalues computed using z_unifact_qr: (real part, imag part, distance to closest exact eigenvalue)"
  do ii=1,N
     diff = sqrt(abs(D(2*ii-1)-rev(1))**2+abs(D(2*ii)-iev(1))**2)
     do jj=2,N
        her = sqrt(abs(D(2*ii-1)-rev(jj))**2+abs(D(2*ii)-iev(jj))**2)
        if (her < diff) then
           diff = her
        end if
     end do
     print*,D(2*ii-1),D(2*ii), diff
  end do
  print*,""

end program example_z_unifact_knowneigs
