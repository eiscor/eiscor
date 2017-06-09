!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_orthhess_knowneigs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine chooses N random eigenvalues on the unit 
! circle, constructs a unitary matrix with prescribed 
! eigenvalues, and computes the N eigenvalues by solving
! a corresponding unitary eigenvalue problem using d_orthfact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_orthhess_knowneigs

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO
  real(8) :: Q(2*(N-1)), D(N), Z
  complex(8) :: E(N), V
  integer :: ITS(N-1), ind, jj
  double precision :: rm,rp,pi = 3.141592653589793239d0
  double precision :: he, diff
  double precision :: nrm,b1(2),b2(2),b3(2),tb(2),tol

  ! real and imag part of eigenvalues
  double precision :: rev(N), iev(N) 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print banner
  print*,""
  print*,"example_d_orthhess_knowneigs:"
  print*,""
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initialize Q and D to be an identity matrix
  Q = 0d0
  do ii=1,n-1
     Q(2*ii-1) = 1d0
  end do
  D = 1d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! random eigenvalues, on the unit circle
  do ii = 1,n,2
     call random_number(rm)
     rev(ii) = sin(2*pi*rm)
     iev(ii) = cos(2*pi*rm)
     ! project on the unit circle
     nrm = rev(ii)**2+iev(ii)**2
     if (abs(nrm-1)<tol) then
        nrm = sqrt(nrm)
        rev(ii) = rev(ii)/nrm
        iev(ii) = iev(ii)/nrm
     end if
     Q(2*ii-1) = rev(ii)
     Q(2*ii)   = iev(ii)
     rev(ii+1) = rev(ii)
     iev(ii+1) = -iev(ii)
  end do

  ! print eigenvalues
  print*,"Eigenvalues:"
  do ii=1,N
     print*,rev(ii),iev(ii)
  end do
  print*,""

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! inverse eigenvalue problem 
  do ii=n-3,1,-2
     ! generate random bulge
     call random_number(rp)     
     call random_number(rm)     
     call d_rot2_vec2gen(sqrt(1-rp**2),rp,b1(1),b1(2),NRM)
     call d_rot2_vec2gen(sqrt(1-rm**2),rm,b2(1),b2(2),NRM)
     ! fuse on the top of the current Q sequence
     tb(1) = b2(1)
     tb(2) = -b2(2)
     b3(1) = b1(1)
     b3(2) = -b1(2)
     ind = 2*(ii-1)
     call d_rot2_turnover(tb,b3,Q((ind+1):(ind+2)))
     call d_rot2_fuse(.FALSE.,b3,Q((ind+3):(ind+4)))
     b3 = Q((ind+1):(ind+2))
     Q((ind+1):(ind+2)) = tb
     ! main chasing loop
     do jj=ii,(n-3)
        ind = 2*(jj-1)
        ! set indices
        call d_rot2_turnover(Q((ind+3):(ind+4)),Q((ind+5):(ind+6)),b1)
        call d_rot2_turnover(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b2)
        call d_rot2_turnover(b3,b1,b2)
        ! update bulges
        tb = b2
        b2 = b3
        b3 = b1
        b1 = tb
     end do
     ind = 2*(n-3)
     ! fusion at bottom
     call d_rot2_fuse(.TRUE.,Q((ind+3):(ind+4)),b1)
     call d_rot2_turnover(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b2)
     call d_rot2_fuse(.TRUE.,b3,b2)
     call d_rot2_fuse(.TRUE.,Q((ind+3):(ind+4)),b3)
  end do

  ! call d_orthfact_qr
  call d_orthfact_qr(.FALSE.,.FALSE.,N,Q,D,N,Z,ITS,INFO)
  
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
  
  ! print computed eigenvalues
  print*,"Eigenvalues computed using d_orthfact_qr: (real part, imag part, distance to closest exact eigenvalue)"
  do ii=1,N
     diff = sqrt(abs(dble(E(ii))-rev(1))**2+abs(aimag(E(ii))-iev(1))**2)
     do jj=2,N
        he = sqrt(abs(dble(E(ii))-rev(jj))**2+abs(aimag(E(ii))-iev(jj))**2)
        if (he < diff) then
           diff = he
        end if
     end do
     print*,dble(E(ii)),aimag(E(ii)), diff
  end do
  print*,""
    
end program example_d_orthhess_knowneigs
