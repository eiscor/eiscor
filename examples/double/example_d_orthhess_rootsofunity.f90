!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_orthhess_rootsofunity
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Nth order roots of unity by solving a
! corresponding unitary eigenvalue problem two different ways.
!
! 1) Form upper hessenberg permutation matrix and compute its 
!    eigenvalues using d_orthhess_qr
!
! 2) Construct the factorization directly and compute its 
!    eigenvalues using d_orthfact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_orthhess_rootsofunity

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO, cpair
  real(8) :: WORK(4*N), Q(2*(N-1)), D(2*N)
  real(8) :: H(N,N), Z
  integer :: ITS(N-1)
  
  ! print banner
  print*,""
  print*,"example_d_orthhess_rootsofunity:"
  print*,""
  
  ! initialize H to be an upper hessenberg permutation matrix
  H = 0d0
  do ii=1,N-1
    H(ii+1,ii) = 1d0
  end do
  H(1,N) = 1d0
  
  ! call dohfqr
  call d_orthhess_qr('N',N,H,Z,ITS,WORK,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"d_orthhess_qr failed."
    print*,"INFO:",INFO
  end if
  
  ! print diag of H
  cpair = 0
  print*,"Roots computed using d_orthhess_qr:"
  do ii=1,N
    if (cpair.EQ.0) then
      if (H(ii,ii+1).EQ.0) then
        print*,H(ii,ii)
      else
        print*,H(ii,ii),H(ii,ii+1)
        print*,H(ii+1,ii+1),H(ii+1,ii)
        cpair = 1
      end if
    else
      cpair = 0
    end if
  end do
  print*,""
  
  ! initialize Q and D to be an upper hessenberg permutation matrix
  do ii=1,N-1
    Q(2*ii-1) = 0d0
    Q(2*ii) = 1d0
  end do
  do ii=1,N-1
    D(2*ii-1) = 1d0
    D(2*ii) = 0d0
  end do
  D(2*N-1) = (-1d0)**(N-1)
  D(2*N) = 0d0
  
  ! call doffqr
  call d_orthfact_qr('N',N,Q,D,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"d_orthfact_qr failed."
    print*,"INFO:",INFO
  end if
  
  ! print D
  print*,"Roots computed using d_orthfact_qr:"
  do ii=1,N
    print*,D(2*ii-1),D(2*ii)
  end do
  print*,""

    
end program example_d_orthhess_rootsofunity
