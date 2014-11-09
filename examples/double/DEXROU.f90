!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DEXROU (Double EXample Roots Of Unity)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Nth order roots of unity by solving a
! corresponding real orthogonal eigenvalue problem
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program DEXROU

  implicit none
  
  ! compute variables
  integer :: ii, N, INFO
  real(8), allocatable :: WORK(:), Hold(:,:), H(:,:), Z(:,:)
  integer, allocatable :: ITS(:)
  
  ! set number of roots
  N = 5
  
  ! allocate memory
  allocate(WORK(4*N),Hold(N,N),H(N,N),Z(N,N),ITS(N-1))
  
  ! initialize H
  H = 0d0
  do ii=1,N-1
    H(ii+1,ii) = 1d0
  end do
!  H(1,N) = (-1d0)**(N-1)
  
  ! store in Hold
  Hold = H
  
  ! call dohfqr
  call DOHFQR('I',N,H,Z,ITS,WORK,INFO)
  
  ! print H
  print*,"H"
  do ii=1,N
    print*,H(ii,:)
  end do
  print*,""
  
  ! print Z
  print*,"Z"
  do ii=1,N
    print*,Z(ii,:)
  end do
  print*,""
  
  ! print error
  Hold = abs(matmul(Hold,Z)-matmul(Z,H))
  print*,"Error: Hold*Z - Z*H"
  do ii=1,N
    print*,Hold(ii,:)
  end do
  print*,""
  
  ! free memory
  deallocate(WORK,Hold,H,Z,ITS)
    
end program DEXROU
