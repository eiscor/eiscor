!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZEXROU (Zomplex EXample Roots Of Unity)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Nth order roots of unity by solving a
! corresponding unitary eigenvalue problem
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program ZEXROU

  implicit none
  
  ! compute variables
  integer :: ii, N, INFO
  real(8), allocatable :: WORK(:)
  complex(8), allocatable :: Hold(:,:), H(:,:), Z(:,:)
  integer, allocatable :: ITS(:)
  
  ! set number of roots
  N = 3
  
  ! allocate memory
  allocate(WORK(5*N),Hold(N,N),H(N,N),Z(N,N),ITS(N-1))
  
  ! initialize H
  H = cmplx(0d0,0d0,kind=8)
  do ii=1,N-1
    H(ii+1,ii) = cmplx(1d0,0d0,kind=8)
  end do
  H(1,N) = cmplx((-1d0)**(N-1),0d0,kind=8)
  
  ! store in Hold
  Hold = H
  
  ! call zuhfqr
  call ZUHFQR('I',N,H,Z,ITS,WORK,INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "ZUHFQR failed."
    write(*,*) "INFO:",INFO
    write(*,*) ""
    stop
  end if
  
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
  Hold = matmul(Hold,Z)-matmul(Z,H)
  print*,"Error: max|Hold*Z - Z*H|"
  print*,maxval(abs(Hold))
  print*,""
  
  ! free memory
  deallocate(WORK,Hold,H,Z,ITS)
    
end program ZEXROU
