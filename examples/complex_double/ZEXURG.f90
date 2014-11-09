!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZEXURG (Zomplex EXample Uniform Random Generators)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Nth order roots of unity by solving a
! corresponding unitary eigenvalue problem
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program ZEXURG

  implicit none
  
  ! compute variables
  integer :: ii, N, INFO
  real(8) :: cr, ci, s, nrm
  complex(8) :: temp(2,2)
  real(8), allocatable :: Q(:), D(:)
  complex(8), allocatable :: Hold(:,:), H(:,:), Z(:,:)
  integer, allocatable :: ITS(:)
  
  ! set number of roots
  N = 3
  
  ! allocate memory
  allocate(Q(3*(N-1)),D(2*N),Hold(N,N),H(N,N),Z(N,N),ITS(N-1))
  
  ! initialize random seed
  call UARIRS(INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    call UARERR(__FILE__,__LINE__,"UARIRS failed",INFO,INFO)
    stop
  end if
  
  ! fill Q
  do ii=1,(N-1)
    call random_number(cr)
    call random_number(ci)
    call random_number(s)
    nrm = sqrt(cr**2 + ci**2 + s**2)
    if (nrm.EQ.0d0) then
      Q(3*(ii-1)+1) = 0d0
      Q(3*(ii-1)+2) = 0d0
      Q(3*(ii-1)+3) = 1d0
    else
      Q(3*(ii-1)+1) = cr/nrm
      Q(3*(ii-1)+2) = ci/nrm
      Q(3*(ii-1)+3) = s/nrm
    end if          
  end do
  
  ! fill D
  do ii=1,N
    call random_number(cr)
    call random_number(ci)
    nrm = sqrt(cr**2 + ci**2)
    if (nrm.EQ.0d0) then
      D(2*(ii-1)+1) = 0d0
      D(2*(ii-1)+2) = 1d0
    else
      D(2*(ii-1)+1) = cr/nrm
      D(2*(ii-1)+2) = ci/nrm
    end if          
  end do
  
  ! initialize H and Hold
  H = cmplx(0d0,0d0,kind=8)
  do ii=1,N
    H(ii,ii) = cmplx(1d0,0d0,kind=8)
  end do
  Hold = H
  
  do ii=1,(N-1)
    temp(1,1) = cmplx(Q(3*(ii-1)+1),Q(3*(ii-1)+2),kind=8)
    temp(2,1) = cmplx(Q(3*(ii-1)+3),0d0,kind=8)
    temp(1,2) = -temp(2,1)
    temp(2,2) = conjg(temp(1,1))
    Hold(:,ii:(ii+1)) = matmul(Hold(:,ii:(ii+1)),temp)
  end do 
  do ii=1,N
    Hold(:,ii) = Hold(:,ii)*cmplx(D(2*(ii-1)+1),D(2*(ii-1)+2),kind=8)
  end do     
  
  ! call zuffqr
  call ZUFFQR('I',N,Q,D,Z,ITS,INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "ZUFFQR failed."
    write(*,*) "INFO:",INFO
    write(*,*) ""
    stop
  end if
  
  ! update H
  do ii=1,N
    H(ii,ii) = cmplx(D(2*(ii-1)+1),D(2*(ii-1)+2),kind=8)
  end do  
  
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
  deallocate(Q,D,Hold,H,Z,ITS)
    
end program ZEXURG
