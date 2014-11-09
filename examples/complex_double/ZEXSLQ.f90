!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZEXSLQ (Zomplex EXample Szego-Labatto Quadrature rule)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Szego-Labatto quadrature rule based on the
! examples in:
!
! Szego-Labatto quadrature rules. Carl Jagels and Lothar Reichel. 
! Journal of Computational and Applied Mathematics, 200, (2007),
! p. 116-126.
!
! See Example 5.2 p. 124. without the fixed nodes.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program ZEXSLQ

  implicit none
  
  ! compute variables
  integer :: ii, N, INFO
  real(8) :: cr, ci, s, nrm
  complex(8) :: temp(2,2)
  real(8), allocatable :: Q(:), D(:)
  complex(8), allocatable :: Hold(:,:), H(:,:), Z(:,:)
  integer, allocatable :: ITS(:)
  
  ! set number of nodes N
  N = 10
  
  ! allocate memory
  allocate(Q(3*(N-1)),D(2*N),Hold(N,N),H(N,N),Z(N,N),ITS(N-1))
  
  ! fill Q and D
  do ii=1,(N-1)
    s = 1d0/(1d0+ii)
    ci = 0d0
    cr = sqrt(1d0-s**2)
    nrm = sqrt(cr**2 + ci**2 + s**2)
    Q(3*(ii-1)+1) = cr/nrm
    Q(3*(ii-1)+2) = ci/nrm
    Q(3*(ii-1)+3) = -s/nrm
    D(2*(ii-1)+1) = -1d0
    D(2*(ii-1)+2) = 0d0        
  end do
  D(2*(N-1)+1) = 1d0
  D(2*(N-1)+2) = 0d0
  
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
  
  ! print nodes and weights
  print*,"Nodes and Weights"
  do ii=1,N
    print*,H(ii,ii),abs(Z(1,ii))**2
  end do
  print*,""
  
  ! print error
  Hold = matmul(Hold,Z)-matmul(Z,H)
  print*,"Error: max|Hold*Z - Z*H|"
  print*,maxval(abs(Hold))
  print*,""
  
  ! free memory
  deallocate(Q,D,Hold,H,Z,ITS)
    
end program ZEXSLQ
