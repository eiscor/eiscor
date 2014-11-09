!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DEXSLQ (Double EXample Szego-Labatto Quadrature rule)
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
program DEXSLQ

  implicit none
  
  ! compute variables
  integer :: ii, cpair, N, INFO
  real(8) :: c, s, nrm, temp(2,2)
  real(8), allocatable :: Q(:), D(:), Hold(:,:), H(:,:), Z(:,:)
  integer, allocatable :: ITS(:)
  
  ! set number of nodes N
  N = 10
  
  ! allocate memory
  allocate(Q(2*(N-1)),D(2*N),Hold(N,N),H(N,N),Z(N,N),ITS(N-1))
  
  ! fill Q and D
  do ii=1,(N-1)
    s = 1d0/(1d0+ii)
    c = sqrt(1d0-s**2)
    nrm = sqrt(c**2 + s**2)
    Q(2*(ii-1)+1) = c/nrm
    Q(2*(ii-1)+2) = -s/nrm
    D(2*(ii-1)+1) = -1d0
    D(2*(ii-1)+2) = 0d0        
  end do
  D(2*(N-1)+1) = 1d0
  D(2*(N-1)+2) = 0d0
  
  ! initialize H and Hold
  H = 0d0
  do ii=1,N
    H(ii,ii) = 1d0
  end do
  Hold = H
  do ii=1,(N-1)
    temp(1,1) = Q(2*(ii-1)+1)
    temp(2,1) = Q(2*(ii-1)+2)
    temp(1,2) = -temp(2,1)
    temp(2,2) = temp(1,1)
    Hold(:,ii:(ii+1)) = matmul(Hold(:,ii:(ii+1)),temp)
  end do 
  do ii=1,N
    Hold(:,ii) = Hold(:,ii)*D(2*(ii-1)+1)
  end do     
  
  ! call doffqr
  call DOFFQR('I',N,Q,D,Z,ITS,INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "DOFFQR failed."
    write(*,*) "INFO:",INFO
    write(*,*) ""
    stop
  end if
  
  ! print nodes and weights
  print*,"Nodes and Weights"
  cpair = 0
  do ii=1,N
    c = D(2*(ii-1)+1)
    s = D(2*(ii-1)+2)
    if (abs(s).EQ.0d0) then
      print*,cmplx(c,s,kind=8),Z(1,ii)**2
    else if (cpair.EQ.0) then
      print*,cmplx(c,s,kind=8),(Z(1,ii)**2+Z(1,ii+1)**2)/2d0
      print*,cmplx(c,-s,kind=8),(Z(1,ii)**2+Z(1,ii+1)**2)/2d0
      cpair = 1
    else
      cpair = 0
    end if
  end do 
  print*,""  
  
  ! update H
  cpair = 0
  H = 0d0
  do ii=1,N
    c = D(2*(ii-1)+1)
    s = D(2*(ii-1)+2)
    if (abs(s).EQ.0d0) then
      H(ii,ii) = c
    else if (cpair.EQ.0) then
      H(ii,ii) = c
      H(ii,ii+1) = s
      H(ii+1,ii) = -s
      H(ii+1,ii+1) = c
      cpair = 1
    else
      cpair = 0
    end if
  end do 
  
  ! print error
  Hold = matmul(Hold,Z)-matmul(Z,H)
  print*,"Error: max|Hold*Z - Z*H|"
  print*,maxval(abs(Hold))
  print*,""
  
  ! free memory
  deallocate(Q,D,Hold,H,Z,ITS)
    
end program DEXSLQ
