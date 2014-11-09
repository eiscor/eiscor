!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DEXURG (Double EXample Uniform Random Generators)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine creates a real orthogonal upper hessenberg matrix by
! constructing generators that have uniformly distributed random 
! coefficients (properly normalized). 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program DEXURG

  implicit none
  
  ! compute variables
  integer :: ii, cpair, N, INFO, rsize
  real(8) :: c, s, nrm, temp(2,2), t_str, t_stp
  real(8), allocatable :: Q(:), D(:), Hold(:,:), H(:,:), Z(:,:)
  integer, allocatable :: ITS(:), seed(:)
  
  ! set number of roots
  N = 2**10
  print*,"number of roots:",N
  print*,""
  
  ! allocate memory
  allocate(Q(2*(N-1)),D(2*N),Hold(N,N),H(N,N),Z(N,N),ITS(N-1))
  
  ! initialize random seed
  call UARIRS(INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    call UARERR(__FILE__,__LINE__,"UARIRS failed",INFO,INFO)
    stop
  end if
  
  ! fill Q
  do ii=1,(N-1)
    call random_number(c)
    call random_number(s)
    nrm = sqrt(c**2 + s**2)
    if (nrm.EQ.0d0) then
      Q(2*(ii-1)+1) = 0d0
      Q(2*(ii-1)+2) = 1d0
    else
      Q(2*(ii-1)+1) = c/nrm
      Q(2*(ii-1)+2) = s/nrm
    end if          
  end do
  
  ! fill D
  do ii=1,N
    call random_number(c)
    D(2*(ii-1)+1) = sign(1d0,c)
    D(2*(ii-1)+2) = 0d0       
  end do
  
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
  
  
  ! print Hold
!  print*,"Hold"
!  do ii=1,N
!    print*,Hold(ii,:)
!  end do
!  print*,""

  ! start timer
  call cpu_time(t_str)
  
  ! call doffqr
  call DOFFQR('I',N,Q,D,Z,ITS,INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    call UARERR(__FILE__,__LINE__,"DOFFQR failed",INFO,INFO)
    stop
  end if
  
  ! stop timer
  call cpu_time(t_stp)
  print*,"time to compute (sec):",t_stp-t_str
  
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
 
  ! print H
!  print*,"H"
!  do ii=1,N
!    print*,H(ii,:)
!  end do
!  print*,""
  
  ! print Z
!  print*,"Z"
!  do ii=1,N
!    print*,Z(ii,:)
!  end do
!  print*,""
  
  ! print error
  Hold = matmul(Hold,Z)-matmul(Z,H)

  print*,"Error: max|Hold*Z - Z*H|"
  print*,maxval(abs(Hold))
  print*,""
  
  ! free memory
  deallocate(Q,D,Hold,H,Z,ITS)
  !deallocate(seed)
    
end program DEXURG
