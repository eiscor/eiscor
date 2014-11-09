!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DEXRUX (Double EXample Roots of Unity eXtreme)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Nth order roots of unity by solving a
! corresponding real orthogonal eigenvalue problem
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program DEXRUX

  implicit none
  
  ! compute variables
  integer :: ii, ind, m, num, N, INFO
  real(8) :: Z(1,1), error, temp, theta, den, t_str, t_stp
  real(8), parameter :: PI = 3.14159265358979323846264338327950d0
  real(8), allocatable :: Q(:), D(:)
  integer, allocatable :: ITS(:)
  
  ! initialize random seed
  call UARIRS(INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    call UARERR(__FILE__,__LINE__,"UARIRS failed",INFO,INFO)
    stop
  end if
  
  ! set number of roots
  N = 2**14
  print*,""
  print*,"number of roots:",N
  
  ! allocate memory
  allocate(Q(2*(N-1)),D(2*N),ITS(N-1))
  
  ! fill Q
  do ii=1,(N-1)
    Q(2*(ii-1)+1) = 0d0
    Q(2*(ii-1)+2) = 1d0       
  end do
  
  ! fill D
  do ii=1,(N-1)
    D(2*(ii-1)+1) = 1d0
    D(2*(ii-1)+2) = 0d0       
  end do
  D(2*(N-1)+1) = (-1d0)**(N-1)
  D(2*(N-1)+2) = 0d0 
  
  ! start timer
  call cpu_time(t_str)
  
  ! call doffqr
  call DOFFQR('N',N,Q,D,Z,ITS,INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    call UARERR(__FILE__,__LINE__,"DOFFQR failed",INFO,INFO)
    stop
  end if
  
  ! stop timer
  call cpu_time(t_stp)
  print*,"time to compute (sec):",t_stp-t_str
  
  ! sort eigenvalues
  call DARSUE('E',N,D,Z,INFO)

  ! check INFO
  if (INFO .NE. 0) then
    call UARERR(__FILE__,__LINE__,"DARSUE failed",INFO,INFO)
    stop
  end if
  
  ! check errors
  error = 0d0
  
  ! root at 1
  temp = abs(cmplx(D(1)-1d0,D(2),kind=8))
  if (temp > error) then
    error = temp
  end if
  
  ! root at -1
  if (mod(N,2).EQ.0) then
    temp = abs(cmplx(D(2*(N-1)+1)+1d0,D(2*(N-1)+2),kind=8))
    if (temp > error) then
      error = temp
    end if
    m = N/2
    num = m-1
    den = dble(m)
  else  
    m = N/2
    num = m
    den = dble(m)+5d-1
  end if
  
  ! conjugate pair roots
  do ii=1,num
    ind = 2+4*(ii-1)
    theta = PI*dble(ii)/den
    temp = abs(cmplx(D(ind+1)-cos(theta),D(ind+2)-sin(theta),kind=8))
    if (temp > error) then
      error = temp
    end if
    temp = abs(cmplx(D(ind+3)-cos(theta),D(ind+4)+sin(theta),kind=8))
    if (temp > error) then
      error = temp
    end if
  end do
  
  print*,"Max error in computed roots:",error
  print*,""
  
  ! free memory
  deallocate(Q,D,ITS)
    
end program DEXRUX
