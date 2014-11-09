!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZEXRUX (Zomplex EXample Roots of Unity eXtreme)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Nth order roots of unity by solving a
! corresponding unitary eigenvalue problem
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program ZEXRUX

  implicit none
  
  ! compute variables
  integer :: ii, ind, N, INFO
  real(8) :: Z(1,1), error, temp, theta, den, t_str, t_stp
  real(8), parameter :: PI = 3.14159265358979323846264338327950d0
  real(8), allocatable :: Q(:), D(:)
  integer, allocatable :: ITS(:)
  
  ! set number of roots
  N = 2**10
  print*,""
  print*,"number of roots:",N
  
  ! allocate memory
  allocate(Q(3*(N-1)),D(2*N),ITS(N-1))
  
  ! fill Q
  do ii=1,(N-1)
    Q(3*(ii-1)+1) = 0d0
    Q(3*(ii-1)+2) = 0d0
    Q(3*(ii-1)+3) = 1d0       
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
 
  ! call zuffqr
  call ZUFFQR('N',N,Q,D,Z,ITS,INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "ZUFFQR failed."
    write(*,*) "INFO:",INFO
    write(*,*) ""
    stop
  end if
  
  ! stop timer
  call cpu_time(t_stp)
  print*,"time to compute (sec):",t_stp-t_str
  
  ! sort eigenvalues
  call ZARSUE('E',N,D,Z,INFO)
  if (INFO.NE.0) then
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "INFO:",INFO
    write(*,*) ""
    stop
  end if
  
  ! check errors
  error = 0d0
  
  do ii=1,N
    ind = 2*(ii-1)
    theta = 2d0*PI*dble(ii-1)/dble(N)
    temp = abs(cmplx(D(ind+1)-cos(theta),D(ind+2)-sin(theta),kind=8))
    if (temp > error) then
      error = temp
    end if
  end do
  
  print*,"Max error in computed roots:",error
  print*,""
  
  ! free memory
  deallocate(Q,D,ITS)
    
end program ZEXRUX
