!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZEXRUX2 (Zomplex EXample Roots of Unity eXtreme 2)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Nth order roots of unity by solving a
! corresponding unitary eigenvalue problem
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program ZEXRUX2

  implicit none
  
  ! compute variables
  integer :: ii, ind, N, INFO
  real(8) :: error, temp, theta, den, t_str, t_stp
  real(8) :: cr, ci, s, nrm
  real(8), parameter :: PI = 3.14159265358979323846264338327950d0
  real(8), allocatable :: Q(:), D(:), C(:), B(:)
  complex(8), allocatable :: T(:,:), Z(:,:)
  integer, allocatable :: ITS(:)
  
  ! set number of roots
  N = 3
  print*,""
  print*,"number of roots:",N
  
  ! allocate memory
  allocate(Q(3*N),D(2*(N+1)),C(3*N),B(3*N),ITS(N-1),Z(N,N),T(N,N))
  
  ! call random number generator
  call UARIRS()
  
  ! fill Q
  do ii=1,(N-1)
    call random_number(cr)
    call random_number(ci)
    call random_number(s)
    nrm = sqrt(cr*cr+ci*ci+s*s)
    Q(3*(ii-1)+1) = cr/nrm
    Q(3*(ii-1)+2) = ci/nrm
    Q(3*(ii-1)+3) = s/nrm    
  end do
  Q(3*(N-1)+1) = 1d0
  Q(3*(N-1)+2) = 0d0
  Q(3*(n-1)+3) = 0d0  
  
  ! fill D
  do ii=1,(N+1)
    call random_number(cr)
    call random_number(ci)
    nrm = sqrt(cr*cr+ci*ci)
    D(2*(ii-1)+1) = cr/nrm
    D(2*(ii-1)+2) = ci/nrm      
  end do
  
  ! fill C
  do ii=1,N
    call random_number(cr)
    call random_number(ci)
    call random_number(s)
    nrm = sqrt(cr*cr+ci*ci+s*s)
    C(3*(ii-1)+1) = cr/nrm
    C(3*(ii-1)+2) = ci/nrm
    C(3*(ii-1)+3) = s/nrm    
  end do
  
  ! fill B
  do ii=1,N
    call random_number(cr)
    call random_number(ci)
    call random_number(s)
    nrm = sqrt(cr*cr+ci*ci+s*s)
    B(3*(ii-1)+1) = cr/nrm
    B(3*(ii-1)+2) = ci/nrm
    B(3*(ii-1)+3) = s/nrm    
  end do
  
  ! print Q
  print*,"Q"
  do ii=1,N
    print*,Q(3*(ii-1)+1),Q(3*(ii-1)+2)
    print*,Q(3*(ii-1)+3)  
  end do
  print*,"" 
 
  ! print D
  print*,"D"
  do ii=1,(N+1)
    print*,D(2*(ii-1)+1),D(2*(ii-1)+2)
  end do
  print*,""  
  
  ! print C
  print*,"C"
  do ii=1,N
    print*,C(3*(ii-1)+1),C(3*(ii-1)+2)
    print*,C(3*(ii-1)+3)  
  end do
  print*,"" 
  
  ! print B
  print*,"B"
  do ii=1,N
    print*,B(3*(ii-1)+1),B(3*(ii-1)+2)
    print*,B(3*(ii-1)+3)  
  end do
  print*,"" 
  
  ! start timer
  call cpu_time(t_str)
  
  ! call zpffqr
  call ZPFFQR('I',N,Q,D,C,B,Z,ITS,INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "ZPFFQR failed."
    write(*,*) "INFO:",INFO
    write(*,*) ""
    return
  end if
  
  ! stop timer
  call cpu_time(t_stp)
  print*,"time to compute (sec):",t_stp-t_str
  
  ! print Q
  print*,"Q"
  do ii=1,N
    print*,Q(3*(ii-1)+1),Q(3*(ii-1)+2)
    print*,Q(3*(ii-1)+3)  
  end do
  print*,"" 
 
  ! print D
  print*,"D"
  do ii=1,(N+1)
    print*,D(2*(ii-1)+1),D(2*(ii-1)+2)
  end do
  print*,""  
  
  ! print C
  print*,"C"
  do ii=1,N
    print*,C(3*(ii-1)+1),C(3*(ii-1)+2)
    print*,C(3*(ii-1)+3)  
  end do
  print*,"" 
  
  ! print B
  print*,"B"
  do ii=1,N
    print*,B(3*(ii-1)+1),B(3*(ii-1)+2)
    print*,B(3*(ii-1)+3)  
  end do
  print*,"" 
  
  ! print Z
  print*,"Z"
  do ii=1,N
    print*,Z(ii,:) 
  end do
  print*,"" 
  
  ! compute upper-triangular part
  call ZPFFET('T',N,Q,D,C,B,T,INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "ZPFFET failed."
    write(*,*) "INFO:",INFO
    write(*,*) ""
    return
  end if
  
  ! print T
  print*,"T"
  do ii=1,N
    print*,T(ii,:) 
  end do
  print*,"" 
  
  ! check errors
!  error = 0d0
  
!  do ii=1,N
!    ind = 2*(ii-1)
!    theta = 2d0*PI*dble(ii-1)/dble(N)
!    temp = abs(complex(D(ind+1)-cos(theta),D(ind+2)-sin(theta)))
!    if (temp > error) then
!      error = temp
!    end if
!  end do
  
!  print*,"Max error in computed roots:",error
!  print*,""
  
  ! free memory
  deallocate(Q,D,C,B,Z,T,ITS)
    
end program ZEXRUX2
