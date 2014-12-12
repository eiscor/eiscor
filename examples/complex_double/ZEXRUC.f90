!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZEXRUC (Zomplex EXample Roots of Unity via Companion matrix)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Nth order roots of unity by solving a
! corresponding companion matrix eigenvalue problem
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program ZEXRUC

  implicit none
  
  ! compute variables
  integer :: ii, ind, N, INFO
  real(8) :: error, temp, theta, den, t_str, t_stp
  real(8) :: cr, ci, s, nrm
  real(8), parameter :: PI = 3.14159265358979323846264338327950d0
  real(8), allocatable :: Q(:), D(:), C(:), B(:)
  complex(8), allocatable :: COEFFS(:), T(:,:), Z(:,:), H(:,:)
  integer, allocatable :: ITS(:)
  
  ! set number of roots
  N = 2**10
  print*,""
  print*,"number of roots:",N
  
  ! allocate memory
  allocate(Q(3*N),D(2*(N+1)),C(3*N),B(3*N),ITS(N-1),COEFFS(N),Z(N,N),T(N,N),H(N,N))
  
  ! call random number generator
  call UARIRS()
  
  ! fill COEFFS
  COEFFS = cmplx(0d0,0d0,kind=8)
  COEFFS(N) = cmplx(-1d0,0d0,kind=8)
  
  ! build companion matrix
  H = cmplx(0d0,0d0,kind=8)
  do ii=1,(N-1)
    H(ii+1,ii) = cmplx(1d0,0d0,kind=8)
  end do
  do ii=1,N
    H(ii,N) = -COEFFS(N+1-ii)
  end do
  
  ! construct factors
  call ZMBCQR('I',N,COEFFS,Q,D,C,B,Z,INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "ZMBCQR failed."
    write(*,*) "INFO:",INFO
    write(*,*) ""
  end if
  
  ! print Q
!  print*,"Q"
!  do ii=1,N
!    print*,Q(3*(ii-1)+1),Q(3*(ii-1)+2)
!    print*,Q(3*(ii-1)+3)  
!  end do
!  print*,"" 
 
  ! print D
!  print*,"D"
!  do ii=1,(N+1)
!    print*,D(2*(ii-1)+1),D(2*(ii-1)+2)
!  end do
!  print*,""  
  
  ! print C
!  print*,"C"
!  do ii=1,N
!    print*,C(3*(ii-1)+1),C(3*(ii-1)+2)
!    print*,C(3*(ii-1)+3)  
!  end do
!  print*,"" 
  
  ! print B
!  print*,"B"
!  do ii=1,N
!    print*,B(3*(ii-1)+1),B(3*(ii-1)+2)
!    print*,B(3*(ii-1)+3)  
!  end do
!  print*,"" 
  
  ! compute upper-triangular part
!  call ZPFFET('T',N,Q,D,C,B,T,INFO)
  
  ! check INFO
!  if (INFO .NE. 0) then
!    write(*,*) "ZPFFET failed."
!    write(*,*) "INFO:",INFO
!    write(*,*) ""
!    return
!  end if
  
  ! print T
!  print*,"T"
!  do ii=1,N
!    print*,T(ii,:) 
!  end do
!  print*,"" 
  
!pause
  
  ! start timer
  call cpu_time(t_str)
  
  ! call zpffqr
  call ZPFFQR('V',N,Q,D,C,B,Z,ITS,INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "ZPFFQR failed."
    write(*,*) "INFO:",INFO
    write(*,*) ""
  end if
  
  ! stop timer
  call cpu_time(t_stp)
  print*,"time to compute (sec):",t_stp-t_str
  
  ! print Q
!  print*,"Q"
!  do ii=1,N
!    print*,Q(3*(ii-1)+1),Q(3*(ii-1)+2)
!    print*,Q(3*(ii-1)+3)  
!  end do
!  print*,"" 
 
  ! print D
!  print*,"D"
!  do ii=1,(N+1)
!    print*,D(2*(ii-1)+1),D(2*(ii-1)+2)
!  end do
!  print*,""  
  
  ! print C
!  print*,"C"
!  do ii=1,N
!    print*,C(3*(ii-1)+1),C(3*(ii-1)+2)
!    print*,C(3*(ii-1)+3)  
!  end do
!  print*,"" 
  
  ! print B
!  print*,"B"
!  do ii=1,N
!    print*,B(3*(ii-1)+1),B(3*(ii-1)+2)
!    print*,B(3*(ii-1)+3)  
!  end do
!  print*,"" 
  
  ! compute upper-triangular part
  call ZPFFET('T',N,Q,D,C,B,T,INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "ZPFFET failed."
    write(*,*) "INFO:",INFO
    write(*,*) ""
  end if
  
  ! compute residual
  H = matmul(H,Z)-matmul(Z,T)
  
  ! print max abs error
  print*,"max abs error:",maxval(abs(H))
  print*,""
  
  ! print H
!  print*,"H"
!  do ii=1,N
!    print*,abs(H(ii,:)) 
!  end do
!  print*,"" 
  
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
  deallocate(Q,D,C,B,COEFFS,Z,T,ITS,H)
    
end program ZEXRUC
