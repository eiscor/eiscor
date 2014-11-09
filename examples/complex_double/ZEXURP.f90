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
  real(8), allocatable :: WORK(:,:)
  complex(8), allocatable :: COEFFS(:), ROOTS(:)
  integer, allocatable :: ITS(:)
  
  ! set number of roots
  N = 2**12
  print*,""
  print*,"number of roots:",N
  
  ! allocate memory
  allocate(WORK(3*N,4),COEFFS(N+1),ROOTS(N),ITS(N-1))
  
  ! call random number generator
  call UARIRS()
  
  ! fill COEFFS
  COEFFS = complex(0d0,0d0)
  COEFFS(1) = complex(1d0,0d0)
  COEFFS(N+1) = complex(-1d0,0d0)
  
  ! start timer
  call cpu_time(t_str)
  
  ! call zmbfrf
  call ZMBFRF(N,COEFFS,ROOTS,WORK,ITS,INFO)
  
  ! check INFO
  if (INFO .NE. 0) then
    write(*,*) "ZMBFRF failed."
    write(*,*) "INFO:",INFO
    write(*,*) ""
    return
  end if
  
  ! stop timer
  call cpu_time(t_stp)
  print*,"time to compute (sec):",t_stp-t_str
  
  ! print ROOTS
!  print*,"ROOTS"
!  do ii=1,N
!    print*,ROOTS(ii)
!  end do
!  print*,""

  ! print ITS
!  print*,"ITS"
!  do ii=1,N-1
!    print*,ITS(ii)
!  end do
!  print*,""
  
  ! free memory
  deallocate(WORK,COEFFS,ROOTS,ITS)
    
end program ZEXRUC
