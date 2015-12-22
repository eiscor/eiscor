#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_symtrid_qr_1dlaplace
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the eigenvalues of the tridiagonal matrix
! [ -1 2 -1 ]. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_symtrid_qr_1dlaplace

  implicit none
  
  ! compute variables
  integer, parameter :: N = 4092
  !integer, parameter :: N = 16
  integer :: ii, jj, INFO
  real(8) :: WORK(14*N), D(N), E(N), Z, t
  !complex(8) :: V
  integer :: ITS(N-1)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  print*,""
  print*,"example_d_orthhess_rootsofunity:"
  print*,""
  
  ! initialize T to be a tridiagonal matrix of the form
  !  2 -1
  ! -1  2 -1
  !     -1 2 ...
  D = 2d0
  E = -1d0
  
  ! call d_orthhess_qr
  call d_symtrid_qr(.FALSE.,.FALSE.,N,D,E,WORK,1,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"d_symtrid_qr failed."
    print*,"INFO:",INFO
  end if
   
  ! print D
  if (N<129) then
     print*,"Roots computed using d_symtrid_qr (and exact roots):"
  end if
  do ii=1,N
     E(ii) = 2d0+2d0*cos(ii*EISCOR_DBL_PI/(N+1d0))
  end do

  Z = 0d0
  do ii=1,N
     t = abs(D(ii)-E(1))
     INFO = 1
     do jj=2,N
        if (abs(D(ii)-E(jj))<t) then
           t = abs(D(ii)-E(jj))
           INFO = jj
        end if
     end do
     
     if (N<129) then
        print*, D(ii), E(INFO), t
     end if
     Z = Z + t**2
  end do
  print*, "forward error (2-norm of the vector) ", dsqrt(Z)
  print*,""
  

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))


  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  print*,""
  print*,"LAPACK"
  print*,""
  
  ! initialize T to be a tridiagonal matrix of the form
  !  2 -1
  ! -1  2 -1
  !     -1 2 ...
  D = 2d0
  E = -1d0

  ! run LAPACK
  call dsteqr ('N', N, D, E, Z, 1, WORK, INFO)

  ! check results
  do ii=1,N
     E(ii) = 2d0+2d0*cos(ii*EISCOR_DBL_PI/(N+1d0))
  end do
  Z = 0d0
  do ii=1,N
     t = abs(D(ii)-E(1))
     INFO = 1
     do jj=2,N
        if (abs(D(ii)-E(jj))<t) then
           t = abs(D(ii)-E(jj))
           INFO = jj
        end if
     end do
     
     if (N<129) then
        print*, D(ii), E(INFO), t
     end if
     Z = Z + t**2
  end do
  print*, "forward error (2-norm of the vector) ", dsqrt(Z)
  print*,""

  ! stop timer
  call system_clock(count=c_stop) 

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program example_d_symtrid_qr_1dlaplace
