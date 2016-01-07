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
  integer, parameter :: N = 6
  !integer, parameter :: N = 16
  !integer, parameter :: N = 200
  !integer, parameter :: N = 6184
  !integer, parameter :: N = 4092
  !integer, parameter :: N = 2048
  !integer, parameter :: N = 1024
  !integer, parameter :: N = 512
  real(8), parameter :: scale = 1d0
  integer :: ii, jj, kk, INFO, IWORK(3+5*N)
  real(8) :: WORK(14*N), D(N), E(N), Z, t, t2
  real(8) :: Ds(N), Es(N)
  !complex(8) :: V
  integer :: ITS(N-1)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  print*,""
  print*,"example_d_symtrid_qr_1dlaplace:"
  print*,""
  print*, "N ", N, " scale ", scale
  print*,""

  ! initialize T to be a tridiagonal matrix of the form
  !  2 -1
  ! -1  2 -1
  !     -1 2 ...
  Ds = 0d0
  Es = -5d-1
  Es = scale*Es

  !call u_randomseed_initialize(INFO)
  if (1==0) then
     do ii=1,N 
        call random_number(t)
        Ds(ii) = t
        call random_number(t)
        Es(ii) = t
     end do
  end if

  ! symtrid_qr
  D = Ds
  E = Es
 
  ! call d_orthhess_qr
  call d_symtrid_qr(.FALSE.,.FALSE.,N,D,E,WORK,1,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"d_symtrid_qr failed."
    print*,"INFO:",INFO
  end if
  ! stop timer
  call system_clock(count=c_stop)
   
  ! print D
  if (N<129) then
     print*,"Roots computed using d_symtrid_qr (and exact roots):"
  end if
  do ii=1,N
     E(ii) = 0d0+scale*cos(ii*EISCOR_DBL_PI/(N+1d0))
  end do

  Z = 0d0
  t2 = 0d0
  do ii=1,N
     t = abs(D(ii)-E(1))
     INFO = 1
     do jj=2,N
        if (abs(D(ii)-E(jj))<t) then
           t = abs(D(ii)-E(jj))
           INFO = jj
        end if
     end do
     
     if (N<519) then
        print*, D(ii), E(INFO), t
     end if
     if (t.GT.t2) then
        t2 = t
        kk = INFO
     end if
     Z = Z + t**2
  end do
  print*, "forward error (2-norm of the vector) ", dsqrt(Z)
  if ((dsqrt(Z)>1d-3).OR.(Z.NE.Z)) then
     call u_test_failed(__LINE__)
  end if
  print*,""
  print*, "largest forward error for eigenvalue ", kk, " error ",t2," eigenvalue ", E(kk)
  print*,""

    
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))


  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  print*,""
  print*,"LAPACK DSTEQR"
  print*,""
  
  ! initialize T to be a tridiagonal matrix of the form
  !  2 -1
  ! -1  2 -1
  !     -1 2 ...
  !D = 0d0
  !E = -5d-1

  D = Ds
  E = Es

  ! run LAPACK
  call dsteqr ('N', N, D, E, Z, 1, WORK, INFO)
  ! stop timer
  call system_clock(count=c_stop) 

  ! check results
  do ii=1,N
     E(ii) = 0d0+scale*cos(ii*EISCOR_DBL_PI/(N+1d0))
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


  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))


  ! print banner
  print*,""
  print*,"LAPACK DSTEVD"
  print*,""
  
  ! initialize T to be a tridiagonal matrix of the form
  !  2 -1
  ! -1  2 -1
  !     -1 2 ...
  D = Ds
  E = Es
  !D = 0d0
  !E = -5d-1

  call system_clock(count=c_start)

  ! run DSTEVD
  call dstevd ('N', N, D, E, Z, 1, WORK, 14*N, IWORK, 3+5*N, INFO) 
  ! stop timer
  call system_clock(count=c_stop) 

  ! check results
  do ii=1,N
     E(ii) = 0d0+scale*cos(ii*EISCOR_DBL_PI/(N+1d0))
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

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

  ! print banner
  print*,""
  print*,"LAPACK DSTERF"
  print*,""
  
  ! initialize T to be a tridiagonal matrix of the form
  !  2 -1
  ! -1  2 -1
  !     -1 2 ...
  !D = 0d0
  !E = -5d-1

  D = Ds
  E = Es

  call system_clock(count=c_start)

  ! run LAPACK
  call dsterf (N, D, E, INFO)

  ! stop timer
  call system_clock(count=c_stop) 

  ! check results
  do ii=1,N
     E(ii) = 0d0+scale*cos(ii*EISCOR_DBL_PI/(N+1d0))
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


  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))



end program example_d_symtrid_qr_1dlaplace
