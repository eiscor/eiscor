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
program example_d_symtrid_qr_1dl_be

  implicit none
  
  ! compute variables
  !integer, parameter :: N = 4
  !integer, parameter :: N = 16
  integer, parameter :: N = 200
  !integer, parameter :: N = 6184
  !integer, parameter :: N = 4092
  !integer, parameter :: N = 2048
  !integer, parameter :: N = 1024
  !integer, parameter :: N = 512
  !integer, parameter :: N = 256
  real(8), parameter :: scale = 1d-6
  integer :: ii, jj, kk, INFO, IWORK(5*(N+1))
  real(8) :: WORK(14*N+N*N), D(N), E(N),  t, t1, t2
  real(8) :: Ds(N), Es(N), Hr(N,N), Zr(N,N), ZrH(N,N)
  complex(8) :: Z(N,N), H(N,N), ZH(N,N), c1
  integer :: ITS(N-1)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! BLAS DNRM2
  double precision :: dnrm2, dznrm2

  
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
  call d_symtrid_qr(.TRUE.,.TRUE.,N,D,E,WORK,N,Z,ITS,INFO)
  
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

  t1 = 0d0
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
     
     if (N<129) then
        print*, D(ii), E(INFO), t
     end if
     if (t.GT.t2) then
        t2 = t
        kk = INFO
     end if
     t1 = t1 + t**2
  end do
  print*, "forward error (2-norm of the vector) ", dsqrt(t1)
  if (dsqrt(t1)>1d-3) then
     call u_test_failed(__LINE__)
  end if
  print*,""
  print*, "largest forward error for eigenvalue ", kk, " error ",t2," eigenvalue ", E(kk)
  print*,""

  ! compute residual matrix
  H = cmplx(0d0,0d0,kind=8)
  do ii=1,N
     H(ii,ii) = cmplx(D(ii),0d0,kind=8)
  end do
  do ii=1,N
     do jj=1,N
        ZH(ii,jj) = conjg(Z(jj,ii))
     end do
  end do
!!$  print*, ""
!!$  print*, "D"
!!$  do ii=1,N
!!$     do jj=1,N
!!$        print*, ii,jj, H(ii,jj)
!!$     end do
!!$  end do

!!$  print*, ""
!!$  print*, "Z"
!!$  do ii=1,N
!!$     do jj=1,N
!!$        print*, ii,jj, Z(ii,jj)
!!$     end do
!!$  end do
  H(1:N,1:N) = matmul(H(1:N,1:N),ZH(1:N,1:N))
  H(1:N,1:N) = matmul(Z(1:N,1:N),H(1:N,1:N))

!!$  print*, ""
!!$  print*, "ZH D Z"
!!$  do ii=1,N
!!$     do jj=1,N
!!$        print*, ii,jj, H(ii,jj)
!!$     end do
!!$  end do
  do ii=1,N
     H(ii,ii) = H(ii,ii) - Ds(ii)
  end do
!!$  print*, ""
!!$  print*, "H-Es"
  do ii=1,N-1
     H(ii+1,ii) = H(ii+1,ii) - Es(ii)
     H(ii,ii+1) = H(ii,ii+1) - Es(ii)
  end do

  t1 = dznrm2(N*N,H,1)
  print*, "backward error Frobenius norm ", t1
  if (t1>1d-3) then
     call u_test_failed(__LINE__)
  end if

!!$  do ii=1,N
!!$     do jj=1,N
!!$        print*, ii,jj, H(ii,jj)
!!$     end do
!!$     !print*, H(ii,:)
!!$  end do

  ! analysis of angles in the first eigenvector
!!$  do jj=1,N
!!$     c1 = Z(1,jj)
!!$     do ii=1,N
!!$        print*, ii,jj, atan(real(Z(ii,jj))/aimag(Z(ii,jj))), Z(ii,jj)/c1
!!$     end do
!!$  end do
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
  call dsteqr ('I', N, D, E, Zr, N, WORK, INFO)
  ! stop timer
  call system_clock(count=c_stop) 

  ! check results
  do ii=1,N
     E(ii) = 0d0+scale*cos(ii*EISCOR_DBL_PI/(N+1d0))
  end do
  t1 = 0d0
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
     t1 = t1 + t**2
  end do
  print*, "forward error (2-norm of the vector) ", dsqrt(t1)
  print*,""

  ! compute residual matrix
  Hr = 0d0
  do ii=1,N
     Hr(ii,ii) = D(ii)
  end do
  do ii=1,N
     do jj=1,N
        ZrH(ii,jj) = Zr(jj,ii)
     end do
  end do
!!$  print*, ""
!!$  print*, "D"
!!$  do ii=1,N
!!$     do jj=1,N
!!$        print*, ii,jj, H(ii,jj)
!!$     end do
!!$  end do

!!$  print*, ""
!!$  print*, "Z"
!!$  do ii=1,N
!!$     do jj=1,N
!!$        print*, ii,jj, Z(ii,jj)
!!$     end do
!!$  end do
  Hr(1:N,1:N) = matmul(Hr(1:N,1:N),ZrH(1:N,1:N))
  Hr(1:N,1:N) = matmul(Zr(1:N,1:N),Hr(1:N,1:N))

!!$  print*, ""
!!$  print*, "ZH D Z"
!!$  do ii=1,N
!!$     do jj=1,N
!!$        print*, ii,jj, Hr(ii,jj)
!!$     end do
!!$  end do
  do ii=1,N
     Hr(ii,ii) = Hr(ii,ii) - Ds(ii)
  end do
  do ii=1,N-1
     Hr(ii+1,ii) = Hr(ii+1,ii) - Es(ii)
     Hr(ii,ii+1) = Hr(ii,ii+1) - Es(ii)
  end do
!!$  print*, ""
!!$  print*, "H-Es"
!!$  do ii=1,N
!!$     do jj=1,N
!!$        print*, ii,jj, Hr(ii,jj)
!!$     end do
!!$  end do

  t1 = dnrm2(N*N,Hr,1)
  print*, "backward error Frobenius norm ", t1

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

  ! run DSTEVD
  call dstevd ('V', N, D, E, Z, N, WORK, 14*N+N*N, IWORK, 5+5*N, INFO) 
  ! stop timer
  call system_clock(count=c_stop) 

  ! check results
  do ii=1,N
     E(ii) = 0d0+scale*cos(ii*EISCOR_DBL_PI/(N+1d0))
  end do
  t1 = 0d0
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
     t1 = t1 + t**2
  end do
  print*, "forward error (2-norm of the vector) ", dsqrt(t1)
  print*,""

  ! compute residual matrix
  Hr = 0d0
  do ii=1,N
     Hr(ii,ii) = D(ii)
  end do
  do ii=1,N
     do jj=1,N
        ZrH(ii,jj) = Zr(jj,ii)
     end do
  end do
!!$  print*, ""
!!$  print*, "D"
!!$  do ii=1,N
!!$     do jj=1,N
!!$        print*, ii,jj, H(ii,jj)
!!$     end do
!!$  end do

!!$  print*, ""
!!$  print*, "Z"
!!$  do ii=1,N
!!$     do jj=1,N
!!$        print*, ii,jj, Z(ii,jj)
!!$     end do
!!$  end do
  Hr(1:N,1:N) = matmul(Hr(1:N,1:N),ZrH(1:N,1:N))
  Hr(1:N,1:N) = matmul(Zr(1:N,1:N),Hr(1:N,1:N))

!!$  print*, ""
!!$  print*, "ZH D Z"
!!$  do ii=1,N
!!$     do jj=1,N
!!$        print*, ii,jj, Hr(ii,jj)
!!$     end do
!!$  end do
  do ii=1,N
     Hr(ii,ii) = Hr(ii,ii) - Ds(ii)
  end do
  do ii=1,N-1
     Hr(ii+1,ii) = Hr(ii+1,ii) - Es(ii)
     Hr(ii,ii+1) = Hr(ii,ii+1) - Es(ii)
  end do
!!$  print*, ""
!!$  print*, "H-Es"
!!$  do ii=1,N
!!$     do jj=1,N
!!$        print*, ii,jj, Hr(ii,jj)
!!$     end do
!!$  end do

  t1 = dnrm2(N*N,Hr,1)
  print*, "backward error Frobenius norm ", t1

  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))


end program example_d_symtrid_qr_1dl_be
