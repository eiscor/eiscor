#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_d_symtrid_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the eigenvalues of the tridiagonal matrix
! [-0.5 0 -0.5] testing the forward error and of a random normally 
! distributed matrix testing the backward error.  
!
! check 1) [-0.5 0 -0.5]
! check 2) random tridiagonal matrix
! check 3) random tridiagonal matrix with first row/column zero
! check 4) random tridiagonal matrix with last row/column zero
! check 5) random tridiagonal matrix with row 5 and column 5 zero
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_d_symtrid_qr  

  implicit none
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! parameters
  integer, parameter :: N = 16
  real(8), parameter :: scale = 1d0
  real(8), parameter :: shift = 0d0
  logical, parameter :: sca = .TRUE.
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! compute variables
  integer :: M
  integer :: ii, jj, INFO
  real(8) :: WORK(5*N), D(N), E(N-1), eig(N), t, t1, t2, t3, nrm
  real(8) :: Ds(N), Es(N-1), pi = EISCOR_DBL_PI
  complex(8) :: Z(N,N), v(N)
  integer :: ITS(N-1)
  logical :: backward
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  call u_test_banner(__FILE__)

  ! check 1) [-0.5 0 -0.5]
  ! initialize T to be a tridiagonal matrix of the form
  !  2 -1
  ! -1  2 -1
  !     -1 2 ...
  Ds = 0d0
  Es = -5d-1

  do ii=1,N
     Ds(ii) = (Ds(ii)-shift)*scale
  end do
  do ii=1,N-1
     Es(ii) = Es(ii)*scale
  end do

  D = Ds
  E = Es

  call d_symtrid_qr(.TRUE.,.TRUE.,sca,N,D,E,WORK,N,Z,ITS,INFO)

  ! check INFO
  if (INFO.NE.0) then
     call u_test_failed(__LINE__)
  end if

  nrm = 0d0
  do ii=1,N
     if (abs(D(ii))>nrm) then
        nrm = abs(D(ii))
     end if
  end do
  
  ! computing forward error
  ! exact eigenvalues
  do ii=1,N
     eig(ii) = (-scale*shift)+scale*cos(ii*EISCOR_DBL_PI/(N+1d0))
  end do
  
  t1 = 0d0
  t2 = 0d0
  
  do ii=1,N
     t = abs(D(ii)-eig(1))
     INFO = 1
     do jj=2,N
        if (abs(D(ii)-eig(jj))<t) then
           t = abs(D(ii)-eig(jj))
           INFO = jj
        end if
     end do
              
     if (t.GT.t2) then
        t2 = t
     end if
     t1 = t1 + t**2
  end do
  t1 = t1/nrm

  if ((dsqrt(t1)>1d-13).OR.(t1.NE.t1)) then
     ! forward error test failed
     call u_test_failed(__LINE__)
  end if

  ! computing backward error
  t2 = 0d0
  do ii=1,N
     ! jj  1
     v(1) = (Ds(1)-D(ii))*Z(1+N*(ii-1),1) + Es(1)*Z(2+N*(ii-1),1)
     do jj=2,N-1
        v(jj) = (Ds(jj)-D(ii))*Z(jj+N*(ii-1),1) + &
             & Es(jj)*Z(jj+1+N*(ii-1),1) + &
             & Es(jj-1)*Z(jj-1+N*(ii-1),1)
     end do
     v(N) = (Ds(N)-D(ii))*Z(N+N*(ii-1),1) + &
          & Es(N-1)*Z(N-1+N*(ii-1),1)

     t3 = 0d0
     do jj=1,N
        t3 = t3 + abs(v(jj))**2
     end do
     t3 = sqrt(t3)
     
     if (t3.GT.t2) then
        t2 = t3
     end if
  end do

  t2 = t2/nrm

  if ((t2>1d-14).OR.(t2.NE.t2)) then
     ! backward error test failed
     call u_test_failed(__LINE__)
  end if


  ! check 2) random tridiagonal matrix
  !call u_randomseed_initialize(INFO)
  call u_fixedseed_initialize(INFO)
  do ii=1,N-1
     call random_number(t)
     Ds(ii) = t
     call random_number(t)
     Es(ii) = t
  end do
  call random_number(t)
  Ds(N) = t

  D = Ds
  E = Es

  call d_symtrid_qr(.TRUE.,.TRUE.,sca,N,D,E,WORK,N,Z,ITS,INFO)

  ! computing backward error
  t2 = 0d0
  do ii=1,N
     ! jj  1
     v(1) = (Ds(1)-D(ii))*Z(1+N*(ii-1),1) + Es(1)*Z(2+N*(ii-1),1)
     do jj=2,N-1
        v(jj) = (Ds(jj)-D(ii))*Z(jj+N*(ii-1),1) + &
             & Es(jj)*Z(jj+1+N*(ii-1),1) + &
             & Es(jj-1)*Z(jj-1+N*(ii-1),1)
     end do
     v(N) = (Ds(N)-D(ii))*Z(N+N*(ii-1),1) + &
          & Es(N-1)*Z(N-1+N*(ii-1),1)

     t3 = 0d0
     do jj=1,N
        t3 = t3 + abs(v(jj))**2
     end do
     t3 = sqrt(t3)

     if (t3.GT.t2) then
        t2 = t3
     end if
  end do

  t2 = t2/nrm

  if ((t2>1d-14).OR.(t2.NE.t2)) then
     ! backward error test failed
     call u_test_failed(__LINE__)
  end if
  
  ! check 3) random tridiagonal matrix with first row/column zero
  call u_fixedseed_initialize(INFO)
  do ii=1,N-1
     call random_number(t)
     Ds(ii) = t
     call random_number(t)
     Es(ii) = t
  end do
  call random_number(t)
  Ds(N) = t

  Ds(1) = 0d0
  Es(1) = 0d0
  
  D = Ds
  E = Es

  call d_symtrid_qr(.TRUE.,.TRUE.,sca,N,D,E,WORK,N,Z,ITS,INFO)

  ! computing backward error
  t2 = 0d0
  do ii=1,N
     ! jj  1
     v(1) = (Ds(1)-D(ii))*Z(1+N*(ii-1),1) + Es(1)*Z(2+N*(ii-1),1)
     do jj=2,N-1
        v(jj) = (Ds(jj)-D(ii))*Z(jj+N*(ii-1),1) + &
             & Es(jj)*Z(jj+1+N*(ii-1),1) + &
             & Es(jj-1)*Z(jj-1+N*(ii-1),1)
     end do
     v(N) = (Ds(N)-D(ii))*Z(N+N*(ii-1),1) + &
          & Es(N-1)*Z(N-1+N*(ii-1),1)

     t3 = 0d0
     do jj=1,N
        t3 = t3 + abs(v(jj))**2
     end do
     t3 = sqrt(t3)

     if (t3.GT.t2) then
        t2 = t3
     end if
  end do

  t2 = t2/nrm

  if ((t2>1d-14).OR.(t2.NE.t2)) then
     ! backward error test failed
     call u_test_failed(__LINE__)
  end if

  ! check 4) random tridiagonal matrix with last row/column zero
  call u_fixedseed_initialize(INFO)
  do ii=1,N-1
     call random_number(t)
     Ds(ii) = t
     call random_number(t)
     Es(ii) = t
  end do
  call random_number(t)
  Ds(N) = t

  Ds(N) = 0d0
  Es(N-1) = 0d0
  
  D = Ds
  E = Es

  call d_symtrid_qr(.TRUE.,.TRUE.,sca,N,D,E,WORK,N,Z,ITS,INFO)

  ! computing backward error
  t2 = 0d0
  do ii=1,N
     ! jj  1
     v(1) = (Ds(1)-D(ii))*Z(1+N*(ii-1),1) + Es(1)*Z(2+N*(ii-1),1)
     do jj=2,N-1
        v(jj) = (Ds(jj)-D(ii))*Z(jj+N*(ii-1),1) + &
             & Es(jj)*Z(jj+1+N*(ii-1),1) + &
             & Es(jj-1)*Z(jj-1+N*(ii-1),1)
     end do
     v(N) = (Ds(N)-D(ii))*Z(N+N*(ii-1),1) + &
          & Es(N-1)*Z(N-1+N*(ii-1),1)

     t3 = 0d0
     do jj=1,N
        t3 = t3 + abs(v(jj))**2
     end do
     t3 = sqrt(t3)

     if (t3.GT.t2) then
        t2 = t3
     end if
  end do

  t2 = t2/nrm

  if ((t2>1d-14).OR.(t2.NE.t2)) then
     ! backward error test failed
     call u_test_failed(__LINE__)
  end if

  ! check 5) random tridiagonal matrix with row 5 and column 5 zero
  call u_fixedseed_initialize(INFO)
  do ii=1,N-1
     call random_number(t)
     Ds(ii) = t
     call random_number(t)
     Es(ii) = t
  end do
  call random_number(t)
  Ds(N) = t

  Ds(5) = 0d0
  Es(4) = 0d0
  Es(5) = 0d0
  
  D = Ds
  E = Es

  call d_symtrid_qr(.TRUE.,.TRUE.,sca,N,D,E,WORK,N,Z,ITS,INFO)

  ! computing backward error
  t2 = 0d0
  do ii=1,N
     ! jj  1
     v(1) = (Ds(1)-D(ii))*Z(1+N*(ii-1),1) + Es(1)*Z(2+N*(ii-1),1)
     do jj=2,N-1
        v(jj) = (Ds(jj)-D(ii))*Z(jj+N*(ii-1),1) + &
             & Es(jj)*Z(jj+1+N*(ii-1),1) + &
             & Es(jj-1)*Z(jj-1+N*(ii-1),1)
     end do
     v(N) = (Ds(N)-D(ii))*Z(N+N*(ii-1),1) + &
          & Es(N-1)*Z(N-1+N*(ii-1),1)

     t3 = 0d0
     do jj=1,N
        t3 = t3 + abs(v(jj))**2
     end do
     t3 = sqrt(t3)

     if (t3.GT.t2) then
        t2 = t3
     end if
  end do

  t2 = t2/nrm

  if ((t2>1d-14).OR.(t2.NE.t2)) then
     ! backward error test failed
     call u_test_failed(__LINE__)
  end if

  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))  

end program test_d_symtrid_qr
