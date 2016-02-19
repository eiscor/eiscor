#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_symtrid_qr_1dlaplace
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the eigenvalues of the tridiagonal Toeplitz
! matrix [ -0.5 0 -0.5 ]. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_symtrid_qr_1dlaplace

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, jj, kk, INFO, IWORK(3+5*N)
  real(8) :: WORK(5*N), D(N), E(N), Z, t, t2
  integer :: ITS(N-1)
  
  ! print banner
  print*,""
  print*,"example_d_symtrid_qr_1dlaplace:"
  print*,""

  ! initialize T to be a tridiagonal matrix of the form
  !   0  -0.5
  ! -0.5   0 -0.5
  !      -0.5  0  ...
  D = 0d0
  E = -5d-1

  ! call d_symtrid_qr
  call d_symtrid_qr(.FALSE.,.FALSE.,.FALSE.,N,D,E,WORK,1,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"d_symtrid_qr failed."
    print*,"INFO:",INFO
  end if
   
  ! print D
  print*,"Roots computed using d_symtrid_qr:"
  print*,"         ii computed roots            closest exact root          difference" 
  do ii=1,N
     E(ii) = 0d0+cos(ii*EISCOR_DBL_PI/(N+1d0))
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
     
     print*, ii, D(ii), E(INFO), t

     if (t.GT.t2) then
        t2 = t
        kk = INFO
     end if
     Z = Z + t**2
  end do

  print*,""
  print*, "forward error (2-norm of the vector) ", dsqrt(Z)
  print*, "largest forward error",t2
  print*, "largest forward error for eigenvalue number", kk 
  print*, "eigenvalue", E(kk)
  print*,""

    

end program example_d_symtrid_qr_1dlaplace
