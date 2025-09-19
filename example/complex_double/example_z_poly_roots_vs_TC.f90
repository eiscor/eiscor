#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_z_poly_roots
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program runs a simple examples with random coefficients.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_z_poly_roots

  implicit none
  
  ! compute variables
  integer :: N, Q
  real(8) :: res, forw, err
  integer :: ii, jj, kk, ll, INFO
  real(8) :: normofp, lambda
  real(8), allocatable :: RESIDUALS(:)
  complex(8), allocatable :: COEFFS(:), ROOTS(:)
  real(8) R1, R2, R3
  real(8) :: t1, t2, L
    
  ! timing variables
  integer:: c_start, c_stop, c_rate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print banner
  print*,""
  print*,"example_z_poly_roots_vs_TC"
  print*,""

  
  N = 3
  Q = 1000000
  
  allocate(RESIDUALS(N),COEFFS(N+1),ROOTS(N))
  ! compute random coeffs
  
  call z_1Darray_random_normal(N+1,COEFFS)  

  do ii=1,N+1
     COEFFS(ii)  = cmplx(dble(COEFFS(ii)),0d0)
  end do
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! call roots
  do ii=1,Q
     call z_poly_roots(N,COEFFS,ROOTS,RESIDUALS,INFO)
  end do
  call system_clock(count=c_stop)
  t1 = dble(c_stop-c_start)/dble(c_rate)/Q



  ! check INFO
  if (INFO.NE.0) then
    !call u_test_failed(__LINE__)
    print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print*, "!   INFO not 0                                   !"
    print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print*, INFO
  end if

  print*, "Roots                                                  Residual"
  do ii=1,N
    print*, ROOTS(ii), RESIDUALS(ii)
  end do
  print*, ""
  
  ! compute normofp
  normofp = 0d0
  do ii=1,N
    normofp = normofp + abs(COEFFS(ii))**2
  end do
  normofp = dsqrt(normofp)
  
  print*, "Norm of the coefficients: ", normofp
  
  ! maximum residuals
  res = 0d0
  do ii=1,N
    if (RESIDUALS(ii) >= res) then
      res = RESIDUALS(ii)
    end if
  end do
     
  print*, "Maximum redidual:         ", res

  print*, "relative residual:        ", res/normofp 


  
  call system_clock(count=c_start)
  
  ! call roots
  do ii=1,Q
     call TC(COEFFS(1),COEFFS(2),COEFFS(3),COEFFS(4),R1,R2,R3,L)
  end do
  call system_clock(count=c_stop)
  t2 = dble(c_stop-c_start)/dble(c_rate)/Q


  if (L.EQ.0) then
     ROOTS(1) = cmplx(R1,0d0,kind=8)
     ROOTS(2) = cmplx(R2,0d0,kind=8)
     ROOTS(3) = cmplx(R3,0d0,kind=8)
  else
     ROOTS(1) = cmplx(R1,0d0,kind=8)
     ROOTS(2) = cmplx(R2,R3,kind=8)
     ROOTS(3) = cmplx(R2,-R3,kind=8)
  end if
     
  print*, "Roots"
  do ii=1,N
    print*, ROOTS(ii)
  end do
  print*, ""

  print*, "eiscor", t1, "TC", t2, t1/t2

  
  deallocate(RESIDUALS,COEFFS,ROOTS)
  
  ! end check 1)


  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  print*, "Runtime of this example:  ", dble(c_stop-c_start)/dble(c_rate)
  print*,""

     
end program example_z_poly_roots
