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
  integer :: N 
  real(8) :: res, forw, err
  integer :: ii, jj, kk, ll, INFO
  real(8) :: normofp, lambda
  real(8), allocatable :: RESIDUALS(:)
  complex(8), allocatable :: COEFFS(:), ROOTS(:)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print banner
  print*,""
  print*,"example_z_poly_roots"
  print*,""

  
  N = 10
  
  allocate(RESIDUALS(N),COEFFS(N+1),ROOTS(N))
  ! compute random coeffs
  
  call z_1Darray_random_normal(N+1,COEFFS)  
        
  ! call roots
  call z_poly_roots(N,COEFFS,ROOTS,RESIDUALS,INFO)
  
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
  
  deallocate(RESIDUALS,COEFFS,ROOTS)
  
  ! end check 1)


  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  print*, "Runtime of this example:  ", dble(c_stop-c_start)/dble(c_rate)
  print*,""

     
end program example_z_poly_roots
