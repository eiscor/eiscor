#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_mergebulge
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_mergebulge. 
! The following tests are run:
!
! 1) Q = [0, -1; 1, 0], P = FALSE, D = I, K = 1, JOB = L and 
!    G = [1+1i,1]/sqrt(3)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_mergebulge

  implicit none
  
  ! compute variables
  integer, parameter :: N = 3
  real(8), parameter :: tol = 1d1*epsilon(1d0)
  integer :: ii
  logical :: P(N-2)
  real(8) :: Q(3*(N-1)),D(2*(N+1)), G(3)
  complex(8) :: temp(2,2), H1(N+1,N+1), H2(N+1,N+1)
  integer :: INFO
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! check 1)
    ! set INFO
    INFO = 0
  
    ! set Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii) = 1d0
    end do
    
    ! set D
    D = 0d0
    do ii=1,(N+1)
      D(2*ii-1) = 1d0
    end do
    
    ! set G
    G(1) = 1d0/sqrt(3d0)
    G(2) = 1d0/sqrt(3d0)
    G(3) = 1d0/sqrt(3d0)
    
    ! set P
    P = .FALSE.
    
    ! initialize H1
    call z_upr1fact_form_hess_matrix(N,P,Q,D,H1)
    
    temp(1,1) = cmplx(G(1),G(2),kind=8)
    temp(2,1) = cmplx(G(3),0d0,kind=8)
    temp(1,2) = -temp(2,1)
    temp(2,2) = conjg(temp(1,1))
    
    H1(1:2,:) = matmul(temp,H1(1:2,:))
    
    ! call merge bulge
    call z_upr1fact_mergebulge('L',N,1,N-1,1,P,Q,D,G,INFO)
  
    ! check INFO
    if (INFO.NE.0) then
      call u_test_failed(__LINE__)
    end if
    
    ! initialize H2
    call z_upr1fact_form_hess_matrix(N,P,Q,D,H2)
    
    ! check difference
    if (maxval(abs(H1-H2)) > tol) then
      call u_test_failed(__LINE__)
    end if 
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))

end program test_z_upr1fact_mergebulge

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This subroutine computes the (N+1)x(N+1) extended hessenberg matrix defined
! by Given's rotations stored in Q whose order is described by P and the
! diagonal matrix D. The output is stored in H.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_form_hess_matrix(N,P,Q,D,H)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  logical, intent(in) :: P(N-2)
  real(8), intent(in) :: Q(3*(N-1)), D(2*(N+1))
  complex(8), intent(inout) :: H(N+1,N+1)
  
  ! compute variables
  integer :: ii
  complex :: temp(2,2)
  
  ! initialize H
  H = cmplx(0d0,0d0,kind=8)
  do ii=1,(N+1)
    H(ii,ii) = cmplx(1d0,0d0,kind=8)
  end do
  H(1,1) = cmplx(Q(1),Q(2),kind=8)
  H(2,1) = cmplx(Q(3),0d0,kind=8)
  H(1,2) = -H(2,1)
  H(2,2) = conjg(H(1,1))
  
  ! apply Q
  do ii=1,(N-2)
    temp(1,1) = cmplx(Q(3*ii+1),Q(3*ii+2),kind=8)
    temp(2,1) = cmplx(Q(3*ii+3),0d0,kind=8)
    temp(1,2) = -temp(2,1)
    temp(2,2) = conjg(temp(1,1))
    
    if (P(ii).EQV..FALSE.) then
      H(:,(ii+1):(ii+2)) = matmul(H(:,(ii+1):(ii+2)),temp)  
    else
      H((ii+1):(ii+2),:) = matmul(temp,H((ii+1):(ii+2),:))  
    end if
  
  end do
  
  ! apply D
  do ii=1,(N+1)
    H(:,ii) = H(:,ii)*cmplx(D(2*ii-1),D(2*ii),kind=8)
  end do
  
end subroutine z_upr1fact_form_hess_matrix
