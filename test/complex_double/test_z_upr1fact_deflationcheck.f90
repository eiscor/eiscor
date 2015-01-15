#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test_z_upr1fact_deflationcheck
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This program tests the subroutine z_upr1fact_deflationcheck. 
! The following tests are run:
!
! 1) upper hessenberg
! 2) inverse hessenberg
! 3) CMV
! 4) inverse CMV
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program test_z_upr1fact_deflationcheck

  implicit none
  
  ! compute variables
  integer, parameter :: N = 6
  integer :: ii, ind, INFO, ITCNT, STR, STP, ZERO, ITS(N-1)
  logical :: P(N-2)
  real(8) :: Q(3*(N-1)), D(2*(N+1))
  complex(8) :: H1(N+1,N+1), H2(N+1,N+1)
  
  ! tolerance
  real(8), parameter :: tol = 10d0*epsilon(1d0)
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! initialize INFO
  INFO = 0
  
  ! initialize ITS
  ITS = 0
  
  ! print banner
  call u_test_banner(__FILE__)
  
  ! check 1)
    ! initialize P
    P = .FALSE.
    
    ! initialize Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii-2) = 1d0/sqrt(2d0)
      Q(3*ii) = 1d0/sqrt(2d0)
    end do
    ind = N/2
    Q(3*ind-2) = 0d0
    Q(3*ind-1) = 1d0
    Q(3*ind) = 0d0
    
    ! initialize D
    D = 0d0
    do ii=1,(N+1)
      D(2*ii-1) = 1d0
    end do
        
    ! initialize H1
    call z_upr1fact_form_hess_matrix(N,P,Q,D,H1)
    
    ! set ITCNT
    ITCNT = 10
    
    ! set STR, STP, ZERO
    STR = 1
    STP = N-1
    ZERO = 0
    
    ! call z_upr1fact_deflationcheck
    call z_upr1fact_deflationcheck(N,STR,STP,ZERO,P,Q,D,ITCNT,ITS,INFO)
    
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
    
    ! check ZERO
    if (ZERO.NE.ind) then
      call u_test_failed(__LINE__)
    end if   
    
    ! check STR
    if (STR.NE.(ind+1)) then
      call u_test_failed(__LINE__)
    end if  
    
    ! check ITCNT
    if (ITCNT.NE.0) then
      call u_test_failed(__LINE__)
    end if  

    ! check ITS
    if (ITS(ind).NE.10) then
      call u_test_failed(__LINE__)
    end if

  ! check 2)
    ! initialize P
    P = .TRUE.
    
    ! initialize Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii-2) = 1d0/sqrt(2d0)
      Q(3*ii) = 1d0/sqrt(2d0)
    end do
    ind = N/2
    Q(3*ind-2) = 0d0
    Q(3*ind-1) = 1d0
    Q(3*ind) = 0d0
    
    ! initialize D
    D = 0d0
    do ii=1,(N+1)
      D(2*ii-1) = 1d0
    end do
        
    ! initialize H1
    call z_upr1fact_form_hess_matrix(N,P,Q,D,H1)
    
    ! set ITCNT
    ITCNT = 10
    
    ! set STR, STP, ZERO
    STR = 1
    STP = N-1
    ZERO = 0
    
    ! call z_upr1fact_deflationcheck
    call z_upr1fact_deflationcheck(N,STR,STP,ZERO,P,Q,D,ITCNT,ITS,INFO)
    
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
    
    ! check ZERO
    if (ZERO.NE.ind) then
      call u_test_failed(__LINE__)
    end if   
    
    ! check STR
    if (STR.NE.(ind+1)) then
      call u_test_failed(__LINE__)
    end if  
    
    ! check ITCNT
    if (ITCNT.NE.0) then
      call u_test_failed(__LINE__)
    end if  

    ! check ITS
    if (ITS(ind).NE.10) then
      call u_test_failed(__LINE__)
 
    end if
    
  ! check 3)
    ! initialize P
    P = .FALSE.
    ii = 2
    do while (ii <= (N-2))
      P(ii) = .TRUE.
      ii = ii + 2
    end do
    
    ! initialize Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii-2) = 1d0/sqrt(2d0)
      Q(3*ii) = 1d0/sqrt(2d0)
    end do
    ind = N/2
    Q(3*ind-2) = 0d0
    Q(3*ind-1) = 1d0
    Q(3*ind) = 0d0
    
    ! initialize D
    D = 0d0
    do ii=1,(N+1)
      D(2*ii-1) = 1d0
    end do
        
    ! initialize H1
    call z_upr1fact_form_hess_matrix(N,P,Q,D,H1)
    
    ! set ITCNT
    ITCNT = 10
    
    ! set STR, STP, ZERO
    STR = 1
    STP = N-1
    ZERO = 0
    
    ! call z_upr1fact_deflationcheck
    call z_upr1fact_deflationcheck(N,STR,STP,ZERO,P,Q,D,ITCNT,ITS,INFO)
    
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
    
    ! check ZERO
    if (ZERO.NE.ind) then
      call u_test_failed(__LINE__)
    end if   
    
    ! check STR
    if (STR.NE.(ind+1)) then
      call u_test_failed(__LINE__)
    end if  
    
    ! check ITCNT
    if (ITCNT.NE.0) then
      call u_test_failed(__LINE__)
    end if  

    ! check ITS
    if (ITS(ind).NE.10) then
      call u_test_failed(__LINE__)
    end if
    
  ! check 4)
    ! initialize P
    P = .TRUE.
    ii = 2
    do while (ii <= (N-2))
      P(ii) = .FALSE.
      ii = ii + 2
    end do
    
    ! initialize Q
    Q = 0d0
    do ii=1,(N-1)
      Q(3*ii-2) = 1d0/sqrt(2d0)
      Q(3*ii) = 1d0/sqrt(2d0)
    end do
    ind = N/2
    Q(3*ind-2) = 0d0
    Q(3*ind-1) = 1d0
    Q(3*ind) = 0d0
    
    ! initialize D
    D = 0d0
    do ii=1,(N+1)
      D(2*ii-1) = 1d0
    end do
        
    ! initialize H1
    call z_upr1fact_form_hess_matrix(N,P,Q,D,H1)
    
    ! set ITCNT
    ITCNT = 10
    
    ! set STR, STP, ZERO
    STR = 1
    STP = N-1
    ZERO = 0
    
    ! call z_upr1fact_deflationcheck
    call z_upr1fact_deflationcheck(N,STR,STP,ZERO,P,Q,D,ITCNT,ITS,INFO)
    
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
    
    ! check ZERO
    if (ZERO.NE.ind) then
      call u_test_failed(__LINE__)
    end if   
    
    ! check STR
    if (STR.NE.(ind+1)) then
      call u_test_failed(__LINE__)
    end if  
    
    ! check ITCNT
    if (ITCNT.NE.0) then
      call u_test_failed(__LINE__)
    end if  

    ! check ITS
    if (ITS(ind).NE.10) then
      call u_test_failed(__LINE__)
    end if
 
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  call u_test_passed(dble(c_stop-c_start)/dble(c_rate))
  
end program test_z_upr1fact_deflationcheck

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
