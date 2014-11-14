#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ZEXROU (Zomplex EXample Known EigenValues)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine chooses N random eigenvalues on the unit 
! circle, constructs a unitary matrix with prescribed 
! eigenvalues, and computes the N eigenvalues by solving
! a corresponding unitary eigenvalue problem two different 
! ways.
!
! 1) Form upper Hessenberg permutation matrix and compute its 
!    eigenvalues using ZUHFQR
!
! 2) Construct the factorization directly and compute its 
!    eigenvalues using ZUFFQR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program ZEXKEV

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO
  real(8) :: WORK(5*N), Q(3*(N-1)), D(2*N)
  complex(8) :: H(N,N), Z, he, temp(2,2)
  integer :: ITS(N-1), ind, jj
  double precision :: rm,ro,rp,pi = 3.141592653589793239d0
  double precision :: diff, her
  double precision :: nrm,b1(3),b2(3),b3(3),tb(3),tol=1e-16

  ! real and imag part of eigenvalues
  double precision :: rev(N), iev(N)



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print banner
  print*,""
  print*,"ZEXKEV: Zomple EXample Known EigenValues"
  print*,""
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initialize Q and D to be an identity matrix
  Q = 0d0
  D = 0d0
  do ii=1,n-1
     Q(3*ii-2) = 1d0
     D(2*ii-1) = 1d0
  end do
  D(2*n-1) = 1d0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! random eigenvalues, on the unit circle
  do ii = 1,n,2
     call random_number(rm)
     rev(ii) = sin(2*pi*rm)
     iev(ii) = cos(2*pi*rm)
     ! project on the unit circle
     nrm = rev(ii)**2+iev(ii)**2
     if (abs(nrm-1)<tol) then
        nrm = sqrt(nrm)
        rev(ii) = rev(ii)/nrm
        iev(ii) = iev(ii)/nrm
     end if
     Q(3*ii-2) = rev(ii)
     Q(3*ii)   = iev(ii)
     rev(ii+1) = rev(ii)
     iev(ii+1) = -iev(ii)
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! inverse eigenvalue problem 
  do ii=n-3,1,-2
     ! generate random bulge
     call random_number(rp)     
     call random_number(ro)     
     call random_number(rm)     
     call ZARCG33(rm,ro,rp,b1(1),b1(2),b1(3),NRM,INFO)
     call random_number(rp)     
     call random_number(ro)     
     call random_number(rm)     
     call ZARCG33(rm,ro,rp,b2(1),b2(2),b2(3),NRM,INFO)
     ! fuse on the top of the current Q sequence
     tb(1) = b2(1)
     tb(2) = -b2(2)
     tb(3) = -b2(3)
     b3(1) = b1(1)
     b3(2) = -b1(2)
     b3(3) = -b1(3)
     ind = 3*(ii-1)
     call ZARGTO(tb,b3,Q((ind+1):(ind+3)),INFO)
     Q((ind+4):(ind+6)) = b3
     b3 = Q((ind+1):(ind+3))
     Q((ind+1):(ind+3)) = tb
     ! main chasing loop
     do jj=ii,(n-3)
        ! set indices
        ind = 3*(jj-1)
        ! through diag
        call ZARGTD('R',D(2*jj+1:2*jj+4),b1,INFO)
        call ZARGTO(Q((ind+4):(ind+6)),Q((ind+7):(ind+9)),b1,INFO)
        call ZARGTD('R',D(2*jj-1:2*jj+2),b2,INFO)
        call ZARGTO(Q((ind+1):(ind+3)),Q((ind+4):(ind+6)),b2,INFO)
        call ZARGTO(b3,b1,b2,INFO)
        ! update bulges
        tb = b2
        b2 = b3
        b3 = b1
        b1 = tb
     end do
     ind = 3*(n-3)
     ! fusion at bottom
     call ZUFFGR('B',N,ii,N-1,Q,D,b1,INFO)
     call ZARGTD('R',D(2*n-5:2*n-2),b2,INFO)
     call ZARGTO(Q((ind+1):(ind+3)),Q((ind+4):(ind+6)),b2,INFO)
     call ZUFFGR('B',N,ii,N-1,Q,D,b3,INFO)
     call ZUFFGR('B',N,ii,N-1,Q,D,b2,INFO)
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initialize H to be an upper Hessenberg permutation matrix
  H = 0d0
  do ii=1,N
     H(ii,ii)=1d0
  end do
  H(1,1) = cmplx(Q(1),Q(2),kind=8)
  H(1,2) = Q(3)
  H(2,1) = -Q(3)
  H(2,2) = cmplx(Q(1),-Q(2),kind=8)
  do ii=2,N-1
     temp(1,1) = cmplx(Q(3*ii-2),Q(3*ii-1),kind=8)
     temp(1,2) = Q(3*ii)
     temp(2,1) = -Q(3*ii)
     temp(2,2) = cmplx(Q(3*ii-2),-Q(3*ii-1),kind=8)
     do jj = 1,ii+1
        he = H(jj,ii)*temp(1,1) + H(jj,ii+1)*temp(2,1)
        H(jj,ii+1) = H(jj,ii)*temp(1,2) + H(jj,ii+1)*temp(2,2)
        H(jj,ii) = he
     end do
  end do
  do ii=1,N
     H(:,ii) = H(:,ii)*cmplx(D(2*ii-1),D(2*ii),kind=8)
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print*,"Eigenvalues:"
  do ii=1,N
     print*, rev(ii), iev(ii)
  end do
  print*,""
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! call dohfqr
  call ZUHFQR('N',N,H,Z,ITS,WORK,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"ZUHFQR failed."
    print*,"INFO:",INFO
  end if
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print diag of H
  print*,"Eigenvalues computed using ZUHFQR: (real part, imag part, distance to closest exact eigenvalue)"
  do ii=1,N
     diff = sqrt(abs(dble(H(ii,ii))-rev(1))**2+abs(aimag(H(ii,ii))-iev(1))**2)
     do jj=2,N
        her = sqrt(abs(dble(H(ii,ii))-rev(jj))**2+abs(aimag(H(ii,ii))-iev(jj))**2)
        if (her < diff) then
           diff = her
        end if
     end do
     print*,dble(H(ii,ii)),aimag(H(ii,ii)),diff
  end do
  print*,""
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! call doffqr
  call ZUFFQR('N',N,Q,D,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"ZUFFQR failed."
    print*,"INFO:",INFO
  end if
  
  ! print D
  print*,"Eigenvalues computed using ZUFFQR: (real part, imag part, distance to closest exact eigenvalue)"
  do ii=1,N
     diff = sqrt(abs(D(2*ii-1)-rev(1))**2+abs(D(2*ii)-iev(1))**2)
     do jj=2,N
        her = sqrt(abs(D(2*ii-1)-rev(jj))**2+abs(D(2*ii)-iev(jj))**2)
        if (her < diff) then
           diff = her
        end if
     end do
     print*,D(2*ii-1),D(2*ii), diff
  end do
  print*,""

    
end program ZEXKEV
