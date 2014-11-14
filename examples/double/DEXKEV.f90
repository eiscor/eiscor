#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DEXROU (Double EXample Known EigenValues)
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
!    eigenvalues using DOHFQR
!
! 2) Construct the factorization directly and compute its 
!    eigenvalues using DOFFQR
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program DEXKEV

  implicit none
  
  ! compute variables
  integer, parameter :: N = 16
  integer :: ii, INFO, cpair
  real(8) :: WORK(4*N), Q(2*(N-1)), D(2*N)
  real(8) :: H(N,N), Z
  integer :: ITS(N-1), ind, jj
  double precision :: rm,rp,pi = 3.141592653589793239d0
  double precision :: temp(2,2), he, diff
  double precision :: nrm,b1(2),b2(2),b3(2),tb(2),tol

  ! real and imag part of eigenvalues
  double precision :: rev(N), iev(N) 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print banner
  print*,""
  print*,"DEXKEV: Double EXample Known EigenValues"
  print*,""
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initialize Q and D to be an identity matrix
  Q = 0d0
  D = 0d0
  do ii=1,n-1
     Q(2*ii-1) = 1d0
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
     Q(2*ii-1) = rev(ii)
     Q(2*ii)   = iev(ii)
     rev(ii+1) = rev(ii)
     iev(ii+1) = -iev(ii)
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! inverse eigenvalue problem 
  do ii=n-3,1,-2
     ! generate random bulge
     call random_number(rp)     
     call random_number(rm)     
     call DARCG22(sqrt(1-rp**2),rp,b1(1),b1(2),NRM,INFO)
     call DARCG22(sqrt(1-rm**2),rm,b2(1),b2(2),NRM,INFO)
     ! fuse on the top of the current Q sequence
     tb(1) = b2(1)
     tb(2) = -b2(2)
     b3(1) = b1(1)
     b3(2) = -b1(2)
     ind = 2*(ii-1)
     call DARGTO(tb,b3,Q((ind+1):(ind+2)),INFO)
     call DARFGR('R',b3,Q((ind+3):(ind+4)),INFO)
     b3 = Q((ind+1):(ind+2))
     Q((ind+1):(ind+2)) = tb
     ! main chasing loop
     do jj=ii,(n-3)
        ind = 2*(jj-1)
        ! set indices
        call DARGTO(Q((ind+3):(ind+4)),Q((ind+5):(ind+6)),b1,INFO)
        call DARGTO(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b2,INFO)
        call DARGTO(b3,b1,b2,INFO)
        ! update bulges
        tb = b2
        b2 = b3
        b3 = b1
        b1 = tb
     end do
     ind = 2*(n-3)
     ! fusion at bottom
     call DARFGR('L',Q((ind+3):(ind+4)),b1,INFO)
     call DARGTO(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b2,INFO)
     call DARFGR('L',b3,b2,INFO)
     call DARFGR('L',Q((ind+3):(ind+4)),b3,INFO)
  end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! initialize H to be an upper Hessenberg permutation matrix
  H = 0d0
  do ii=1,N
     H(ii,ii)=1d0
  end do
  H(1,1) = Q(1)
  H(1,2) = Q(2)
  H(2,1) = -Q(2)
  H(2,2) = Q(1)
  do ii=2,N-1
     temp(1,1) = Q(2*ii-1)
     temp(1,2) = Q(2*ii)
     temp(2,1) = -Q(2*ii)
     temp(2,2) = Q(2*ii-1)
     do jj = 1,ii+1
        he = H(jj,ii)*temp(1,1) + H(jj,ii+1)*temp(2,1)
        H(jj,ii+1) = H(jj,ii)*temp(1,2) + H(jj,ii+1)*temp(2,2)
        H(jj, ii) = he
     end do
  end do


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! call dohfqr
  call DOHFQR('N',N,H,Z,ITS,WORK,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"DOHFQR failed."
    print*,"INFO:",INFO
  end if
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print*,"Eigenvalues:"
  do ii=1,N
     print*, rev(ii), iev(ii)
  end do
  print*,""
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print diag of H
  cpair = 0
  print*,"Eigenvalues computed using DOHFQR: (real part, imag part, distance to closest exact eigenvalue)"
  do ii=1,N
    if (cpair.EQ.0) then
      if (H(ii,ii+1).EQ.0) then
         diff = abs(cmplx(H(ii,ii)-rev(1),iev(1),kind=8))
         do jj=2,N
            he = abs(cmplx(H(ii,ii)-rev(jj),iev(jj),kind=8))
            if (he < diff) then
               diff = he
            end if
         end do
         print*,H(ii,ii), 0d0, diff
      else
         diff = sqrt(abs(H(ii,ii)-rev(1))**2+abs(H(ii,ii+1)-iev(1))**2)
         do jj=2,N
            he = sqrt(abs(H(ii,ii)-rev(jj))**2+abs(H(ii,ii+1)-iev(jj))**2)
            if (he < diff) then
               diff = he
            end if
         end do
         print*,H(ii,ii),H(ii,ii+1),diff
         diff = sqrt(abs(H(ii+1,ii+1)-rev(1))**2+abs(H(ii+1,ii)-iev(1))**2)
         do jj=2,N
            he = sqrt(abs(H(ii+1,ii+1)-rev(jj))**2+abs(H(ii+1,ii)-iev(jj))**2)
            if (he < diff) then
               diff = he
            end if
         end do
         print*,H(ii+1,ii+1),H(ii+1,ii),diff
         cpair = 1
      end if
    else
      cpair = 0
    end if
  end do
  print*,""
  
  
  ! call doffqr
  call DOFFQR('N',N,Q,D,Z,ITS,INFO)
  
  ! check INFO
  if (INFO.NE.0) then
    print*,"DOFFQR failed."
    print*,"INFO:",INFO
  end if
  
  ! print D
  print*,"Eigenvalues computed using DOFFQR: (real part, imag part, distance to closest exact eigenvalue)"
  do ii=1,N
     diff = sqrt(abs(D(2*ii-1)-rev(1))**2+abs(D(2*ii)-iev(1))**2)
     do jj=2,N
        he = sqrt(abs(D(2*ii-1)-rev(jj))**2+abs(D(2*ii)-iev(jj))**2)
        if (he < diff) then
           diff = he
        end if
     end do
     print*,D(2*ii-1),D(2*ii), diff
  end do
  print*,""

    
end program DEXKEV
