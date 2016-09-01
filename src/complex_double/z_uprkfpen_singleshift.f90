#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkfpen_singleshift
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes a shift for a factored unitary plus rank one
! (uprkfpen) matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D1,D2           REAL(8) arrays of dimension (2*N*K)
!                    arrays of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1,C2,B1,B2     REAL(8) arrays of dimension (3*N*K)
!                    arrays of generators for unitary plus rank one
!                    upper-trinagular matrices
!
! OUTPUT VARIABLES:
!
!  SHFT            COMPLEX(8)
!                    shift
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkfpen_singleshift(N,K,STR,STP,P,Q,D1,C1,B1,D2,C2,B2,SHFT)
  
  implicit none
  
  ! input variables
  integer, intent(in) :: N, K, STR, STP
  logical, intent(in) :: P(N-2)
  real(8), intent(in) :: Q(3*(N-1)), D1(2*N*K), C1(3*N*K), B1(3*N*K)
  real(8), intent(in) :: D2(2*N*K), C2(3*N*K), B2(3*N*K)
  complex(8), intent(inout) :: SHFT
  
  ! compute variables
  complex(8) :: rho, R1(3,3), R2(3,3), H(2,2), HK(2,2)
  integer :: ind

  R1 = cmplx(0d0,0d0,kind=8)
  R2 = cmplx(0d0,0d0,kind=8)

  if (STR.EQ.STP) then

     ! extract R1
     call z_uprkutri_decompress(.FALSE.,N,K,STP,STP,D1,C1,B1,H)
     
     ! extract R2
     call z_uprkutri_decompress(.FALSE.,N,K,STP,STP,D2,C2,B2,HK)
    
     R1(2:3,2:3) = H
     R2(2:3,2:3) = HK
     
     ! apply second Q       
     ind = 3*STP-2
     
     ! update R1 ! Q(2)
     H(1,1) = cmplx(Q(ind),Q(ind+1),kind=8)
     H(2,1) = cmplx(Q(ind+2),0d0,kind=8)
     H(1,2) = cmplx(-Q(ind+2),0d0,kind=8)
     H(2,2) = cmplx(Q(ind),-Q(ind+1),kind=8)
     R1(2:3,:) = matmul(H,R1(2:3,:))

  else
     ! extract R1
     call z_uprkutri_decompress(.FALSE.,N,K,STP-1,STP,D1,C1,B1,R1)
     
     ! extract R2
     call z_uprkutri_decompress(.FALSE.,N,K,STP-1,STP,D2,C2,B2,R2)

     ! apply second Q
     if (P(STP-1)) then 

        ! inv hess        
        ind = 3*STP-2
        ! update R2 ! inv(Q(2))
        H(1,1) = cmplx(Q(ind),-Q(ind+1),kind=8)
        H(2,1) = cmplx(-Q(ind+2),0d0,kind=8)
        H(1,2) = cmplx(Q(ind+2),0d0,kind=8)
        H(2,2) = cmplx(Q(ind),Q(ind+1),kind=8)
        R2(2:3,:) = matmul(H,R2(2:3,:))
        
     else 

        ! hess        
        ind = 3*STP-2
        ! update R1 ! Q(2)
        H(1,1) = cmplx(Q(ind),Q(ind+1),kind=8)
        H(2,1) = cmplx(Q(ind+2),0d0,kind=8)
        H(1,2) = cmplx(-Q(ind+2),0d0,kind=8)
        H(2,2) = cmplx(Q(ind),-Q(ind+1),kind=8)
        R1(2:3,:) = matmul(H,R1(2:3,:))
        
     end if
     
     if (STP.EQ.2) then
        
        ind = 3*STP-5        
        ! update R1 ! Q(1)
        H(1,1) = cmplx(Q(ind),Q(ind+1),kind=8)
        H(2,1) = cmplx(Q(ind+2),0d0,kind=8)
        H(1,2) = cmplx(-Q(ind+2),0d0,kind=8)
        H(2,2) = cmplx(Q(ind),-Q(ind+1),kind=8)
        R1(1:2,:) = matmul(H,R1(1:2,:))
        
     else
        ! apply first Q
        if (P(STP-2)) then 

           ! inv hess
           ind = 3*STP-5           
           ! update R2 ! inv(Q(1))
           H(1,1) = cmplx(Q(ind),-Q(ind+1),kind=8)
           H(2,1) = cmplx(-Q(ind+2),0d0,kind=8)
           H(1,2) = cmplx(Q(ind+2),0d0,kind=8)
           H(2,2) = cmplx(Q(ind),Q(ind+1),kind=8)
           R2(1:2,:) = matmul(H,R2(1:2,:))
           
        else

           ! hess
           ind = 3*STP-5
           ! update R1 ! Q(1)
           H(1,1) = cmplx(Q(ind),Q(ind+1),kind=8)
           H(2,1) = cmplx(Q(ind+2),0d0,kind=8)
           H(1,2) = cmplx(-Q(ind+2),0d0,kind=8)
           H(2,2) = cmplx(Q(ind),-Q(ind+1),kind=8)
           R1(1:2,:) = matmul(H,R1(1:2,:))

        end if
     end if

  end if

  !print*, "R1", R1(1,:)
  !print*, "R1", R1(2,:)
  !print*, "R1", R1(3,:)
  !print*, "R2", R2(1,:)
  !print*, "R2", R2(2,:)
  !print*, "R2", R2(3,:)
  
  ! store ratio of bottom right entries
  rho = R1(3,3)/R2(3,3) 
  
  ! compute eigenvalues and eigenvectors
  call z_2x2array_eig(.TRUE.,R1(2:3,2:3),R2(2:3,2:3),H,HK)
  
  ! wilkinson shift
  if(abs(R1(3,3)/R2(3,3)-rho) < abs(R1(2,2)/R2(2,2)-rho))then
     SHFT = R1(3,3)/R2(3,3)
  else
     SHFT = R1(2,2)/R2(2,2)
  end if
  
  ! avoid INFs and NANs
  if ((SHFT.NE.SHFT).OR.(abs(SHFT) > EISCOR_DBL_INF)) then
     SHFT = cmplx(0d0,0d0,kind=8) ! not sure if this is a good idea?
  end if
  
  !print*, "shift", shft

end subroutine z_uprkfpen_singleshift
