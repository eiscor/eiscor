#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthfact_doublestep
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' doubleshift 
! algorithm on an orthogonal upper hessenberg matrix that is stored as a 
! product of givens rotations and a diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: update eigenvectors
!                    .FALSE.: no eigenvectors
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (2*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (N)
!                    array of generators for complex diagonal matrix
!
!  M               INTEGER
!                    leading dimension of Z
!
!  Z               REAL(8) array of dimension (M,N)
!                    if VEC = .TRUE. updated
!                    if VEC = .FALSE. unused
!
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthfact_doublestep(STR,STP,VEC,N,Q,D,M,Z,ITCNT)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: N, M
  integer, intent(in) :: STR, STP
  integer, intent(inout) :: ITCNT
  real(8), intent(inout) :: Q(2*(N-1)), D(N), Z(M,N)
  
  ! compute variables
  integer :: ii, ind
  real(8) :: s1, s2
  real(8) :: b1(2), b2(2), b3(2), temp(2), nrm
  real(8) :: block(2,2) 
  complex(8) :: eigs(2), h(2,2)

  ! compute a nonzero shift
  ! shifts +1
  if((mod(ITCNT+1,11) == 0).AND.(ITCNT.NE.0))then
    eigs(1) = cmplx(1d0,0d0,kind=8)
    eigs(2) = cmplx(1d0,0d0,kind=8)
    
  ! shifts -1
  else if((mod(ITCNT+1,16) == 0).AND.(ITCNT.NE.0))then
    eigs(1) = cmplx(-1d0,0d0,kind=8)
    eigs(2) = cmplx(-1d0,0d0,kind=8)
    
  ! random shift
  else if((mod(ITCNT+1,21) == 0).AND.(ITCNT.NE.0))then
    call random_number(s1)
    call random_number(s2)
    eigs(1) = cmplx(s1,s2,kind=8)
    eigs(2) = cmplx(s1,-s2,kind=8)
          
  ! wilkinson shifts
  else

    ! get 2x2 block
    call d_orthfact_2x2diagblock(.FALSE.,Q((2*N-5):(2*N-2)),D((N-1):N),block)
      
    ! compute eigenvalues and eigenvectors
    call d_2x2array_eig(block,eigs,h)
      
  end if
  
  ! two real shifts
  if ((aimag(eigs(1)).EQ.0d0).AND.(aimag(eigs(2)).EQ.0d0)) then
 
    ! first shift
    eigs(1) = cmplx(sign(1d0,dble(eigs(1))),0d0,kind=8)
    
    ! build first bulge
    temp(1) = dble(eigs(1))
    temp(2) = aimag(eigs(1))
    call d_orthfact_buildbulge('S',N,STR,Q,D,temp,b1,b2)

    ! fusion to initialize first bulge
    ind = 2*(STR-1)
    b3(1) = b1(1)
    b3(2) = -b1(2)
    call d_orthfact_mergebulge('R',b3,Q((ind+1):(ind+2)))
     
    ! first bulge update eigenvectors
    if (VEC)then
      h(1,1) = b1(1)
      h(2,1) = b1(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,STR:(STR+1)) = matmul(Z(:,STR:(STR+1)),h)
    end if 
    
    ! first bulge through D
    call d_rot2_swapdiag('R',D((ind+1):(ind+4)),b1)
      
    ! first bulge through Q
    call d_rot2_turnover(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b1)
      
    ! second shift
    eigs(2) = cmplx(sign(1d0,dble(eigs(2))),0d0,kind=8)
  
    ! build bulge
    temp(1) = dble(eigs(2))
    temp(2) = aimag(eigs(2))
    call d_orthfact_buildbulge('S',N,STR,Q,D,temp,b2,b3)
          
    ! fusion to initialize second bulge
    b3(1) = b2(1)
    b3(2) = -b2(2)
    call d_orthfact_mergebulge('R',b3,Q((ind+1):(ind+2)))
     
    ! main chasing loop
    do ii=STR,(STP-2)
       
      ! set ind
      ind = 2*(ii-1)
       
      ! first bulge update eigenvectors
      if (VEC)then
        h(1,1) = b1(1)
        h(2,1) = b1(2)
        h(1,2) = -h(2,1)
        h(2,2) = h(1,1)
        Z(:,(ii+1):(ii+2)) = matmul(Z(:,(ii+1):(ii+2)),h)
      end if 
        
      ! first bulge through D
      call d_rot2_swapdiag('R',D((ind+3):(ind+6)),b1)
      
      ! first bulge through Q
      call d_rot2_turnover(Q((ind+3):(ind+4)),Q((ind+5):(ind+6)),b1)
      
      ! second bulge update eigenvectors
      if (VEC)then
        h(1,1) = b2(1)
        h(2,1) = b2(2)
        h(1,2) = -h(2,1)
        h(2,2) = h(1,1)
        Z(:,ii:(ii+1)) = matmul(Z(:,ii:(ii+1)),h)
      end if 
        
      ! second bulge through D
      call d_rot2_swapdiag('R',D((ind+1):(ind+4)),b2)
      
      ! second bulge through Q
      call d_rot2_turnover(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b2)
      
    end do
    
    ! set ind
    ind = 2*(stp-2)
    
    ! first bulge update eigenvectors
    if (VEC)then
      h(1,1) = b1(1)
      h(2,1) = b1(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,STP:(STP+1)) = matmul(Z(:,STP:(STP+1)),h)
    end if 
       
    ! first bulge through D
    call d_rot2_swapdiag('R',D((ind+3):(ind+6)),b1)
    
    ! first bulge fuse with Q
    call d_orthfact_mergebulge('L',Q((ind+3):(ind+4)),b1)
    
    ! second bulge update eigenvectors
    if (VEC)then
      h(1,1) = b2(1)
      h(2,1) = b2(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,(STP-1):STP) = matmul(Z(:,(STP-1):STP),h)
    end if
        
    ! second bulge through D
    call d_rot2_swapdiag('R',D((ind+1):(ind+4)),b2)
    
    ! second bulge through Q  
    call d_rot2_turnover(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b2)
    
    ! set index
    ind = 2*(stp-1)
    
    ! last bulge update eigenvectors
    if (VEC)then
      h(1,1) = b2(1)
      h(2,1) = b2(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,STP:(STP+1)) = matmul(Z(:,STP:(STP+1)),h)
    end if 
  
    ! last bulge through D
    call d_rot2_swapdiag('R',D((ind+1):(ind+4)),b2)
    
    ! last bulge fuse with Q
    call d_orthfact_mergebulge('L',Q((ind+1):(ind+2)),b2)
  
    ! update ITCNT
    ITCNT = ITCNT + 1

  ! complex conjugate pair
  else
  
    ! normalize shifts
    call d_rot2_vec2gen(dble(eigs(1)),aimag(eigs(1)),temp(1),temp(2),nrm) 

    ! build bulge
    call d_orthfact_buildbulge('D',N,STR,Q,D,temp,b1,b2)

    ! turnover to initialize bulge
    ind = 2*(STR-1)
    temp(1) = b2(1)
    temp(2) = -b2(2)
    b3(1) = b1(1)
    b3(2) = -b1(2)
    call d_rot2_turnover(temp,b3,Q((ind+1):(ind+2)))
      
    ! fusion to finish initialization
    call d_orthfact_mergebulge('R',b3,Q((ind+3):(ind+4)))
     
    ! update b3
    b3 = Q((ind+1):(ind+2))
    
    ! update Q
    Q((ind+1):(ind+2)) = temp
    
    ! main chasing loop
    do ii=str,(stp-2)
       
      ! set ind
      ind = 2*(ii-1)
       
      ! first bulge update eigenvectors
      if (VEC)then
        h(1,1) = b1(1)
        h(2,1) = b1(2)
        h(1,2) = -h(2,1)
        h(2,2) = h(1,1)
        Z(:,(ii+1):(ii+2)) = matmul(Z(:,(ii+1):(ii+2)),h)
      end if 
        
      ! first bulge through D
      call d_rot2_swapdiag('R',D((ind+3):(ind+6)),b1)
      
      ! first bulge through Q
      call d_rot2_turnover(Q((ind+3):(ind+4)),Q((ind+5):(ind+6)),b1)
      
      ! second bulge update eigenvectors
      if (VEC)then
        h(1,1) = b2(1)
        h(2,1) = b2(2)
        h(1,2) = -h(2,1)
        h(2,2) = h(1,1)
        Z(:,ii:(ii+1)) = matmul(Z(:,ii:(ii+1)),h)
      end if 
        
      ! second bulge through D
      call d_rot2_swapdiag('R',D((ind+1):(ind+4)),b2)
      
      ! second bulge through Q
      call d_rot2_turnover(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b2)
      
      ! push b3 down
      call d_rot2_turnover(b3,b1,b2)
      
    ! update bulges
      temp = b2
      b2 = b3
      b3 = b1
      b1 = temp
       
    end do
    
    ! set ind
    ind = 2*(stp-2)
    
    ! first bulge update eigenvectors
    if (VEC)then
      h(1,1) = b1(1)
      h(2,1) = b1(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,STP:(STP+1)) = matmul(Z(:,STP:(STP+1)),h)
    end if 
       
    ! first bulge through D
    call d_rot2_swapdiag('R',D((ind+3):(ind+6)),b1)
    
    ! first bulge fuse with Q
    call d_orthfact_mergebulge('L',Q((ind+3):(ind+4)),b1)
    
    ! second bulge update eigenvectors
    if (VEC)then
      h(1,1) = b2(1)
      h(2,1) = b2(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,(STP-1):STP) = matmul(Z(:,(STP-1):STP),h)
    end if
        
    ! second bulge through D
    call d_rot2_swapdiag('R',D((ind+1):(ind+4)),b2)
    
    ! second bulge through Q  
    call d_rot2_turnover(Q((ind+1):(ind+2)),Q((ind+3):(ind+4)),b2)
    
    ! fuse b2 and b3
    call d_orthfact_mergebulge('L',b3,b2)
    
    ! set index
    ind = 2*(stp-1)
    
    ! last bulge update eigenvectors
    if (VEC)then
      h(1,1) = b3(1)
      h(2,1) = b3(2)
      h(1,2) = -h(2,1)
      h(2,2) = h(1,1)
      Z(:,STP:(STP+1)) = matmul(Z(:,STP:(STP+1)),h)
    end if 
  
    ! last bulge through D
    call d_rot2_swapdiag('R',D((ind+1):(ind+4)),b3)
    
    ! last bulge fuse with Q
    call d_orthfact_mergebulge('L',Q((ind+1):(ind+2)),b3)
  
    ! update ITCNT
    ITCNT = ITCNT + 1
  
  end if ! end of doubleshift step

end subroutine d_orthfact_doublestep
