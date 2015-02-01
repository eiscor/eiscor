#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_deflationcheck 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine checks for deflations in a unitary upper hessenberg 
! matrix that is stored as a product of givens rotations and a complex 
! diagonal matrix. When a deflation occurs the corresponding rotation
! is set to the identity matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!                    generators must be orthogonal to working precision
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    generators must be orthogonal to working precision
!
! OUTPUT VARIABLES:
!
!  ZERO            INTEGER
!                     index of the last known deflation
!                     on output contains index of newest deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_deflationcheck(N,P,Q,D,ZERO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N)
  integer, intent(inout) :: ZERO

  ! compute variables
  integer :: ii, jj, up, down
  real(8), parameter :: tol = EISCOR_DBL_EPS
  real(8) :: qr, qi, dr, di, cr, ci, s, nrm
  
  ! check for deflation
  do ii=1,(N-1)
  
    ! deflate if subdiagonal is small enough
    nrm = abs(Q(3*ii))
    if(nrm < tol)then
      
      ! extract diagonal
      qr = Q(3*ii-2)
      qi = Q(3*ii-1)
                
      ! set rotation to identity
      Q(3*ii-2) = 1d0
      Q(3*ii-1) = 0d0
      Q(3*ii) = 0d0
      
      ! initialize up
      up = ii
        
      ! deflate upward
      do jj = 1,(N-1-ii)
        
        ! set upward index
        up = ii-jj
        
        ! exit loop if P == .FALSE.
        if (.NOT.P(up)) then
          up = up + 1
          exit    
        end if
   
        ! update Q
        cr = Q(3*up-2)
        ci = Q(3*up-1)
        s = Q(3*up)
                
        nrm = qr*cr + qi*ci
        ci = qr*ci - qi*cr
        cr = nrm
        
        call z_rot3_vec3gen(cr,ci,s,Q(3*up-2),Q(3*up-1),Q(3*up),nrm)
        
      end do
       
      ! update upward diagonal
      dr = D(2*up-1)
      di = D(2*up)
            
      nrm = qr*dr - qi*di
      di = qr*di + qi*dr
      dr = nrm

      call d_rot2_vec2gen(dr,di,D(2*up-1),D(2*up),nrm)

      ! initialize downward index
      down = ii
        
      ! deflate downward
      do jj = 1,(N-1-ii)
        
        ! exit if P == .TRUE.
        if (P(down)) then
          exit
        end if
                
        ! set downward index
        down = ii + jj
        
        ! update Q
        cr = Q(3*down-2)
        ci = Q(3*down-1)
        s = Q(3*down)
                
        nrm = qr*cr + qi*ci
        ci = qr*ci - qi*cr
        cr = nrm

        call z_rot3_vec3gen(cr,ci,s,Q(3*down-2),Q(3*down-1),Q(3*down),nrm)
        
      end do
      
      ! update downward index
      down = down + 1
      
      ! update downward diagonal
      dr = D(2*down-1)
      di = D(2*down)
            
      nrm = qr*dr + qi*di
      di = qr*di - qi*dr
      dr = nrm
      
      call d_rot2_vec2gen(dr,di,D(2*down-1),D(2*down),nrm)

      ! update ZERO index
      ZERO = max(0,N-ii)
      
      ! exit loop 
      exit
        
    end if
    
  end do
  
end subroutine z_upr1fact_deflationcheck
