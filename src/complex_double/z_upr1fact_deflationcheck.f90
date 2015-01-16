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
!  STR             INTEGER
!                    index of the top most givens rotation where 
!                    deflation should be checked
!
!  STP             INTEGER
!                    index of the bottom most givens rotation where 
!                    deflation should be checked
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!                    generators must be orthogonal to working precision
!
!  D               REAL(8) array of dimension (2*(N+1))
!                    array of generators for complex diagonal matrix
!                    generators must be orthogonal to working precision
!
!  ITCNT           INTEGER
!                    number of iterations since last deflation
!
! OUTPUT VARIABLES:
!
!  ZERO            INTEGER
!                     index of the last known deflation
!                     on output contains index of newest deflation
!
!  ITS             INTEGER array of dimension (N-1)
!                    Contains the number of iterations per deflation
!
!  INFO            INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies N is invalid
!                   INFO = -2 implies STR is invalid
!                   INFO = -3 implies STP is invalid
!                   INFO = -4 implies ZERO is invalid
!                   INFO = -6 implies Q is invalid
!                   INFO = -7 implies D is invalid
!                   INFO = -8 implies ITCNT is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fact_deflationcheck(N,STR,STP,ZERO,P,Q,D,ITCNT,ITS,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N, STP
  integer, intent(inout) :: STR, ZERO, ITCNT, INFO, ITS(N-1)
  logical, intent(in) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D(2*(N+1))

  ! compute variables
  integer :: ii, jj, ind, up, down
  real(8), parameter :: tol = epsilon(1d0)
  real(8) :: qr, qi, dr, di, cr, ci, s, nrm
  
  ! initialize info
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
    
    ! check N
    if (N < 2) then
      INFO = -1
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
      return
    end if
    
    ! check STR
    if ((STR < 1).OR.(STR > N-1)) then
      INFO = -2
      call u_infocode_check(__FILE__,__LINE__,"STR must 1 <= STR <= N-1",INFO,INFO)
      return
    end if 
    
    ! check STP
    if ((STP < STR).OR.(STP > N-1)) then
      INFO = -3
      call u_infocode_check(__FILE__,__LINE__,"STP must STR <= STP <= N-1",INFO,INFO)
      return
    end if  
    
    ! check ZERO
    if ((ZERO >= STR).OR.(ZERO < 0)) then
      INFO = -4
      call u_infocode_check(__FILE__,__LINE__,"ZERO must 0 <= ZERO < STR",INFO,INFO)
      return
    end if  
    
    ! check Q
    call d_1Darray_check(3*(N-1),Q,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,-6)
      return
    end if

    ! check D
    call d_1Darray_check(2*(N+1),D,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,-7)
      return
    end if

    call z_unifact_factorcheck(N,Q,D(1:2*N),INFO)
    if (INFO.EQ.-1) then
      call u_infocode_check(__FILE__,__LINE__,"N must be at least 2",INFO,INFO)
      return
    end if
    if (INFO.EQ.-2) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,-6)
      return
    end if
    if (INFO.EQ.-3) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,-7)
      return
    end if
    
    ! check ITCNT
    if (ITCNT < 0) then
      INFO = -8
      call u_infocode_check(__FILE__,__LINE__,"ITCNT must be non-negative",INFO,INFO)
      return
    end if 

  end if
  
  ! check for deflation
  do ii=1,(STP-STR+1)
  
    ! index of the rotaion being checked
    ind = (STP-ii+1)
   
    ! deflate if subdiagonal is small enough
    nrm = abs(Q(3*ind))
    if(nrm < tol)then
      
      ! extract diagonal
      qr = Q(3*ind-2)
      qi = Q(3*ind-1)
                
      ! set rotation to identity
      Q(3*ind-2) = 1d0
      Q(3*ind-1) = 0d0
      Q(3*ind) = 0d0
      
      ! initialize up
      up = ind
        
      ! deflate upward
      do jj = 1,(ind-STR)
        
        ! set upward index
        up = ind-jj
        
        ! exit loop if P == .FALSE.
        if (P(up).EQV..FALSE.) then
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
        
        call z_rot3_vec3gen(cr,ci,s,Q(3*up-2),Q(3*up-1),Q(3*up),nrm,INFO)
        
        ! check INFO in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec3gen failed",INFO,INFO)
          if (INFO.NE.0) then 
            return 
          end if 
        end if
        
      end do
       
      ! update upward diagonal
      dr = D(2*up-1)
      di = D(2*up)
            
      nrm = qr*dr - qi*di
      di = qr*di + qi*dr
      dr = nrm

      call d_rot2_vec2gen(dr,di,D(2*up-1),D(2*up),nrm,INFO)

      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if    
    
      ! initialize downward index
      down = ind
        
      ! deflate downward
      do jj = 1,(STP-ind)
        
        ! exit if P == .TRUE.
        if (P(down).EQV..TRUE.) then
          exit
        end if
                
        ! set downward index
        down = ind + jj
        
        ! update Q
        cr = Q(3*down-2)
        ci = Q(3*down-1)
        s = Q(3*down)
                
        nrm = qr*cr + qi*ci
        ci = qr*ci - qi*cr
        cr = nrm

        call z_rot3_vec3gen(cr,ci,s,Q(3*down-2),Q(3*down-1),Q(3*down),nrm,INFO)
        
        ! check INFO in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec3gen failed",INFO,INFO)
          if (INFO.NE.0) then 
            return 
          end if 
        end if
                    
      end do
      
      ! update downward index
      down = down + 1
      
      ! update downward diagonal
      dr = D(2*down-1)
      di = D(2*down)
            
      nrm = qr*dr + qi*di
      di = qr*di - qi*dr
      dr = nrm
      
      call d_rot2_vec2gen(dr,di,D(2*down-1),D(2*down),nrm,INFO)

      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"d_rot2_vec2gen failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if  
      
      ! update indices
      ZERO = STP+1-ii
      STR = ZERO + 1
    
      ! store it_count
      ITS(ZERO) = ITCNT
      ITCNT = 0
          
      ! exit loop 
      exit
        
    end if
    
  end do
  
end subroutine z_upr1fact_deflationcheck
