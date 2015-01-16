#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_orthfact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the eigenvalues and optionally eigen(Schur)vectors of 
! an orthogonal upper hessenberg matrix that is stored as a product of 
! givens rotations and a diagonal matrix.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  COMPZ           CHARACTER
!                    'N': no eigenvectors
!                    'I': eigenvectors, initializing Z to the identity
!                    'V': eigenvectors, assume Z already initialized
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (2*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
! OUTPUT VARIABLES:
!
!  Z               COMPLEX(8) array of dimension (N,N)
!                    if COMPZ = 'N' unused
!                    if COMPZ = 'I' stores eigenvectors in Z 
!                    if COMPZ = 'V' update Z to store eigenvectors 
!
!  ITS             INTEGER array of dimension (N-1)
!                    Contains the number of iterations per deflation
!
!  INFO            INTEGER
!                    INFO >= 2 implies failure in one of the subroutines
!                    INFO = 1 implies failure to converge 
!                    INFO = 0 implies successful computation
!                    INFO = -1 implies COMPZ is invalid
!                    INFO = -2 implies N is invalid
!                    INFO = -3 implies Q is invalid
!                    INFO = -4 implies D is invalid
!                    INFO = -5 implies Z is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_orthfact_qr(COMPZ,N,Q,D,Z,ITS,INFO)
  
  implicit none
  
  ! input variables
  character, intent(in) :: COMPZ
  integer, intent(in) :: N
  integer, intent(inout) :: INFO, ITS(N-1)
  real(8), intent(inout) :: Q(2*N-1), D(2*N), Z(N,N)
  
  ! compute variables
  integer :: ii, jj, kk, ind, m
  integer :: start_index, stop_index, zero_index, it_max, it_count
  real(8) :: block(2,2) 
  complex(8) :: temp(2,2), eigs(2)

  ! initialize INFO
  INFO = 0
  
  ! check COMPZ
  if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
    INFO = -1
      
    ! print error in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"COMPZ must be 'N', 'I' or 'V'",INFO,INFO)
    end if  
      
    return
  end if
  
  ! check factorization
  call d_orthfact_factorcheck(N,Q,D,INFO)
    
  ! print error in debug mode
  if (DEBUG) then
    call u_infocode_check(__FILE__,__LINE__,"Invalid factorization",INFO,INFO)
  end if
  
  if (INFO.EQ.-1) then 
    INFO = -2
    return 
  end if 
  if (INFO.EQ.-2) then 
    INFO = -3
    return 
  end if 
  if (INFO.EQ.-3) then 
    INFO = -4
    return 
  end if 
  
  ! check Z
  if (COMPZ.EQ.'V') then
     
    call d_2Darray_check(N,N,Z,INFO)
     
    ! print error in debug mode
    if (DEBUG) then
       call u_infocode_check(__FILE__,__LINE__,"Z is invalid",INFO,INFO)
    end if
      
    if (INFO.NE.0) then 
       INFO = -5
       return 
    end if
      
  end if
  
  ! initialize storage
  ITS = 0
    
  if (COMPZ.EQ.'I') then
    Z = 0d0
    do ii=1,n
      Z(ii,ii) = 1d0
    end do
  end if
  
  ! initialize indices
  start_index = 1
  stop_index = N-1
  zero_index = 0
  it_max = 20*N
  it_count = 0

  ! loop for bulgechasing
  do kk=1,it_max

    ! check for completion
    if(stop_index <= 0)then
      exit
    end if
        
    ! check for deflation
    call d_orthfact_deflationcheck(N,start_index,stop_index,zero_index,Q,D,it_count,ITS,INFO)
      
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"d_orthfact_deflationcheck failed",INFO,2)
      if (INFO.NE.0) then 
        return 
      end if 
    end if

    ! if 1x1 block remove and check again 
    if(stop_index == zero_index)then
          
      ! update indices
      stop_index = stop_index - 1
      zero_index = 0
      start_index = 1
           
    ! if 2x2 block remove and check again
    else if(stop_index == start_index)then
    
      ! get 2x2 block
      call d_orthfact_2x2diagblock(N,stop_index,Q,D,block,INFO)
        
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"d_orthfact_2x2diagblock failed",INFO,3)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
      
      ! detm > 0
      if ((block(1,1)*block(2,2) - block(2,1)*block(1,2)) > 0) then

        ! store eigenvalues
        D(2*(stop_index-1)+1) = block(1,1)
        D(2*(stop_index-1)+2) = block(1,2)
        D(2*(stop_index-1)+3) = block(2,2)
        D(2*(stop_index-1)+4) = block(2,1)
        
      ! detm < 0
      else
      
        ! compute eigenvalues and eigenvectors
        call d_2x2array_eig(block,eigs,temp,INFO)
          
        ! check INFO in debug mode
        if (DEBUG) then
          call u_infocode_check(__FILE__,__LINE__,"d_2x2array_eig failed",INFO,4)
          if (INFO.NE.0) then 
            return 
          end if 
        end if
        
        ! store eigenvalues
        D(2*(stop_index-1)+1) = sign(1d0,dble(eigs(1)))
        D(2*(stop_index-1)+2) = 0d0
        D(2*(stop_index-1)+3) = sign(1d0,dble(eigs(2)))
        D(2*(stop_index-1)+4) = 0d0
        
        ! update eigenvectors
        if (COMPZ.NE.'N') then

          block(1,1) = dble(temp(1,1))
          block(2,1) = dble(temp(2,1))
          block(1,2) = dble(temp(1,2))
          block(2,2) = dble(temp(2,2))
          
          Z(:,stop_index:(stop_index+1)) = matmul(Z(:,stop_index:(stop_index+1)),block)
        end if
        
      end if
      
      ! update Q
      Q(2*(stop_index-1)+1) = 1d0
      Q(2*(stop_index-1)+2) = 0d0    
             
      ! update indices
      stop_index = stop_index-2
      zero_index = 0
      start_index = 1 
        
    ! if greater than 2x2 chase a bulge and check again
    else
     
      ! chase bulge
      call d_orthfact_doublestep(COMPZ,N,start_index,stop_index,Q,D,Z,it_count,INFO)

      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"d_orthfact_doublestep failed",INFO,5)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
      
    end if
    
    ! if it_max hit
    if (kk == it_max) then
      INFO = 1
      ITS(stop_index) = it_count
      return
    end if
     
  end do 
  
end subroutine d_orthfact_qr
