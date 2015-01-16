#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the eigenvalues and optionally eigenvectors of 
! a unitary upper hessenberg matrix that is stored as a product of 
! givens rotations and a complex diagonal matrix. 
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
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
! OUTPUT VARIABLES:
!
!  Z              COMPLEX(8) array of dimension (N,N)
!                   if COMPZ = 'N' unused
!                   if COMPZ = 'I' stores eigenvectors in Z 
!                   if COMPZ = 'V' update Z to store eigenvectors 
!
!  ITS            INTEGER array of dimension (N-1)
!                   Contains the number of iterations per deflation
!
!  INFO           INTEGER
!                   INFO = 1 implies no convergence
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies COMPZ is invalid
!                   INFO = -2 implies N is invalid
!                   INFO = -3 implies Q is invalid
!                   INFO = -4 implies D is invalid
!                   INFO = -5 implies Z is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_qr(COMPZ,N,Q,D,Z,ITS,INFO)

  implicit none
  
  ! input variables
  character, intent(in) :: COMPZ
  integer, intent(in) :: N
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N)
  integer, intent(inout) :: INFO, ITS(N-1)
  complex(8), intent(inout) :: Z(N,N)
  
  ! compute variables
  integer :: ii, jj, kk, ind1, ind2, ll, strt, k
  integer :: start_index, stop_index, zero_index, it_max, it_count
  real(8) :: nrm
  complex(8) :: block(2,2), eigs(2), temp(2,2)
  
  ! initialize info
  INFO = 0
  
  ! check COMPZ
  if ((COMPZ.NE.'N').AND.(COMPZ.NE.'I').AND.(COMPZ.NE.'V')) then
    INFO = -1
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"COMPZ must be 'N', 'I' or 'V'",INFO,INFO)
    end if
    return
  end if
  
  ! check factorization
  call z_unifact_factorcheck(N,Q,D,INFO)
  if (INFO.EQ.-1) then
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N is invalid",INFO,INFO)
    end if
    INFO = -2
    return
  end if
  if (INFO.EQ.-2) then
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"Q is invalid",INFO,INFO)
    end if
    INFO = -3
    return
  end if
  if (INFO.EQ.-3) then
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"D is invalid",INFO,INFO)
    end if
    INFO = -4
    return
  end if     
  
  ! check Z
  if (COMPZ.EQ.'V') then
    call z_2Darray_check(N,N,Z,INFO)
    if (INFO.NE.0) then
      INFO = -5
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"Z is invalid",INFO,INFO)
      end if
      return
    end if
  end if   
 
  ! initialize storage
  ITS = 0
  
  if (COMPZ.EQ.'I') then
    Z = cmplx(0d0,0d0,kind=8)
    do ii=1,n
      Z(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if
  
  ! initialize indices
  start_index = 1
  stop_index = N-1
  zero_index = 0
  it_max = 20*N
  it_count = 0
  
  ! iteration loop
  do kk=1,it_max
  
    ! check for completion
    if(stop_index <= 0)then
      exit
    end if

    ! check for deflation
    call z_unifact_deflationcheck(N,start_index,stop_index,zero_index,Q,D,it_count,ITS,INFO)
   
    ! check INFO in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"z_unifact_deflationcheck failed",INFO,INFO)
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
    else if(stop_index-1 == zero_index)then
     
      ! get 2x2 block
      call z_unifact_2x2diagblock(N,stop_index,Q,D,block,INFO) 
      
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_unifact_2x2diagblock failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
        
      ! compute eigenvalues and eigenvectors
      call z_2x2array_eig(block,eigs,temp,INFO)
      
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_2x2array_eig failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
      
      ! store eigenvalues
      D(2*(stop_index-1)+1) = dble(eigs(1))
      D(2*(stop_index-1)+2) = aimag(eigs(1))
      D(2*(stop_index-1)+3) = dble(eigs(2))
      D(2*(stop_index-1)+4) = aimag(eigs(2)) 
      
      ! update Q
      Q(3*(stop_index-1)+1) = 1d0
      Q(3*(stop_index-1)+2) = 0d0
      Q(3*(stop_index-1)+3) = 0d0     
           
      ! update eigenvectors
      if ((COMPZ.EQ.'I').OR.(COMPZ.EQ.'V')) then
        Z(:,stop_index:(stop_index+1)) = matmul(Z(:,stop_index:(stop_index+1)),temp)
      end if
           
      ! update indices
      stop_index = stop_index-2
      zero_index = 0
      start_index = 1 
         
    ! if greater than 2x2 chase a bulge and check again
    else
       
      ! chase bulge
      call z_unifact_singlestep(COMPZ,N,start_index,stop_index,Q,D,Z,it_count,INFO)
 
      ! check INFO in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"z_unifact_singlestep failed",INFO,INFO)
        if (INFO.NE.0) then 
          return 
        end if 
      end if
                
    end if
    
    ! if it_max hit
    if (kk == it_max) then
      INFO = 1
      its(stop_index) = it_count
    end if
    
  end do

end subroutine z_unifact_qr
