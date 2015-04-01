#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! d_2x2array_eig 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the real Schur decomposition of a general
! 2x2 double matrix pencil.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  FLAG            LOGICAL
!                    .TRUE. implies B is not the identity
!                    .FALSE. implies B is the identity matrix
!
!  A,B             REAL(8) array of dimension (2,2)
!                    Contains the 2x2 matrix. On exit (A,B) are real and
!                    quasi-uppertriangular. If FLAG=.FALSE. B is
!                    unused.
!
! OUTPUT VARIABLES:
!
!  Q,Z             COMPLEX(8) array of dimension (2,2)
!                    On exit the columns of Q and Z contain the left and right
!                    Schur vectors of the pencil (A,B). If FLAG=.FALSE. Z is
!                    unused.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine d_2x2array_eig(FLAG,A,B,Q,Z)
  
  implicit none
  
  ! input variables
  logical, intent(inout) :: FLAG
  real(8), intent(inout) :: A(2,2), B(2,2), Q(2,2), Z(2,2)
  
  ! compute variables
  integer :: ii, id
  real(8) :: trace, disc, detm, temp
  
  ! B not identity
  if (FLAG) then
  
  
  ! B identity
  else  

    ! compute intermediate values
    trace = A(1,1) + A(2,2)
    detm = A(1,1)*A(2,2) - A(2,1)*A(1,2)
    temp = (A(1,1)-A(2,2))**2 + 4d0*A(1,2)*A(2,1)
    
    ! imaginary roots
    if (temp < 0) then
      
      ! move A to standard form (A(1,1) = A(2,2))
      ! real part of lambda
      trace = trace/2d0
    
      ! imaginary part of lambda 
      disc = sqrt(-temp)/2d0
   
      ! compute Q
      temp = (A(2,2)-A(1,1))
      if ((temp.NE.0).AND.(abs(temp)>1)) then
        temp = (A(1,2)+A(2,1))/temp
        temp = temp*(1d0-sqrt(1d0+1d0/temp/temp))
      else if (temp.NE.0) then 
        temp = (A(1,2)+A(2,1))/temp
        temp = temp-sqrt(1d0+temp*temp)
      end if
      call d_rot2_vec2gen(1d0,temp,Q(1,1),Q(2,1),temp)
      Q(2,2) = Q(1,1)
      Q(1,2) = -Q(2,1)
    
      ! update A
      A = matmul(transpose(Q),matmul(A,Q))

    ! real roots
    else
     
      ! sqrt of discriminant 
      disc = sqrt(temp)
      
      ! compute most accurate eigenvalue
      if(abs(trace+disc) > abs(trace-disc))then
        temp = (trace+disc)/2d0
      else
        temp = (trace-disc)/2d0
      end if

      ! compute Schur vectors
      call d_rot2_vec2gen(A(1,2),temp-A(1,1),Q(1,1),Q(2,1),trace)
      Q(2,2) = Q(1,1)
      Q(1,2) = -Q(2,1)
      
      ! update A
      A(1,1) = temp
      A(1,2) = 0d0
      A(2,2) = detm/temp
      A(2,1) = 0d0
      
    end if
  
  end if
  
end subroutine d_2x2array_eig
