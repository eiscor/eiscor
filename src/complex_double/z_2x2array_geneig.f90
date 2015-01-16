#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_2x2array_geneig
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generalized Schur decomposition of a 
! 2x2 matrix pencil (A,B). On exit the pencil (A,B) will be in upper-
! triangular form.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  JOB            CHARACTER
!                    'S': assumes B = I (Standard)
!                    'G': makes no assumptions on B (General)
!
!  A, B           COMPLEX(8) array of dimension (2,2)
!                   The 2x2 matrix pencil. Upper-triangular on exit.
!
! OUTPUT VARIABLES:
!
!  Q, Z           COMPLEX(8) array of dimension (2,2)
!                   On exit unitary transformations such that
!                   Q*(A,B)Z is upper-trinagular
!
! INFO            INTEGER
!                   INFO = 1 subroutine failed
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies JOB is invalid
!                   INFO = -2 implies A is invalid
!                   INFO = -3 implies B is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_2x2array_geneig(JOB,A,B,Q,Z,INFO)
  
  implicit none
  
  ! input variables
  character, intent(in) :: JOB
  integer, intent(inout) :: INFO
  complex(8), intent(inout) :: A(2,2), B(2,2), Q(2,2), Z(2,2)
  
  ! compute variables
  integer :: ii
  real(8) :: nrm, cr, ci, s
  complex(8) :: temp(2,2), E(2), p0, p1, p2, lambda, disc
  
  ! initialize info
  INFO = 0
  
  ! check input in debug mode
  if (DEBUG) then
  
    ! check JOB
    if ((JOB.NE.'S').AND.(JOB.NE.'G')) then
      INFO = -1
      call u_infocode_check(__FILE__,__LINE__,"JOB must be 'S' or 'G'",INFO,INFO)
      return
    end if
  
    ! check A
    call z_2Darray_check(2,2,A,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"A is invalid",INFO,-2)
      return
    end if
    
    ! check B
    call z_2Darray_check(2,2,B,INFO)
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"B is invalid",INFO,-3)
      return
    end if
  
  end if
  
  ! call z_2x2array_eig if JOB == S
  if (JOB.EQ.'S') then
  
    ! call z_2x2array_eig
    temp = A
    call z_2x2array_eig(temp,E,Q,INFO)
    
    ! check info in debug mode
    if (DEBUG) then
      if (INFO.NE.0) then
        call u_infocode_check(__FILE__,__LINE__,"z_2x2array_eig",INFO,1)
        return
      end if
    end if
    
    ! set A
    A = matmul(A,Q)
    A = matmul(transpose(conjg(Q)),A)
    A(2,1) = cmplx(0d0,0d0,kind=8)
  
    ! set B
    B(1,1) = cmplx(1d0,0d0,kind=8)
    B(2,1) = cmplx(0d0,0d0,kind=8)
    B(2,2) = cmplx(1d0,0d0,kind=8)
    B(1,2) = cmplx(0d0,0d0,kind=8)
    
    ! set Z
    Z = Q
    
    ! return
    return
    
  end if
  
  ! make B upper-triangular
  if (abs(B(2,1)).NE.0d0) then
 
    ! choose row/column with best relative scale
    ! first column is best
    if (abs(abs(B(1,1))-abs(B(2,1))) < abs(abs(B(2,2))-abs(B(2,1)))) then
    
      ! create eliminator
      call z_rot3_vec4gen(dble(B(1,1)),aimag(B(1,1)),dble(B(2,1)),aimag(B(2,1)),cr,ci,s,nrm,INFO)
      
      ! check info in debug mode
      if (DEBUG) then
        if (INFO.NE.0) then
          call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec4gen failed",INFO,1)
          return
        end if
      end if
      
      ! store in Q
      Q(1,1) = cmplx(cr,ci,kind=8)
      Q(2,1) = cmplx(s,0d0,kind=8)
      Q(2,2) = conjg(Q(1,1))
      Q(1,2) = -Q(2,1)
      
      ! update B
      B = matmul(conjg(transpose(Q)),B)
      B(2,1) = cmplx(0d0,0d0,kind=8)
      
      ! update A
      A = matmul(conjg(transpose(Q)),A)
      
      ! set Z
      Z(1,1) = cmplx(1d0,0d0,kind=8)
      Z(2,1) = cmplx(0d0,0d0,kind=8)
      Z(2,2) = cmplx(1d0,0d0,kind=8)
      Z(1,2) = cmplx(0d0,0d0,kind=8)
    
    ! second row is best
    else
    
      ! create eliminator
      call z_rot3_vec4gen(dble(B(2,2)),aimag(B(2,2)),dble(B(2,1)),aimag(B(2,1)),cr,ci,s,nrm,INFO)
      
      ! check info in debug mode
      if (DEBUG) then
        if (INFO.NE.0) then
          call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec4gen failed",INFO,1)
          return
        end if
      end if
      
      ! store in Z
      Z(1,1) = cmplx(cr,-ci,kind=8)
      Z(2,1) = cmplx(s,0d0,kind=8)
      Z(2,2) = conjg(Z(1,1))
      Z(1,2) = -Z(2,1)
      
      ! update B
      B = matmul(B,Z)
      B(2,1) = cmplx(0d0,0d0,kind=8)
      
      ! update A
      A = matmul(A,Z)
      
      ! set Q
      Q(1,1) = cmplx(1d0,0d0,kind=8)
      Q(2,1) = cmplx(0d0,0d0,kind=8)
      Q(2,2) = cmplx(1d0,0d0,kind=8)
      Q(1,2) = cmplx(0d0,0d0,kind=8)
      
    end if
    
  ! set Q and Z to identity if already upper-triangular  
  else
  
    ! set Q
    Q(1,1) = cmplx(1d0,0d0,kind=8)
    Q(2,1) = cmplx(0d0,0d0,kind=8)
    Q(2,2) = cmplx(1d0,0d0,kind=8)
    Q(1,2) = cmplx(0d0,0d0,kind=8)
      
    ! set Z
    Z = Q
  
  end if
  
  ! check for infinite eigenvalues
  if (abs(B(1,1)).EQ.0d0) then
  
      ! create eliminator
      call z_rot3_vec4gen(dble(A(1,1)),aimag(A(1,1)),dble(A(2,1)),aimag(A(2,1)),cr,ci,s,nrm,INFO)
      
      ! check info in debug mode
      if (DEBUG) then
        if (INFO.NE.0) then
          call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec4gen failed",INFO,1)
          return
        end if
      end if
      
      ! store in temp
      temp(1,1) = cmplx(cr,ci,kind=8)
      temp(2,1) = cmplx(s,0d0,kind=8)
      temp(2,2) = conjg(temp(1,1))
      temp(1,2) = -temp(2,1)
      
      ! update A
      A = matmul(conjg(transpose(temp)),A)
      A(2,1) = cmplx(0d0,0d0,kind=8)
      
      ! update B
      B = matmul(conjg(transpose(Q)),B)
      
      ! update Q
      Q = matmul(Q,temp)
      
      ! return
      return
  
  end if
  
  ! check for infinite eigenvalues
  if (abs(B(2,2)).EQ.0d0) then
  
      ! create eliminator
      call z_rot3_vec4gen(dble(A(2,2)),aimag(A(2,2)),dble(A(2,1)),aimag(A(2,1)),cr,ci,s,nrm,INFO)
      
      ! check info in debug mode
      if (DEBUG) then
        if (INFO.NE.0) then
          call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec4gen failed",INFO,1)
          return
        end if
      end if
      
      ! store in temp
      temp(1,1) = cmplx(cr,-ci,kind=8)
      temp(2,1) = cmplx(s,0d0,kind=8)
      temp(2,2) = conjg(temp(1,1))
      temp(1,2) = -temp(2,1)
      
      ! update A
      A = matmul(A,temp)
      A(2,1) = cmplx(0d0,0d0,kind=8)
      
      ! update B
      B = matmul(B,temp)
      
      ! update Z
      Z = matmul(Z,temp)
      
      ! return
      return
  
  end if
  
  ! compute polynomial coefficients
  p2 = B(1,1)*B(2,2)
  p1 = A(2,1)*B(1,2) - (A(1,1)*B(2,2)+A(2,2)*B(1,1))
  p0 = A(1,1)*A(2,2) - A(2,1)*A(1,2)  
  
  ! compute intermediate values
  disc = sqrt(p1**2 - 4d0*p2*p0)
  
  ! compute most accurate lambda
  if(abs(-p1+disc) > abs(-p1-disc))then
     lambda = (-p1+disc)/2d0/p2
  else
     lambda = (-p1-disc)/2d0/p2
  end if
  
  ! compute right eigenvector 
  temp = A - lambda*B

  ! create eliminator
  call z_rot3_vec4gen(dble(temp(1,1)),-aimag(temp(1,1)),dble(temp(1,2)),-aimag(temp(1,2)),cr,ci,s,nrm,INFO)
      
  ! check info in debug mode
  if (DEBUG) then
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec4gen failed",INFO,1)
      return
    end if
  end if
      
  ! store in temp
  temp(1,1) = cmplx(-s,0d0,kind=8)
  temp(2,1) = cmplx(cr,-ci,kind=8)
  temp(1,2) = cmplx(cr,ci,kind=8)
  temp(2,2) = cmplx(s,0d0,kind=8)
  
  ! update A
  A = matmul(A,temp)
  
  ! update B
  B = matmul(B,temp)
  
  ! update Z
  Z = matmul(Z,temp)
  
  ! return to upper-triangular form  
  ! create eliminator
  call z_rot3_vec4gen(dble(A(1,1)),aimag(A(1,1)),dble(A(2,1)),aimag(A(2,1)),cr,ci,s,nrm,INFO)
      
  ! check info in debug mode
  if (DEBUG) then
    if (INFO.NE.0) then
      call u_infocode_check(__FILE__,__LINE__,"z_rot3_vec4gen failed",INFO,1)
      return
    end if
  end if
      
  ! store in temp
  temp(1,1) = cmplx(cr,ci,kind=8)
  temp(2,1) = cmplx(s,0d0,kind=8)
  temp(1,2) = cmplx(-s,0d0,kind=8)
  temp(2,2) = cmplx(cr,-ci,kind=8)
  
  ! update A
  A = matmul(conjg(transpose(temp)),A)
  A(2,1) = cmplx(0d0,0d0,kind=8)
  
  ! update B
  B = matmul(conjg(transpose(temp)),B)
  B(2,1) = cmplx(0d0,0d0,kind=8)
  
  ! update Q 
  Q = matmul(Q,temp)
  
end subroutine z_2x2array_geneig
