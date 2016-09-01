!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkutri_decompress
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine decompresses a unitary plus rank one upper-triangular 
! matrix (uprkutri). Such a matrix is stored as a product 
! of two sequences of Givens' rotations and a complex unimodular 
! diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  DIAG            CHARACTER
!                    .TRUE.: diagonals only
!                    .FALSE.: entire triangular part
!
!  N               INTEGER
!                    dimension of matrix
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  C               REAL(8) array of dimension (3*N)
!                    array of generators for first sequence of rotations
!
!  B               REAL(8) array of dimension (3*N)
!                    array of generators for third sequence of rotations
!
! OUTPUT VARIABLES:
!
!  T               COMPLEX(8) array of dimension (N,N)
!                    if DIAG == .TRUE. diagonal of T stored in first column
!                    if DIAG == .FALSE. entire triangular part is constructed
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkutri_decompress(DIAG,N,K,STR,STP,D,C,B,T)

  implicit none
  
  ! input variables
  logical, intent(in) :: DIAG
  integer, intent(in) :: N, K, STR, STP
  real(8), intent(in) :: D(2*N*K), C(3*N*K), B(3*N*K)
  complex(8), intent(inout) :: T(STP-STR+2,STP-STR+2)
  
  ! compute variables
  integer :: ii, ll, ind, ind2, dind, dind2
  complex(8) :: T2(STP-STR+2,STP-STR+2)
  
  ! diagonals only
  if (DIAG) then

     ! initialize T
     T(:,1) = cmplx(0d0,0d0,kind=8)

     ! compute diagonals
     ind = STR*3-2
     ind2 = STP*3+3
     dind = STR*2-1
     dind2 = STP*2+2

     call z_upr1utri_decompress(.TRUE.,STP-STR+2,&
          &D(dind:dind2),C(ind:ind2),B(ind:ind2),T)

     do ll = 2,K
        ind   = 3*N*(ll-1) + STR*3-2
        ind2  = 3*N*(ll-1) + STP*3+3
        dind  = 2*N*(ll-1) + STR*2-1
        dind2 = 2*N*(ll-1) + STP*2+2
        
        call z_upr1utri_decompress(.TRUE.,STP-STR+2,&
             &D(dind:dind2),C(ind:ind2),B(ind:ind2),T2)
        do ii = 1, STP-STR+2
           T(ii,1) = T(ii,1)*T2(ii,1)
        end do
     end do
   
  ! entire upper-triangular part
  else

     ! initialize T
     !T = cmplx(0d0,0d0,kind=8)
     !do ii=1,STP-STR+2
     !   T(ii,ii) = cmplx(1d0,0d0,kind=8)
     !end do

     ! compute diagonals
     ind = STR*3-2
     ind2 = STP*3+3
     dind = STR*2-1
     dind2 = STP*2+2

     call z_upr1utri_decompress(.FALSE.,STP-STR+2,&
          &D(dind:dind2),C(ind:ind2),B(ind:ind2),T)
     do ll = 2,K
        ind   = 3*N*(ll-1) + STR*3-2
        ind2  = 3*N*(ll-1) + STP*3+3
        dind  = 2*N*(ll-1) + STR*2-1
        dind2 = 2*N*(ll-1) + STP*2+2
        
        call z_upr1utri_decompress(.FALSE.,STP-STR+2,&
             &D(dind:dind2),C(ind:ind2),B(ind:ind2),T2)

        T = matmul(T,T2)
     end do
  
  end if  
  
end subroutine z_uprkutri_decompress
