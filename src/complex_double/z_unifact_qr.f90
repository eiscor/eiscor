#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_qr
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Schur factorization of
! a unitary upper hessenberg matrix that is stored as a product of 
! N-1 Givens rotations and a complex diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvectors
!                    .FALSE.: no schurvectors
!
!  ID              LOGICAL
!                    .TRUE.: initialize to Z to identity
!                    .FALSE.: assume Z is already initialized
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for Givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    on output contains the eigenvalues
!
!  M               INTEGER
!                    leading dimension of Z
!
! OUTPUT VARIABLES:
!
!  Z              COMPLEX(8) array of dimension (M,N)
!                   components of schurvectors
!                   if VEC = .FALSE. unused
!                   if VEC = .TRUE. and ID = .TRUE. initializes Z to I 
!                   if VEC = .TRUE. and ID = .FALSE. assumes Z initialized
!
!  ITS            INTEGER array of dimension (N-1)
!                   contains the number of iterations per deflation
!
!  INFO           INTEGER
!                   INFO = 1 implies no convergence
!                   INFO = 0 implies successful computation
!                   INFO = -3 implies N is invalid
!                   INFO = -4 implies Q is invalid
!                   INFO = -5 implies D is invalid
!                   INFO = -6 implies M is invalid
!                   INFO = -7 implies Z is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_qr(VEC,ID,N,Q,D,M,Z,ITS,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, ID
  integer, intent(in) :: N, M
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N)
  integer, intent(inout) :: INFO, ITS(N-1)
  complex(8), intent(inout) :: Z(M,N)
  
  ! compute variables
  logical :: flg
  integer :: ii, jj, kk, ind1, ind2, ll, strt, k
  integer :: STR, STP, ZERO, ITMAX, ITCNT
  real(8) :: nrm
  real(8) :: Qs(3*(N-1)), Ds(2*N), t2, t3
  complex(8) :: v(N), temp(2,2), Zb(N,M), Z2(N,N)

  ! BLAS
  double precision :: dnrm2, dznrm2
  
  ! initialize info
  INFO = 0
  
  ! check factorization
  call z_unifact_factorcheck(N,Q,D,INFO)
  if (INFO.NE.0) then
    INFO = INFO - 2
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"N, Q, or D is invalid",INFO,INFO)
    end if
    return
  end if
  
  Qs = Q
  Ds = D

  do ii=1,N
     do jj=1,M
        Zb(ii,jj) = conjg(Z(jj,ii))
     end do
  end do

  ! check M
  if (VEC.AND.(M < 1)) then
    INFO = -6
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__,"M must be at least 1",INFO,INFO)
    end if
    return
  end if
  
  ! check Z
  if (VEC.AND..NOT.ID) then
    call z_2Darray_check(M,N,Z,flg)
    if (.NOT.flg) then
      INFO = -7
      ! print error message in debug mode
      if (DEBUG) then
        call u_infocode_check(__FILE__,__LINE__,"Z is invalid",INFO,INFO)
      end if
      return
    end if
  end if   
 
  ! initialize storage
  ITS = 0
  
  if (VEC.AND.ID) then
    Z = cmplx(0d0,0d0,kind=8)
    do ii=1,min(M,N)
      Z(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
    Zb = Z
  end if
  
  ! initialize indices
  STR = 1
  STP = N-1
  ZERO = 0
  ITMAX = 20*N
  ITCNT = 0

  ! iteration loop
  do kk=1,ITMAX

    ! check for completion
    if(STP <= 0)then    
      exit
    end if
    
    ! check for deflation
    call z_unifact_deflationcheck(STP-STR+2,Q((3*STR-2):(3*STP)) &
    ,D((2*STR-1):(2*STP+2)),ZERO)
    
    !print*, STP, STR, ZERO

    ! if 1x1 block remove and check again 
    if(STP == (STR+ZERO-1))then
    
      ! update indices
      ITS(STR+STP-1) = ITCNT
      ITCNT = 0
      STP = STP - 1
      ZERO = 0
      STR = 1
    
    ! if greater than 2x2 chase a bulge
    else

      ! check ZERO
      if (ZERO.GT.0) then
         STR = STR+ZERO
      end if

      ! perform singleshift iteration
      call z_unifact_singlestep(VEC,STP-STR+2,Q((3*STR-2):(3*STP)),D((2*STR-1):(2*STP+2)) &
      ,M,Z(:,STR:(STP+1)),ITCNT)
     
      ! update indices
      ITCNT = ITCNT + 1
 
    end if
    
    ! if ITMAX hit
    if (kk == ITMAX) then
      INFO = 1
      !print*, ITCNT
      ITS(STR+STP-1) = ITCNT
      !do ii=1,N
      !   print*, ITS(ii)
      !end do
      !do ii=1,N
      !   print*, D(2*ii-1), D(2*ii)
      !end do
    end if
    
  end do
  
  !ITCNT = 0
  !do ii=1,N
  !  print*, ITS(ii)
  !  ITCNT = ITCNT + ITS(ii)
  !end do
  !print*, ""
  !print*, "it/N", 1d0*ITCNT/N
  

  ! check backward errors
  if (VEC) then
     Z2 = matmul(Zb,Z)
     t2 = 0d0
     do ii=1,M
        do jj=1,N
           v(jj) = cmplx(Ds(2*jj-1),Ds(2*jj),kind=8)*Z2(jj,ii)
           !print*, jj, v(jj), Z2(jj,ii)
        end do
        !print*, ""
        do jj=N-1,1,-1
           temp(1,1) = cmplx(Qs(3*jj-2),Qs(3*jj-1),kind=8)
           temp(2,1) = cmplx(Qs(3*jj),0d0,kind=8)
           temp(1,2) = -temp(2,1)
           temp(2,2) = conjg(temp(1,1))
           v(jj:jj+1) = matmul(temp,v(jj:jj+1))
        end do

        !do jj=1,N
        !   print*, jj, v(jj)
        !end do
        !print*, ""
        do jj=1,N
           v(jj) = v(jj)-cmplx(D(2*ii-1),D(2*ii),kind=8)*Z2(jj,ii)
           !print*, jj, v(jj)
        end do           
        !print*, ""
        t3 =  dznrm2(N,v,1)

        if (t3.GT.t2) then
           t2 = t3
           !print*, ii, t3, " ev", D(2*ii-1), D(2*ii), "phi-1(ev) ", D(2*ii)/(1d0+D(2*ii-1))
        end if
     end do
  end if
        
  
  !do ii=1,N-1
  !   print*, D(2*ii-1), D(2*ii), ITS(ii), Q(3*ii) 
  !end do
  !print*, D(2*N-1), D(2*N), ITS(N) 
     
end subroutine z_unifact_qr
