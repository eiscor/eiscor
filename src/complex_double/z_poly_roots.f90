!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_poly_roots 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the roots of a polynomial expressed in the 
! monomial basis using the fast algorithm described in:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    degree of the polynomial
!
!  COEFFS          COMPLEX(8) array of dimension (N+1)
!                    coefficients of polynomial ordered from highest
!                    degree coefficient to lowest degree
!
!  WORK            REAL(8) array of dimension (3*N,4)
!                    array of generators for first sequence of Givens' 
!                    rotations
!
!  ITS             INTEGER array of dimension (N-1)
!                    array that stores the number of iterations per 
!                    computed root
!
! OUTPUT VARIABLES:
!
!  ROOTS           COMPLEX(8) array of dimension (N)
!                    computed roots
!
!  INFO            INTEGER
!                    INFO = 7 sort failed
!                    INFO = 6 root extraction failed
!                    INFO = 5 eigensolver failed
!                    INFO = 4 factorization failed
!                    INFO = 3 implies constant function with no roots
!                    INFO = 2 implies degree is less than N
!                    INFO = 1 implies degree is less than 1
!                    INFO = 0 implies successful computation.
!                    INFO negative implies that an input variable has
!                    an improper value, i.e. INFO=-2 => COEFFS is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_poly_roots(N,COEFFS,ROOTS,WORK,ITS,INFO)

  implicit none
  
  ! input variables
  integer, intent(in) :: N
  complex(8), intent(in) :: COEFFS(N+1)
  complex(8), intent(inout) :: ROOTS(N)
  real(8), intent(inout) :: WORK(3*N,4)
  integer, intent(inout) :: ITS(N-1), INFO
  
  ! compute variables
  integer :: ii, ntmp
  real(8) :: t_stp, t_str
  complex(8) :: Z(1,1)
  
  ! initialize INFO
  INFO = 0
  
  ! prune zero roots
  ITS = 0
  ntmp = N
  do ii=1,N
    if (abs(COEFFS(N+2-ii)).NE.0d0) then
      exit
    end if
    ROOTS(N+1-ii) = cmplx(0d0,0d0,kind=8)
    ntmp = ntmp-1
  end do

  ! return if the constant polynomial
  if (ntmp.EQ.0) then
    INFO = 3
    write(*,*) "Error in "//__FILE__//" line:",__LINE__
    write(*,*) "Polynomial is a non-zero constant with no roots."
    write(*,*) ""
    return
    
  ! if linear factor remains
  else if (ntmp.EQ.1) then
    ROOTS(1) = -COEFFS(2)/COEFFS(1) 
    
  ! use companion QR
  else
  
    ! balance

! start timer
call cpu_time(t_str)    
    ! factor companion matrix
    call z_poly_factorcomp('N',ntmp,COEFFS(2:(ntmp+1))/COEFFS(1),WORK(:,1),WORK(:,2), &
      WORK(:,3),WORK(:,4),Z,INFO)  
    
    ! check INFO
    if (INFO.NE.0) then
      INFO = 4
      write(*,*) "Error in "//__FILE__//" line:",__LINE__
      write(*,*) "ZMBCQR failed."
      write(*,*) ""
      return
    end if
! stop timer
call cpu_time(t_stp)
print*,"time to factor (sec):",t_stp-t_str

! start timer
call cpu_time(t_str)     
    ! call ZPFFQR
    call z_upr1fact_twistedqz('N',ntmp,WORK(:,1),WORK(:,2),WORK(:,3),WORK(:,4),Z,ITS,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      INFO = 5
      write(*,*) "Error in "//__FILE__//" line:",__LINE__
      write(*,*) "ZPFFQR failed."
      write(*,*) ""
      return
    end if
! stop timer
call cpu_time(t_stp)
print*,"time to solve (sec):",t_stp-t_str

! start timer
call cpu_time(t_str)     
    ! extract roots
    call z_upr1fact_extracttri('E',ntmp,WORK(:,1),WORK(:,2),WORK(:,3),WORK(:,4),Z,INFO)
    
    ! check INFO
    if (INFO.NE.0) then
      INFO = 6
      write(*,*) "Error in "//__FILE__//" line:",__LINE__
      write(*,*) "ZPFFET failed."
      write(*,*) ""
      return
    end if
! stop timer
call cpu_time(t_stp)
print*,"time to extract (sec):",t_stp-t_str
    
    ! sort roots by argument
!    call ZARSUE('E',ntmp,WORK(:,2),Z,INFO)
    
    ! check INFO
!    if (INFO.NE.0) then
!      INFO = 7
!      write(*,*) "Error in "//__FILE__//" line:",__LINE__
!      write(*,*) "ZARSUE failed."
!      write(*,*) ""
!      return
!    end if
    
    ! update ROOTS
    do ii=1,ntmp
      ROOTS(ii) = cmplx(WORK(2*(ii-1)+1,2),WORK(2*(ii-1)+2,2),kind=8)
    end do
    
    ! adjust for balancing
  
  end if   

end subroutine z_poly_roots
