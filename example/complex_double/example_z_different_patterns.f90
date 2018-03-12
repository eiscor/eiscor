#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_z_different_patterns
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_z_different_patterns

  implicit none
  
  ! compute variables
  integer :: N 
  real(8) :: res, forw, err
  integer :: ii, jj, kk, ll, INFO, TYP, noits
  integer, allocatable :: ARGS(:)
  real(8) :: normofp, lambda
  real(8), allocatable :: RESIDUALS(:), FORWARD(:)
  complex(8), allocatable :: COEFFS(:),COEFFSS(:), EROOTS(:), ROOTS(:)
  logical :: eigsknown, EK(48)
  integer, parameter :: N1 = 200
  integer, parameter :: N2 = 200
  
  ! timing variables
  integer:: c_start, c_stop, c_rate
  integer:: c_a2, c_b2
  
  ! z_poly_roots variable
  real(8) :: scl
  logical, allocatable :: P(:)
  integer, allocatable :: ITS(:)
  real(8), allocatable :: Q(:),D1(:),C1(:),B1(:)
  real(8), allocatable :: D2(:),C2(:),B2(:)
  complex(8), allocatable :: V(:),W(:),T(:,:)
  interface
    function l_upr1fact_hess(m,flags)
      logical :: l_upr1fact_hess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_hess
  end interface
  interface
    function l_upr1fact_inversehess(m,flags)
      logical :: l_upr1fact_inversehess
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags

    end function l_upr1fact_inversehess
  end interface
  interface
    function l_upr1fact_cmv(m,flags)
      logical :: l_upr1fact_cmv
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_cmv
  end interface
  interface
    function l_upr1fact_random(m,flags)
      logical :: l_upr1fact_random
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function l_upr1fact_random
  end interface

  
  ! initalize random generator seed
  call u_randomseed_initialize(INFO)
  !call u_fixedseed_initialize(INFO)


  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! print banner
  print*,""
  print*,"example_z_different_patterns"
  print*,"N2", N2
  print*,""

  N = N1

  do while (N.LE.N2) 
     allocate(RESIDUALS(N),FORWARD(N),COEFFS(N+1),COEFFSS(N+1),EROOTS(N),ROOTS(N),ARGS(N))
     allocate(P(N-2),ITS(N-1),Q(3*(N-1)),D1(2*(N+1)),C1(3*N),B1(3*N))    
     allocate(V(N),W(N),D2(2*(N+1)),C2(3*N),B2(3*N))    
     
     
     ! fill coefficients
     call z_1Darray_random_normal(N+1,COEFFS)
     do ii=1,N+1
        !COEFFS(ii) = cmplx(1d0,0d0,kind=8)/COEFFS(ii)
        COEFFS(ii) = cmplx(1d0,0d0,kind=8)
     end do
     COEFFSS=COEFFS
     
     do TYP=1,4
        select case (TYP)
        case (1)
           print*, "Hessenberg"
        case (2)
           print*, "inverse Hessenberg"
        case (3)
           print*, "CMV"
        case (4)
           print*, "random"
        end select
        COEFFS=COEFFSS   

        ! call roots
        call system_clock(count=c_a2)
        
        ! fill P
        P = .FALSE.
        
        ! fill V and W
        scl = maxval(abs(COEFFS))
        V(N) = ((-1d0)**(N))*COEFFS(N+1)/scl
        do ii=1,(N-1)
           V(ii) = -COEFFS(N+1-ii)/scl
        end do
        W = cmplx(0d0,0d0,kind=8)
        W(N) = COEFFS(1)/scl
        
        ! factor companion matrix
        call z_comppen_compress(N,P,V,W,Q,D1,C1,B1,D2,C2,B2)
        
        ! call z_upr1fpen_qz
        select case (TYP)
        case (1)
           call z_upr1fpen_qz(.FALSE.,.FALSE.,l_upr1fact_hess,N,P,Q,D1,C1,B1,D2,C2,B2,N,V,W,ITS,INFO)
        case (2)
           call z_upr1fpen_qz(.FALSE.,.FALSE.,l_upr1fact_inversehess,N,P,Q,D1,C1,B1,D2,C2,B2,N,V,W,ITS,INFO)
        case (3)
           call z_upr1fpen_qz(.FALSE.,.FALSE.,l_upr1fact_cmv,N,P,Q,D1,C1,B1,D2,C2,B2,N,V,W,ITS,INFO)
        case (4)
           call z_upr1fpen_qz(.FALSE.,.FALSE.,l_upr1fact_random,N,P,Q,D1,C1,B1,D2,C2,B2,N,V,W,ITS,INFO)
        end select

        ! extract roots
        call z_upr1utri_decompress(.TRUE.,N,D1,C1,B1,V)
        call z_upr1utri_decompress(.TRUE.,N,D2,C2,B2,W)
        do ii=1,N
           ROOTS(ii) = V(ii)/W(ii)
        end do
        
        ! compute residuals
        call z_poly_residuals(N,COEFFS,ROOTS,0,RESIDUALS)
        
        call system_clock(count=c_b2)
        
        !print*, ROOTS

        noits = 0
        do ii=1,N-1
           noits = noits + ITS(ii)
        end do
        
        
        ! check INFO
        if (INFO.NE.0) then
           !call u_test_failed(__LINE__)
           print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
           print*, "!   INFO not 0                                   !"
           print*, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
           print*, INFO
        end if
        
        ! compute normofp
        normofp = abs(COEFFS(1))**2
        do ii=1,N
           normofp = normofp + abs(COEFFS(ii+1))**2
        end do
        normofp = dsqrt(normofp)
        
        ! maximum residuals
        res = 0d0
        do ii=1,N
           if (RESIDUALS(ii) >= res) then
              res = RESIDUALS(ii)
           end if
        end do
        
        write (*,"(A,1x,I5,1x,A,1x,ES14.4E2,1x,A,ES14.4E2)") "(",N, "&", res, ") \\%", normofp

        write (*,"(A,1x,I5,1x,A,1x,ES14.4E2,1x,A,I6)") "(",N, "&", dble(c_b2-c_a2)/dble(c_rate), ") \\%", noits
        
        ! maximum forward error
        if (eigsknown) then
           forw = 0d0
           ARGS = 0
           do ii = 1,N
              
              err = EISCOR_DBL_INF
              ll = 0
              do jj = 1,N
                 !if ((ARGS(ll).EQ.0).AND.(err .GT. abs(ROOTS(ii)-EROOTS(jj)))) then
                 if (err .GT. abs(ROOTS(ii)-EROOTS(jj))) then
                    err = abs(ROOTS(ii)-EROOTS(jj))
                    ll = jj
                 end if
              end do
              print*, ii, ll, ROOTS(ii), EROOTS(ll), abs(ROOTS(ii)-EROOTS(ll)), err
              
              ARGS(ll) = ARGS(ll)+1
        
              if (err/abs(EROOTS(ll)) .GT. forw) then
                 forw = err/abs(EROOTS(ll))
              end if
           end do
           
        else
        end if
        EK(kk) = eigsknown
        
     end do
     ! free memory
     deallocate(P,ITS,Q,D1,C1,B1,D2,C2,B2,V,W)
     deallocate(RESIDUALS,FORWARD,COEFFS,COEFFSS,EROOTS,ROOTS,ARGS)
     N = N*2
  end do
  
  ! stop timer
  call system_clock(count=c_stop)
  
  ! print success
  print*, "Runtime of this example: ", dble(c_stop-c_start)/dble(c_rate)
  print*,""

end program example_z_different_patterns
