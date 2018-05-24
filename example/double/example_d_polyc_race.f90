#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_polyc_race
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the roots of a polynomial in Chebyshev basis.
! The polynomials are of dimension N1 ... N2.
!
! type 1: random coefficients
! type 2: T_n(x) - 1d0  = 0d0
! type 3: T_n(x) - 1d-1 = 0d0
! type 4: \sum  i T_i(x) = 0d0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_polyc_race

  implicit none
  
  ! compute variables
  integer, parameter :: ttype = 4
  integer, parameter :: ires = 3
  integer, parameter :: NEWTONSTEPS = 0
  integer, parameter :: N1 = 4
  integer, parameter :: isolv_max = 3
  !integer, parameter :: N2 = 128
  integer, parameter :: N2 = 1024
  !integer, parameter :: N2 = 4096
  !integer, parameter :: N2 = 16384
  !integer, parameter :: N2 = 4
  integer :: N, N3, ii, jj, ij, isolv , info, type
  integer, allocatable :: ITS(:)
  real(8), allocatable :: COEFFS(:), RESIDUALS(:), RES(:,:), RRECUR(:,:)
  real(8), allocatable :: RROOTS(:), IROOTS(:), PCOEFFS(:)
  real(8) :: a, b
  complex(8), allocatable :: ROOTS(:)
  complex(8), allocatable :: CCOEFFS(:), RECUR(:,:), ALLROOTS (:,:)

  ! timing variables
  integer:: c_start, c_stop, c_rate
  integer:: c_str2, c_stp2

  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)
  
  ! print banner
  print*,""
  print*,"example_d_polyc_race:"
  print*,""

  
  do type = 1,ttype

    ! output
    open (unit=21, file='tim_eiscor.txt', status='unknown', position="append")
    open (unit=31, file='err_eiscor.txt', status='unknown', position="append")
    open (unit=22, file='tim_avw.txt', status='unknown', position="append")
    open (unit=32, file='err_avw.txt', status='unknown', position="append")
    open (unit=23, file='tim_eiscor_scale.txt', status='unknown', position="append")
    open (unit=33, file='err_eiscor_scale.txt', status='unknown', position="append")
  
    
    ! allocate memory
    allocate(COEFFS(N2),RESIDUALS(N2),ROOTS(N2),RROOTS(N2),IROOTS(N2))
    allocate(CCOEFFS(N2),ITS(N2),PCOEFFS(N2+1))
    !,C(N+3))
    

    COEFFS = 0d0

  
  select case (type)
  case (1) 
     do isolv=1,isolv_max
        write (20+isolv,*) "% type 1: random coefficients"
        write (30+isolv,*) "% type 1: random coefficients"
        write (*,*) "% type 1: random coefficients"
        write (*,*) "% type 1: random coefficients"
     end do

  case (2) 
     do isolv=1,isolv_max
        write (20+isolv,*) "% type 2: T_n(x) - 1d0  = 0d0"
        write (30+isolv,*) "% type 2: T_n(x) - 1d0  = 0d0"
        write (*,*) "% type 2: T_n(x) - 1d0  = 0d0"
        write (*,*) "% type 2: T_n(x) - 1d0  = 0d0"
     end do
     
  case (3) 
     do isolv=1,isolv_max
        write (20+isolv,*) "% type 3: T_n(x) - 1d-1 = 0d0"
        write (30+isolv,*) "% type 3: T_n(x) - 1d-1 = 0d0"
        write (*,*) "% type 3: T_n(x) - 1d-1 = 0d0"
        write (*,*) "% type 3: T_n(x) - 1d-1 = 0d0"
     end do
     
  case (4) 
     do isolv=1,isolv_max
        write (20+isolv,*) "% type 4: \sum  i T_i(x) = 0d0"
        write (30+isolv,*) "% type 4: \sum  i T_i(x) = 0d0"
        write (*,*) "% type 4: \sum  i T_i(x) = 0d0"
        write (*,*) "% type 4: \sum  i T_i(x) = 0d0"
     end do

  end select

  do isolv=1,isolv_max
     select case (isolv)
     case (1)
        write (20+isolv,*) "\addplot coordinates{ % eiscor"
        write (30+isolv,*) "\addplot coordinates{ % eiscor"

     case (2)
        write (20+isolv,*) "\addplot coordinates{ % avw"
        write (30+isolv,*) "\addplot coordinates{ % avw"

     case (3)
        write (20+isolv,*) "\addplot coordinates{ % eiscor scale"
        write (30+isolv,*) "\addplot coordinates{ % eiscor scale"

     end select
  end do

  N = N1
  do while (N.LE.N2)
     allocate(RECUR(N,3),ALLROOTS(N,(NEWTONSTEPS+1)),RES(N,(3*(NEWTONSTEPS+1))))
     allocate(RRECUR(N,3))

     select case (type)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! type 1: random coefficients
     case (1) 
        call d_1Darray_random_normal(N,COEFFS)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! type 2: T_n(x) - 1d0  = 0d0
     case (2) 
        COEFFS = 0d0        
        COEFFS(N) = 1d0

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! type 3: T_n(x) - 1d-1 = 0d0
     case (3) 
        COEFFS = 0d0        
        COEFFS(N) = 1d-1

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! type 4: \sum  i T_i(x) = 0d0
     case (4) 
        do ii = 1,N
           COEFFS(ii) = 1d0*(N+1-ii)
        end do
        
     end select
     PCOEFFS(1)=1d0
     do ii=1,N
        CCOEFFS(ii) = cmplx(COEFFS(ii),0d0,kind=8)
        PCOEFFS(1+ii) = COEFFS(ii)
     end do

     RECUR = cmplx(0d0,0d0,kind=8)
     RECUR(:,1) = cmplx(.5d0,0d0,kind=8)
     RECUR(:,3) = cmplx(.5d0,0d0,kind=8)
     RECUR(N,1) = cmplx(1d0,0d0,kind=8)
     RECUR(N,2) = cmplx(0d0,0d0,kind=8)    
     RECUR(N,3) = cmplx(0d0,0d0,kind=8)    

     RRECUR = 0d0
     RRECUR(:,1) = .5d0
     RRECUR(:,3) = .5d0
     RRECUR(N,1) = 1d0
     !RRECUR(N,2) = 0d0
     RRECUR(N,3) = 0d0
     
     N3 = N2/N

     !N3 = N3/16
     !if (N3.LT.1) then
     !  N3 = 1
     !end if
     !N3 = 1
     
     do isolv=1,isolv_max
        RES = -1d0
        ! start timer
        call system_clock(count=c_str2)
        
        select case (isolv)
        case (1)
           do ij=1,N3
              call d_polyc_roots(N,COEFFS,ROOTS,RESIDUALS,.TRUE.)
           end do
           
        case (2)
           do ij=1,N3
              call DAVW2(1,N,3,PCOEFFS,RRECUR,RROOTS,IROOTS,ITS,INFO)
           end do

        case (3)
           do ij=1,N3
              call d_polyc_roots(N,COEFFS,ROOTS,RESIDUALS,.FALSE.)
           end do

        end select

        ! stop timer
        call system_clock(count=c_stp2)        
        if (isolv.EQ.2) then
           do ii=1,N
              ROOTS(ii) = cmplx(RROOTS(ii),IROOTS(ii),kind=8)
           end do
        end if

        call z_polyc_residuals(N,3,NEWTONSTEPS,CCOEFFS,RECUR,ROOTS,ALLROOTS,RES)
        
        a = 0d0
        b = 0d0
        do ii=1,N
           if (abs(RES(ii,ires)).GT.a) then
              a = abs(RES(ii,ires))
           end if
           b = b + abs(RES(ii,ires))
        end do
        
        !print*, "max residual", a, "sum of residuals", b
        
        print*, "(",N,",",dble(c_stp2-c_str2)/dble(c_rate)/N3,")%", isolv
        write (20+isolv,*) "(",N,",",dble(c_stp2-c_str2)/dble(c_rate)/N3,")%"
        print*, "(",N,",",a,")%", isolv
        write (30+isolv,*) "(",N,",",a,")%"

     end do

     deallocate(RECUR,ALLROOTS,RES)
     deallocate(RRECUR)
     N = N*2
  end do     
  
  deallocate(COEFFS,RESIDUALS,ROOTS,RROOTS,IROOTS)
  deallocate(CCOEFFS,ITS,PCOEFFS)


  do isolv=1,isolv_max
     write (20+isolv,*) "};"
     write (20+isolv,*) ""
     close(20+isolv)

     write (30+isolv,*) "};"
     write (30+isolv,*) ""
     close(30+isolv)
  end do

end do

  ! stop timer
  call system_clock(count=c_stop)
  print*, "Example took ", dble(c_stop-c_start)/dble(c_rate), "s"

end program example_d_polyc_race
