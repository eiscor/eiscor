#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_symtrid_qr_1dlaplace
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the eigenvalues of the tridiagonal matrix
! [ -1 2 -1 ]. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_symtrid_qr_race

  implicit none
  
  ! compute variables
  integer, parameter :: N1 = 512
  integer, parameter :: N2 = 512
  real(8), parameter :: scale = 1d0
  !logical, parameter :: backward = .TRUE.
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: N, M, MM, N3
  integer :: ii, jj, kk, ll, ij, INFO, IWORK(3+5*N2)
  real(8) :: WORK(14*N2+N2*N2), D(N2), E(N2), t, t1, t2, t3
  real(8) :: Ds(N2), Es(N2), Hr(N2,N2), Zr(N2,N2), pi = EISCOR_DBL_PI
  complex(8) :: Z(N2,N2), H(N2,N2), v(N2), c1
  integer :: ITS(N2-1)
  logical :: backward, random, gauss
  
  ! timing variables
  integer:: c_start, c_start2, c_stop, c_stop2, c_rate
  
  ! BLAS
  double precision :: dnrm2, dznrm2

  random = .TRUE.
  gauss = .TRUE.
  !gauss = .FALSE.
  !random = .FALSE.

  if (.NOT.random) then
     ! initialize T to be a tridiagonal matrix of the form
     !  2 -1
     ! -1  2 -1
     !     -1 2 ...
     Ds = 0d0
     Es = -5d-1
     Es = scale*Es
  else
     call u_fixedseed_initialize(INFO)
     do ii=1,N2
        if (gauss) then
           call d_scalar_random_normal(Ds(ii))
           call d_scalar_random_normal(Es(ii))
         else
           call random_number(t)
           Ds(ii) = t
           call random_number(t)
           Es(ii) = t
        end if
        Ds(ii) = Ds(ii)/3d0
        Es(ii) = Es(ii)/3d0
    end do
  end if
  print*, Ds(1), Ds(N2), EISCOR_DBL_EPS
  
  do ll=2,2
     if (ll.EQ.1) then
        backward = .FALSE.
        print*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print*, "Timings"
     else
        backward = .TRUE.
        print*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print*, "Backward error"
     end if
     
     call system_clock(count_rate=c_rate)
     ! start timer
     call system_clock(count=c_start2)
     
     if (backward) then
        MM = 3
     else
        MM = 4
     end if
     

  do M=1,1!,MM
     ! print banner
     if (M.EQ.1) then
        print*,""
        print*,"example_d_symtrid_qr_1dlaplace:"
        print*,""
        print*, "N ", N1, " ... ", N2, " scale ", scale
        print*,""
     elseif (M.EQ.2) then
       print*,""
       print*,"LAPACK DSTEQR"
       print*,""
     elseif (M.EQ.3) then
        print*,""
        print*,"LAPACK DSTEVD"
        print*,""
     elseif (M.EQ.4) then
         print*,""
         print*,"LAPACK DSTERF"
         print*,""         
     else
        exit
     end if
     N = N1
     do while (N.LE.N2)
        
        if (backward) then
           N3 = 1
        else
           N3 = N2/N
        end if
        
        ! start timer
        call system_clock(count=c_start)
        
        do ij = 1,N3
           ! symtrid_qr
           D = Ds
           E = Es
           if (M.EQ.1) then
              ! call d_orthhess_qr
              if (backward) then
                 call d_symtrid_qr(.TRUE.,.TRUE.,N,D,E,WORK,N,Z,ITS,INFO)
              else
                 call d_symtrid_qr(.FALSE.,.FALSE.,N,D,E,WORK,1,Z,ITS,INFO)
              end if
           elseif (M.EQ.2) then
              ! run LAPACK
              if (backward) then
                 call dsteqr ('I', N, D, E, Zr, N, WORK, INFO)
              else
                 call dsteqr ('N', N, D, E, Z, 1, WORK, INFO)
              end if
           elseif (M.EQ.3) then
              ! run DSTEVD
              if (backward) then
                 call dstevd ('V', N, D, E, Zr, N, WORK, 14*N+N*N, IWORK, 5+5*N, INFO) 
              else
                 call dstevd ('N', N, D, E, Z, 1, WORK, 14*N, IWORK, 3+5*N, INFO) 
              end if
           elseif (M.EQ.4) then
              ! run LAPACK
              call dsterf (N, D, E, INFO)
           end if
        end do

        ! check INFO
        if (INFO.NE.0) then
           print*,"d_symtrid_qr failed."
           print*,"INFO:",INFO
        end if
        
        ! stop timer
        call system_clock(count=c_stop)
        
        if (.NOT.random) then
           ! computing forward error
           do ii=1,N
              E(ii) = 0d0+scale*cos(ii*EISCOR_DBL_PI/(N+1d0))
           end do
           
           t1 = 0d0
           t2 = 0d0
           
           do ii=1,N
              t = abs(D(ii)-E(1))
              INFO = 1
              do jj=2,N
                 if (abs(D(ii)-E(jj))<t) then
                    t = abs(D(ii)-E(jj))
                    INFO = jj
                 end if
              end do
              
              if (t.GT.t2) then
                 t2 = t
                 kk = INFO
              end if
              t1 = t1 + t**2
           end do
           if ((dsqrt(t1)>1d-3).OR.(t1.NE.t1)) then
              call u_test_failed(__LINE__)
           end if
        else 
           t1 = 0d0
           t1 = 0d0/t1
        end if
        
!!$        do ii=1,N
!!$           print*, ii, D(ii)
!!$        end do

        ! computing backward error
        if (backward) then
           if (M==1) then
              t2 = 0d0
              do ii=1,N
                 ! jj  1
                 v(1) = (Ds(1)-D(ii))*Z(1+N*(ii-1),1) + Es(1)*Z(2+N*(ii-1),1)
                 do jj=2,N-1
                    v(jj) = (Ds(jj)-D(ii))*Z(jj+N*(ii-1),1) + &
                         & Es(jj)*Z(jj+1+N*(ii-1),1) + &
                         & Es(jj-1)*Z(jj-1+N*(ii-1),1)
                 end do
                 v(N) = (Ds(N)-D(ii))*Z(N+N*(ii-1),1) + &
                      & Es(N-1)*Z(N-1+N*(ii-1),1)
                 t3 =  dznrm2(N,v,1)
                 if (t3.GT.t2) then
                    t2 = t3
                 end if
              end do
           else
              t2 = 0d0
              do ii=1,N
                 ! jj = 1
                 v(1) = (Ds(1)-D(ii))*cmplx(Zr(1+N*(ii-1),1),0d0,kind=8) + &
                      & Es(1)*cmplx(Zr(2+N*(ii-1),1),0d0,kind=8)
                 do jj=2,N-1
                    v(jj) = (Ds(jj)-D(ii))*cmplx(Zr(jj+N*(ii-1),1),0d0,kind=8) + &
                         & Es(jj)*cmplx(Zr(jj+1+N*(ii-1),1),0d0,kind=8) + &
                         & Es(jj-1)*cmplx(Zr(jj-1+N*(ii-1),1),0d0,kind=8)
                 end do
                 v(N) = (Ds(N)-D(ii))*cmplx(Zr(N+N*(ii-1),1),0d0,kind=8) + &
                      & Es(N-1)*cmplx(Zr(N-1+N*(ii-1),1),0d0,kind=8)
                 t3 =  dznrm2(N,v,1)
                 if (t3.GT.t2) then
                    t2 = t3
                 end if
              end do
           end if
        end if



        if (backward) then
           print*, "(",N, ",",t2,")% for ", dsqrt(t1), " time ", dble(c_stop-c_start)/dble(c_rate)/N3 
        else
           print*, "(",N, ",",dble(c_stop-c_start)/dble(c_rate)/N3,")% for err ", dsqrt(t1)
        end if

        N = 2*N
     end do
  end do
  end do
  
  ! stop timer
  call system_clock(count=c_stop2)
  
  ! print success
  call u_test_passed(dble(c_stop2-c_start2)/dble(c_rate))

end program example_d_symtrid_qr_race
