#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_symtrid_qr_race_split
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the eigenvalues of the tridiagonal matrix
! [ -1 2 -1 ]. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_symtrid_qr_race_split

  implicit none
  
  ! compute variables
  integer, parameter :: problem = 1 ! [-0.5, 0, -0.5]  
  !integer, parameter :: problem = 2 ! osipov
  !integer, parameter :: problem = 3 ! random uniform  
  !integer, parameter :: problem = 4 ! random normal
  !integer, parameter :: problem = 5 ! random a exp(10 b), a,b normally distributed
  !integer, parameter :: problem = 6 ! cluster of [-0.5, 0, -0.5] with size=cluster
  integer, parameter :: cluster = 32 
  integer, parameter :: N1 = 2
  integer, parameter :: N2 = 4096
  real(8), parameter :: scale1 = 1d0
  real(8), parameter :: scale2 = 1d0
  real(8), parameter :: shift = 0d0
  logical, parameter :: sca = .TRUE.
  !logical, parameter :: sca = .FALSE.
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer :: N, M, MM, N3
  integer :: ii, jj, kk, ll, ij, INFO, IWORK(3+5*N2)
  real(8) :: WORK(14*N2+N2*N2), D(N2), E(N2), eig(N2), t, t1, t2, t3, nrm, scale, scaleqr
  real(8) :: WORKs(14*N2+N2*N2), Ds(N2), Es(N2), Hr(N2,N2), Zr(N2,N2), pi = EISCOR_DBL_PI
  complex(8) :: Z(N2,N2), H(N2,N2), v(N2), c1
  integer :: ITS(N2-1)
  logical :: backward
  
  ! timing variables
  integer:: c_start, c_start2, c_stop, c_stop2, c_rate, c_mid, c_mid2
  
  ! BLAS
  double precision :: dnrm2, dznrm2

  call system_clock(count_rate=c_rate)
  ! start timer
  call system_clock(count=c_start2)

  scale = scale1
  
  do while (scale.LE.scale2)
     select case (problem)
     case (1) 
        ! initialize T to be a tridiagonal matrix of the form
        !  2 -1
        ! -1  2 -1
        !     -1 2 ...
        Ds = 0d0
        Es = -5d-1
     case (2)
        do ii=1,N2
           Ds(ii) = 2d0 + ii**2/1d6
           Es(ii) = -1d0
        end do
     case (3)  
        call u_randomseed_initialize(INFO)
        !call u_fixedseed_initialize(INFO)
        do ii=1,N2
           call random_number(t)
           Ds(ii) = t
           call random_number(t)
           Es(ii) = t
        end do
     case (4)
        do ii=1,N2
           call d_scalar_random_normal(Ds(ii))
           call d_scalar_random_normal(Es(ii))
        end do
     case (5)
        do ii=1,N2
           call random_number(t)
           call random_number(t1)
           Ds(ii) = t * exp(10*t1)
           call random_number(t)
           call random_number(t1)
           Es(ii) = t * exp(10*t1)
        end do        
     case (6) 
        ! initialize T to be a block tridiagonal matrix of the form
        !  2 -1
        ! -1  2 -1
        !     -1 2  eps
        !        eps 2  -1
        ! ...
        Ds = 0d0
        Es = -5d-1
        do ii=cluster-1,N2,cluster
           Es(ii) = 1e-10
        end do
     end select

     do ii=1,N2
        Ds(ii) = (Ds(ii)-shift)*scale
        Es(ii) = Es(ii)*scale
     end do

     print*, Ds(1), Ds(N2), Es(1), EISCOR_DBL_EPS
  
  do ll=1,2
     if (ll.EQ.1) then
        backward = .FALSE.
        print*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print*, "Timings"
     else
        backward = .TRUE.
        print*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
        print*, "Backward error"
     end if
          
     if (backward) then
        MM = 3
     else
        MM = 4
     end if
     

  do M=1,4
    if (M.EQ.1) then
        print*,""
        print*,"example_d_symtrid_qr_1dlaplace:"
        print*,""
        print*, "N", N1, " ... ", N2, "problem", problem, "scale", scale, "shift", shift
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
                 call d_symtrid_qr(.TRUE.,.TRUE.,sca,N,D,E,WORK,N,Z,ITS,INFO)
              else
                 ! inline
                 !call d_symtrid_qr(.FALSE.,.FALSE.,sca,N,D,E,WORK,1,Z,ITS,INFO)

                 
                 ! factorize \Phi(T) and reduction to unitary Hessenberg form
                 call d_symtrid_factor3(.FALSE.,.FALSE.,SCA,N,D,E,WORK(1:(3*(N-3))),WORK((3*N+1):(5*N)),scaleqr,1,Z,INFO)
                 
                 ! check info
                 if (INFO.NE.0) then 
                    ! print error in debug mode
                    if (DEBUG) then
                       call u_infocode_check(__FILE__,__LINE__,"d_symtrid_factor failed",INFO,INFO)
                    end if
                    INFO = 1
                    ! no eigenvalues found, no back transform necessary
                 end if
                 
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

        
        ! stop timer
        call system_clock(count=c_mid)
        WORKs = WORK
        do ij = 1,N3
           ! symtrid_qr
           WORK = WORKs
           if (M.EQ.1) then
              ! call d_orthhess_qr
              if (backward) then
                 ! compute eigenvalues
              else
                 call z_unifact_qr(.FALSE.,.FALSE.,N,WORK(1:3*N-3),WORK((3*N+1):(5*N-1)),1,Z,ITS,INFO)
                 
                 ! check info
                 if (INFO.NE.0) then 
                    ! print error in debug mode
                    if (DEBUG) then
                       call u_infocode_check(__FILE__,__LINE__,"z_unifact_qr failed",INFO,INFO)
                    end if
                    INFO = 2
                    ! since some of the eigenvalues have been found, the back transform is performed for all
                    ! return
                 end if
              end if
           end if
        end do
       
        ! stop timer
        call system_clock(count=c_mid2)
        WORKs = WORK
        do ij = 1,N3
           ! symtrid_qr
           WORK = WORKs
           if (M.EQ.1) then
              ! call d_orthhess_qr
              if (backward) then
                 ! back transformation
              else
                 ! back transformation
                 do ii=1,N
                    D(ii) = WORK(3*N+2*ii)/(1d0+WORK(3*N+2*ii-1))
                 end do
                 
                 ! reverse scaling
                 if (SCA) then
                    do ii=1,N
                       D(ii) = D(ii)*scaleqr
                    end do
                 end if
              end if
           end if
        end do
        
        ! check INFO
        if (INFO.NE.0) then
           print*,"d_symtrid_qr failed."
           print*,"INFO:",INFO
        end if
        
        ! stop timer
        call system_clock(count=c_stop)

        nrm = 0d0
        do ii=1,N
           if (abs(D(ii))>nrm) then
              nrm = abs(D(ii))
           end if
        end do

        if (problem.EQ.1) then
           ! computing forward error
           do ii=1,N
              eig(ii) = (-scale*shift)+scale*cos(ii*EISCOR_DBL_PI/(N+1d0))
           end do
           
           t1 = 0d0
           t2 = 0d0
           
           do ii=1,N
              t = abs(D(ii)-eig(1))
              INFO = 1
              do jj=2,N
                 if (abs(D(ii)-eig(jj))<t) then
                    t = abs(D(ii)-eig(jj))
                    INFO = jj
                 end if
              end do
              
              if (t.GT.t2) then
                 t2 = t
                 kk = INFO
              end if
              t1 = t1 + t**2
           end do
           t1 = t1/nrm
           !if ((dsqrt(t1)>1d-3*max(scale,1d0/scale)**2).OR.(t1.NE.t1)) then
           !   do ii=1,N
           !      print*, ii, D(ii),eig(ii)
           !   end do
           !   call u_test_failed(__LINE__)
           !end if
        else 
           t1 = 0d0
           t1 = 0d0/t1
        end if
        
!        do ii=1,N
!           print*, ii, D(ii)
!        end do

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
                 !if (ii<20) then
                 ! print*, ii, D(ii)!, t3
                 !end if
                 if (t3.GT.t2) then
                    !if (ii>=20) then
                    !   print*, ii, D(ii), t3
                    !end if
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

        t2 = t2/nrm

!!$        if (N>110) then
!!$           ITS = 0
!!$           !<-5 ii = 1
!!$           !>5
!!$           do jj=1,N
!!$              INFO = int(floor((D(jj)+5d0)*10d0))
!!$              ITS(INFO) = ITS(INFO) + 1
!!$           end do
!!$           do ii=1,105
!!$              print*, ii, ITS(ii)
!!$           end do
!!$        end if

        if (backward) then
           print*, "(",N, ",",t2,")% for ", dsqrt(t1), " time ", dble(c_stop-c_start)/dble(c_rate)/N3 
        else
           print*, "(",N, ",",dble(c_stop-c_start)/dble(c_rate)/N3,")% for err ", dsqrt(t1)
           print*, "(",N, ",",dble(c_mid-c_start)/dble(c_rate)/N3,")% reduction"
           print*, "(",N, ",",dble(c_mid2-c_mid)/dble(c_rate)/N3,")% qr"
           print*, "(",N, ",",dble(c_stop-c_mid2)/dble(c_rate)/N3,")% backtransform"
        end if

        N = 2*N
     end do
  end do
  end do
  
  scale = 1d1*scale
  end do
  ! stop timer
  call system_clock(count=c_stop2)
  
  ! print success
  call u_test_passed(dble(c_stop2-c_start2)/dble(c_rate))

end program example_d_symtrid_qr_race_split
