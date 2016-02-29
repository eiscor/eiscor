!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! DAVW2
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file has been first published as part of 
! AVW --  Polynomial root finders 
! based on fast companion matrix eigensolvers on:
! https://github.com/jaurentz/avw
! It has been modified (line 423).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! The MIT License (MIT)
! 
! Copyright (c) 2015 Jared Aurentz, Raf Vandebril and David S. Watkins
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ********************************************** 
!  Last updated: December 10, 2012
! ********************************************** 
! 
! This subroutine computes the roots of a real polynomial of degree N
! in any recursively defined polynomial basis. The general 
! algorithm is outlined in: 
!
! Fast Computation of Eigenvalues of Companion, Comrade and Related Matrices, 
! BIT Numerical Mathematics 2012
!
! **************************************************************
!
! Input variables:
!
! BANDSWITCH   0 indicates no COEFFS are input, 1 indicates COEFFS are given
!
! N            degree of polynomial
!
! K            maximum length of the recursion
!
! POLY         array of dimension N+1 containing the coefficients 
!              of the polynomial in terms of the basis elements
!
! COEFFS       array of dimension (N x K) containing the recursion
!	       coefficients of the basis elements, ignored if BANDSWITCH = 0 
!
! Output variables:
!
! RROOTS       contains the real part of the computed roots, dimension N
!
! IROOTS       contains the imaginary part of the computed roots, dimension N
!
! ITERATIONS   integer array storing the number of iterations per root calculation
!
! FLAG         reports the position of any instance where the algorithm 
!              tried to divide by 0
!
! ***************************************************************


subroutine DAVW2(BANDSWITCH,N,K,POLY,COEFFS,RROOTS,IROOTS,ITERATIONS,FLAG)


   
  implicit none 

  integer, intent(in) :: BANDSWITCH,N,K
  integer, intent(inout) :: ITERATIONS(N),FLAG
  double precision, intent(in) :: POLY(N+1),COEFFS(N,K) 
  double precision, intent(inout) :: RROOTS(N),IROOTS(N)

  FLAG = 0

  if(BANDSWITCH == 0)then
    
    call DLR(N,POLY,RROOTS,IROOTS,ITERATIONS,FLAG)

  else if(BANDSWITCH == 1)then

    call DBLR(N,K,POLY,COEFFS,RROOTS,IROOTS,ITERATIONS,FLAG,0,0)

  else

    write(*,*) "Not a valid argument for BANDSWITCH"
    write(*,*) ""
    return
  end if

end subroutine

! ************************************************************** 
! 
! This subroutine calculates the roots using the companion matrix
!
! **************************************************************

subroutine DLR(N,POLY,RROOTS,IROOTS,ITERATIONS,FLAG)

   
  implicit none 

  integer, intent(in) :: N
  double precision, intent(in) :: poly(N) 
  double precision, intent(out) :: rroots(N),iroots(N)
  integer, intent(out) :: iterations(N),flag
  integer :: stop_index,start_index,zero_index,ii,jj,kk,mm,it_count,it_max,info
  double precision :: tolerance=1d-16,zero_point,diagonals,diag1,diag2,norm
  double precision :: trace,detm
  double precision :: A(2,2),L1,L2,L3,K,G(3),Blocks(2,2,N)

! Initialize and normalize blocks
  if(dabs(poly(1)) == 0)then
  write(*,*) "Polynomial is of degree less than N!"
  return
  end if

  do ii=2,N
  Blocks(1,1,ii) = -poly(ii)/poly(1)
  Blocks(1,2,ii) = 1d0
  Blocks(2,1,ii) = 1d0
  Blocks(2,2,ii) = 0d0
  end do

  Blocks(1,2,N) = -poly(N+1)/poly(1)

  Blocks(1,1,1) = 0d0
  Blocks(1,2,1) = 0d0
  Blocks(2,1,1) = 0d0
  Blocks(2,2,1) = 1d0

  
! Main loop

  stop_index = N
  start_index = 2
  zero_index = 1
  it_max = 20*N
  it_count = 0
  flag = 0

  do mm=1,it_max

  ! Zero check
  do ii=1,stop_index
    zero_point = dabs(Blocks(2,1,stop_index+1-ii))

    if(ii>1)then
      diag2 = dabs(Blocks(1,1,stop_index+2-ii)*Blocks(2,2,stop_index+1-ii))
    else
      diag2 = dabs(Blocks(2,2,stop_index+1-ii))
    end if

    if(ii<stop_index)then
      diag1 = dabs(Blocks(1,1,stop_index+1-ii)*Blocks(2,2,stop_index-ii))
    else
      diag1 = dabs(Blocks(1,1,stop_index+1-ii))
    end if

    diag1 = max(dabs(Blocks(2,1,stop_index-ii)),diag1)

    diag2 = max(dabs(Blocks(2,1,stop_index+2-ii)),diag2)

    diagonals = diag1 + diag2

    if(zero_point < tolerance*diagonals)then
      zero_index = stop_index+1-ii
      Blocks(2,1,zero_index) = 0d0

      ! Deflate down if not at bottom
      if(zero_index<stop_index)then
        Blocks(1,1,zero_index+1) = Blocks(1,1,zero_index+1)*Blocks(2,2,zero_index)
        Blocks(1,2,zero_index+1) = Blocks(1,2,zero_index+1)*Blocks(2,2,zero_index)
        Blocks(2,2,zero_index) = 1d0
        Blocks(2,1,zero_index) = 0d0
        start_index = zero_index+1
      end if
      
      exit

    end if

  end do

  ! zero at bottom
  if(zero_index==stop_index)then
    rroots(zero_index) = Blocks(2,2,zero_index)
    iroots(zero_index) = 0d0
    
    iterations(N-zero_index+1) = it_count
    it_count = 0

    if(stop_index<=1)then
      exit
    end if

    ! deflate up
    Blocks(1,2,zero_index-1) = Blocks(1,2,zero_index-1)*Blocks(1,1,zero_index)
    Blocks(2,2,zero_index-1) = Blocks(2,2,zero_index-1)*Blocks(1,1,zero_index)
    Blocks(1,1,zero_index) = 1d0
    Blocks(2,1,zero_index) = 0d0

    stop_index = stop_index-1

  ! 2x2 case
  else if(stop_index==start_index)then

    ! calculate eigenvalues using modified quadratic
    trace = Blocks(1,1,stop_index) + Blocks(2,2,stop_index)
    detm = Blocks(1,1,stop_index)*Blocks(2,2,stop_index) - Blocks(1,2,stop_index)*Blocks(2,1,stop_index)

    if(dabs(trace)**2 < 4d0*detm)then
      rroots(stop_index) = trace/2d0
      iroots(stop_index) = dsqrt(4d0*detm - dabs(trace)**2)/2d0
      rroots(stop_index-1) = trace/2d0
      iroots(stop_index-1) = -dsqrt(4d0*detm - dabs(trace)**2)/2d0
    else
      rroots(stop_index) = (trace + dsqrt(dabs(trace)**2 - 4d0*detm))/2d0
      iroots(stop_index) = 0d0
      rroots(stop_index-1) = (trace - dsqrt(dabs(trace)**2 - 4d0*detm))/2d0
      iroots(stop_index-1) = 0d0
    end if

    iterations(N-stop_index+1) = it_count
    iterations(N-stop_index+2) = 0
    it_count = 0

    stop_index = stop_index-2

    ! Breakout when done
    if(stop_index==0)then
      exit
    end if

    ! Deflate up
    Blocks(1,2,stop_index) = Blocks(1,2,stop_index)*Blocks(1,1,stop_index+1)
    Blocks(2,2,stop_index) = Blocks(2,2,stop_index)*Blocks(1,1,stop_index+1)

    Blocks(1,1,stop_index+1) = 1d0
  
  ! more than one block
  else
    it_count = it_count + 1

    ! Build G
    if(mm == 0)then
      call dnormalpoly(2,G(1:2))
      trace = G(1) + G(2)
      detm = G(1)*G(2)
    else if(mod(it_count,15) == 0)then
      call dnormalpoly(2,G(1:2))
      G(2) = G(1)*G(1) + G(2)*G(2)
      G(2) = dsqrt(G(2))
      trace = 2d0*G(1)/G(2)
      detm = 1d0
      write(*,*) "Exceptional shift in DLR!"
    else
      trace = Blocks(2,2,stop_index-1)*Blocks(1,1,stop_index) + Blocks(2,2,stop_index)
      detm = Blocks(1,1,stop_index)*Blocks(2,2,stop_index) - Blocks(1,2,stop_index)*Blocks(2,1,stop_index)
      detm = Blocks(2,2,stop_index-1)*detm
    end if



    G(1) = Blocks(1,1,start_index)*Blocks(1,1,start_index) - trace*Blocks(1,1,start_index) + detm
    G(1) = G(1) + Blocks(1,2,start_index)*Blocks(2,1,start_index)*Blocks(1,1,start_index+1)
    G(2) = Blocks(1,1,start_index) + Blocks(2,2,start_index)*Blocks(1,1,start_index+1) - trace
    G(2) = Blocks(2,1,start_index)*G(2)
    G(3) = Blocks(2,1,start_index)*Blocks(2,1,start_index+1)

    ! Build L3 L2 and L1
    L2 = G(3)/G(2);
    L3 = G(2)/G(1);
    
    ! Build bulge
    A = Blocks(:,:,start_index)
    Blocks(2,:,start_index) = Blocks(2,:,start_index) - L3*Blocks(1,:,start_index)
  
    L1 = -L2*A(2,1)/Blocks(2,1,start_index)
    K = -L2*A(2,2) - L1*Blocks(2,2,start_index) 

    Blocks(2,:,start_index+1) = Blocks(2,:,start_index+1) + K*Blocks(1,:,start_index+1)

    ! Bulgechase
    do kk=(start_index+1),(stop_index)

      ! Exit if at bottom
      if(kk==stop_index)then
  
        Blocks(:,1,kk) = Blocks(:,1,kk) + L2*Blocks(:,2,kk)

        call DGATO(Blocks(:,:,kk-1),Blocks(:,:,kk),L3,flag)

        Blocks(:,1,kk) = Blocks(:,1,kk) + L1*Blocks(:,2,kk)
        Blocks(:,1,kk) = Blocks(:,1,kk) + L3*Blocks(:,2,kk)

        exit
      end if  

      call DGATO(Blocks(:,:,kk),Blocks(:,:,kk+1),L2,flag)
      call DGATO(Blocks(:,:,kk-1),Blocks(:,:,kk),L3,flag)

      call DGLTO(L1,L2,L3,flag)  
  
    end do

  end if

  end do

  ! set flag for maxing out loop
  if(mm == it_max)then
    flag = 1
  end if

end subroutine 

! ************************************************************** 
! 
! This subroutine calculates the roots using any congenial matrix
!
! **************************************************************

subroutine DBLR(N,K,POLY,COEFFS,RROOTS,IROOTS,ITERATIONS,FLAG,INIT,CHASE)

  implicit none 

  integer, intent(in) :: N,K,INIT,CHASE
  double precision, intent(in) :: POLY(N+1),COEFFS(N,K)
  double precision, intent(inout) :: RROOTS(N),IROOTS(N)
  integer, intent(out) :: ITERATIONS(N),flag

  integer :: stop_index,start_index,zero_index,ii,jj,hh,mm,it_count=0,it_max
  integer :: length,switch=0
  double precision :: tolerance=1d-16,zero_point,diagonals,diag1,diag2
  double precision :: L1,L2,L3,W(2,2),A(3,2),M,bulge,trace,detm,G(3)
  double precision, allocatable :: spike(:),Blocks(:,:,:),Bands(:,:) 


! Initialize Bands and Blocks
  
  allocate(spike(N),Blocks(2,2,N),Bands(N,K))

  if(dabs(poly(1)) == 0 .OR. dabs(Coeffs(1,1)) == 0)then
    write(*,*) "Polynomial is of degree less than N!"
    flag = 1
    return
  end if

  spike = -Poly(2:(N+1))/Poly(1)*Coeffs(1,1)
  spike(1:(K-1)) = spike(1:(K-1)) + Coeffs(1,2:K)

  Bands = 0d0
  Bands(1:(N-1),1:K) = Coeffs(2:N,1:K)

  Blocks(:,:,1) = 0d0
  Blocks(2,2,1) = 1d0

  do ii=1,(N-1)
    if(abs(Bands(ii,1)) == 0)then
      write(*,*) "Matrix is not properly upper-hessenberg!"
      flag = 2
      return
    end if

    L1 = spike(ii)/Bands(ii,1)

    length = min(N+1-ii,K)
    spike(ii:(ii+length-1)) = spike(ii:(ii+length-1)) - L1*Bands(ii,1:length)

    Blocks(1,1,ii+1) = L1
    Blocks(1,2,ii+1) = 1d0  
    Blocks(2,1,ii+1) = 1d0
    Blocks(2,2,ii+1) = 0d0
  end do

  Bands(N,1) = spike(N)

! Main loop

  stop_index = N
  start_index = 2
  zero_index = 1
  it_max = 20*N
  it_count = 0
  flag = 0

  do mm=1,it_max

  ! Zero check
  do ii=1,stop_index
    zero_point = dabs(Blocks(2,1,stop_index+1-ii))

    if(ii>1)then
      diag2 = dabs(Blocks(1,1,stop_index+2-ii)*Blocks(2,2,stop_index+1-ii)*Bands(stop_index+1-ii,1))
    else
      diag2 = dabs(Blocks(2,2,stop_index+1-ii)*Bands(stop_index+1-ii,1))
    end if

    if(ii<stop_index)then
      diag1 = dabs(Blocks(1,1,stop_index+1-ii)*Blocks(2,2,stop_index-ii)*Bands(stop_index-ii,1))
    else
       !print*,stop_index
       !print*,ii
       !diag1 = dabs(Blocks(1,1,stop_index+1-ii)*Bands(stop_index-ii,1))
       diag1 = 0d0
    end if

    diagonals = diag1 + diag2

    if(zero_point < tolerance*diagonals)then
      zero_index = stop_index+1-ii
      Blocks(2,1,zero_index) = 0d0

      ! Deflate down if not at bottom
      if(zero_index<stop_index)then
        Blocks(1,1,zero_index+1) = Blocks(1,1,zero_index+1)*Blocks(2,2,zero_index)
        Blocks(1,2,zero_index+1) = Blocks(1,2,zero_index+1)*Blocks(2,2,zero_index)
        Blocks(1,2,zero_index) = Blocks(1,2,zero_index)/Blocks(2,2,zero_index)
        Blocks(2,2,zero_index) = 1d0
        start_index = zero_index+1
      end if

      exit

    end if
  end do

  ! zero at bottom
  if(zero_index==stop_index)then
    rroots(zero_index) = Blocks(2,2,zero_index)*Bands(zero_index,1)
    iroots(zero_index) = 0d0
    
    iterations(N-zero_index+1) = it_count
    it_count = 0

    if(stop_index<=1)then
      exit
    end if

    ! deflate up
    Blocks(:,2,zero_index-1) = Blocks(1,1,zero_index)*Blocks(:,2,zero_index-1)
    Blocks(1,1,zero_index) = 1d0

    stop_index = stop_index-1

  ! 2x2 case
  else if(stop_index==start_index)then
    ! calculate eigenvalues using modified quadratic
    if(K == 1)then
    W(1,1) = Blocks(1,1,stop_index)*Bands(stop_index-1,1)
    W(1,2) = Blocks(1,2,stop_index)*Bands(stop_index,1)
    W(2,1) = Blocks(2,1,stop_index)*Bands(stop_index-1,1)
    W(2,2) = Blocks(2,2,stop_index)*Bands(stop_index,1)
    else
    W(1,1) = Blocks(1,1,stop_index)*Bands(stop_index-1,1)
    W(1,2) = Blocks(1,1,stop_index)*Bands(stop_index-1,2) + Blocks(1,2,stop_index)*Bands(stop_index,1)
    W(2,1) = Blocks(2,1,stop_index)*Bands(stop_index-1,1)
    W(2,2) = Blocks(2,1,stop_index)*Bands(stop_index-1,2) + Blocks(2,2,stop_index)*Bands(stop_index,1)
    end if

    trace = W(1,1) + W(2,2)
    detm = W(1,1)*W(2,2) - W(1,2)*W(2,1)

    if(dabs(trace)**2 < 4d0*detm)then
      rroots(stop_index) = trace/2d0
      iroots(stop_index) = dsqrt(4d0*detm - dabs(trace)**2)/2d0
      rroots(stop_index-1) = trace/2d0
      iroots(stop_index-1) = -dsqrt(4d0*detm - dabs(trace)**2)/2d0
    else
      rroots(stop_index) = (trace + dsqrt(dabs(trace)**2 - 4d0*detm))/2d0
      iroots(stop_index) = 0d0
      rroots(stop_index-1) = (trace - dsqrt(dabs(trace)**2 - 4d0*detm))/2d0
      iroots(stop_index-1) = 0d0
    end if

    iterations(N-stop_index+1) = it_count
    iterations(N-stop_index+2) = 0
    it_count = 0

    stop_index = stop_index-2

    ! Breakout when done
    if(stop_index==0)then
      exit
    end if

    ! Deflate up
    Blocks(1,2,stop_index) = Blocks(1,2,stop_index)*Blocks(1,1,stop_index+1)
    Blocks(2,2,stop_index) = Blocks(2,2,stop_index)*Blocks(1,1,stop_index+1)

    Blocks(1,1,stop_index+1) = 1d0

  ! more than one block
  else
    it_count = it_count + 1

    ! shift calc
    if(K == 1)then
    W(1,1) = Blocks(1,1,stop_index)*Bands(stop_index-1,1)*Blocks(2,2,stop_index-1)
    W(1,2) = Blocks(1,2,stop_index)*Bands(stop_index,1)
    W(1,2) = W(1,2)*Blocks(2,2,stop_index-1)
    W(2,1) = Blocks(2,1,stop_index)*Bands(stop_index-1,1)
    W(2,2) = Blocks(2,2,stop_index)*Bands(stop_index,1)
    else
    W(1,1) = Blocks(1,1,stop_index)*Bands(stop_index-1,1)*Blocks(2,2,stop_index-1)
    W(1,2) = Blocks(1,1,stop_index)*Bands(stop_index-1,2) + Blocks(1,2,stop_index)*Bands(stop_index,1)
    W(1,2) = W(1,2)*Blocks(2,2,stop_index-1)
    W(2,1) = Blocks(2,1,stop_index)*Bands(stop_index-1,1)
    W(2,2) = Blocks(2,1,stop_index)*Bands(stop_index-1,2) + Blocks(2,2,stop_index)*Bands(stop_index,1)
    end if

    if(mm == 0)then
      trace = 0d0
      detm = -1d0
    else if(mod(it_count,15) == 0)then
      call dnormalpoly(2,G(1:2))
      G(2) = G(1)*G(1) + G(2)*G(2)
      G(2) = dsqrt(G(2))
      trace = 2d0*G(1)/G(2)
      detm = 1d0
      write(*,*) "Exceptional shift in DBLR!"
    else
      trace = W(1,1) + W(2,2)
      detm = W(1,1)*W(2,2) - W(1,2)*W(2,1)
    end if

    ! Build G
    A = 0d0
    if(K == 1)then
      A(1,1) = Bands(start_index-1,1)
      A(2,2) = Bands(start_index,1)
    else
      A(1,1:2) = Bands(start_index-1,1:2)
      A(2,2) = Bands(start_index,1)
    end if

    A(2:3,:) = matmul(Blocks(:,:,start_index+1),A(2:3,:))
    A(1:2,:) = matmul(Blocks(:,:,start_index),A(1:2,:))

    G(1) = A(1,1)**2 - trace*A(1,1) + detm + A(1,2)*A(2,1)
    G(2) = A(2,1)*(A(1,1) + A(2,2) - trace)
    G(3) = A(2,1)*A(3,2)

    ! Build L1 L2 and L3
    L2 = G(3)/G(2);
    L3 = G(2)/G(1);
    
    ! Build bulge
    call DBB(L1,L2,L3,start_index,N,Blocks,init,flag)

    ! Move L2 through bands
    call DBP(L2,start_index+1,start_index,stop_index,N,K,Bands)

    ! Move L3 through bands
    call DBP(L3,start_index,start_index,stop_index,N,K,Bands)

    ! Turnover iterations
    if(chase == 0)then
    do hh=(start_index+2),(stop_index)
      ! Pass through Blocks and turnover L's over
      call DGATO(Blocks(:,:,hh-1),Blocks(:,:,hh),L2,flag)
      call DGATO(Blocks(:,:,hh-2),Blocks(:,:,hh-1),L3,flag)

      call DGLTO(L1,L2,L3,flag)

      ! Move L2 through bands
      call DBP(L2,hh,start_index,stop_index,N,K,Bands)

      ! Move L3 through bands
      call DBP(L3,hh-1,start_index,stop_index,N,K,Bands)

    end do

    ! Dump the Bulge      
    Blocks(:,1,stop_index) = Blocks(:,1,stop_index) + L2*Blocks(:,2,stop_index)
    
    call DGATO(Blocks(:,:,stop_index-1),Blocks(:,:,stop_index),L3,flag)

    ! Move L1 through bands
    call DBP(L1,stop_index,start_index,stop_index,N,K,Bands)

    Blocks(:,1,stop_index) = Blocks(:,1,stop_index) + L1*Blocks(:,2,stop_index)

    ! Move L3 through bands
    call DBP(L3,stop_index,start_index,stop_index,N,K,Bands)

    Blocks(:,1,stop_index) = Blocks(:,1,stop_index) + L3*Blocks(:,2,stop_index)

    else if(chase == 1)then
    do hh=(start_index+2),(stop_index)
      ! Pass through Blocks and turnover L's over
      call DGATO(Blocks(:,:,hh-1),Blocks(:,:,hh),L2,flag)
      call DGATO(Blocks(:,:,hh-2),Blocks(:,:,hh-1),L3,flag)

      ! Move L1 through bands
      call DBP(L1,hh-1,start_index,stop_index,N,K,Bands)

      ! Move L2 through bands
      call DBP(L2,hh,start_index,stop_index,N,K,Bands)

      ! Move L3 through bands
      call DBP(L3,hh-1,start_index,stop_index,N,K,Bands)

      ! Pass L1 through Blocks
      call DGATO(Blocks(:,:,hh-1),Blocks(:,:,hh),L1,flag)

    end do

    ! Dump the Bulge      
    Blocks(:,1,stop_index) = Blocks(:,1,stop_index) + L2*Blocks(:,2,stop_index)
    
    call DGATO(Blocks(:,:,stop_index-1),Blocks(:,:,stop_index),L3,flag)

    ! Move L1 through bands
    call DBP(L1,stop_index,start_index,stop_index,N,K,Bands)

    Blocks(:,1,stop_index) = Blocks(:,1,stop_index) + L1*Blocks(:,2,stop_index)

    ! Move L3 through bands
    call DBP(L3,stop_index,start_index,stop_index,N,K,Bands)

    Blocks(:,1,stop_index) = Blocks(:,1,stop_index) + L3*Blocks(:,2,stop_index)

    end if
  end if

  end do

  ! set flag for maxing out loop
  if(mm == it_max)then
    flag = 1
  end if

  ! free memory
  deallocate(spike,Blocks,Bands)

end subroutine

! **********************************************************************
!
!
!
! **********************************************************************

subroutine DGATO(A,B,L,flag)

implicit none

integer :: flag
double precision :: A(2,2),B(2,2),L

  A(:,1) = A(:,1) + L*B(1,1)*A(:,2)

  if(dabs(A(2,1)) == 0)then
    flag = 1
    return
  end if

  L = L*B(2,1)/A(2,1)

  B(2,:) = B(2,:) - L*A(2,2)*B(1,:)

end subroutine

! **********************************************************************
!
!
!
! **********************************************************************

subroutine DGLTO(L1,L2,L3,flag)

implicit none

integer :: flag
double precision :: L1,L2,L3,temp1,temp2,temp3

  temp1 = L1
  temp2 = L2
  temp3 = L3

  L3 = temp1 + temp3

  if(dabs(L3) == 0)then
    write(*,*) "Tried to divide by 0. Algorithm failed!"
    flag = 1
    return
  end if

  L2 = temp2*temp3/L3
  L1 = temp2 - L2

end subroutine

! **********************************************************************
!
! Double Precision Band Pass
!
! **********************************************************************

subroutine DBP(L,start,top,bottom,N,K,Bands)

implicit none

integer :: start,top,bottom,N,K,jj,length
double precision :: L,Bands(N,K),bulge

    length = min(start-top+1,K-1)      
    do jj=1,length
      Bands(start-jj,jj) = Bands(start-jj,jj) + L*Bands(start-jj,jj+1)
    end do
      
    bulge = L*Bands(start,1)

    L = bulge/Bands(start-1,1)

    length = min(K-1,bottom+1-start)
    do jj=1,length
      Bands(start,jj) = Bands(start,jj) - L*Bands(start-1,jj+1)
    end do   

end subroutine

! **********************************************************************
!
! Double Precision Bulge Build
!
! **********************************************************************

subroutine DBB(L1,L2,L3,start,N,Blocks,switch,flag)

implicit none

integer :: start,N,switch,flag
double precision :: L1,L2,L3,Blocks(2,2,N),M,W(2,2)

    if(switch == 0)then
      W = Blocks(:,:,start)
      Blocks(2,:,start) = Blocks(2,:,start) - L3*Blocks(1,:,start)

      if(dabs(Blocks(2,1,start)) == 0)then
        flag = 1
        return
      end if
  
      L1 = -L2*W(2,1)/Blocks(2,1,start)
      M = -L2*W(2,2) - L1*Blocks(2,2,start) 

      Blocks(2,:,start+1) = Blocks(2,:,start+1) + M*Blocks(1,:,start+1)

    else if(switch == 1)then
      Blocks(2,:,start+1) = Blocks(2,:,start+1) - L2*Blocks(2,2,start)*Blocks(1,:,start+1)

      if(dabs(Blocks(2,1,start+1)) == 0)then
        flag = 1
        return
      end if
    
      L1 = -L2*Blocks(2,1,start)/Blocks(2,1,start+1)

      Blocks(:,1,start)= Blocks(:,1,start) - L1*Blocks(1,1,start+1)*Blocks(:,2,start)

      Blocks(2,:,start) = Blocks(2,:,start) - L3*Blocks(1,:,start)

      call DGATO(Blocks(:,:,start),Blocks(:,:,start+1),L1,flag)

    end if

end subroutine

! **********************************************************************
!
!
!
! **********************************************************************

subroutine dnormalpoly(N,poly)

  implicit none
  
  integer, intent(in) :: N
  real(kind(1d0)), intent(inout) :: poly(N) 
  
  double precision :: u,v,s,pi = 3.141592653589793239d0
  integer :: ii,jj

  poly = 0d0

  do ii=1,N
    do jj=1,20
      
      call random_number(u)
      call random_number(v)
  
      s = u**2 + v**2
  
      if(s > 0 .and. s < 1)then        
        poly(ii) = dcos(2d0*pi*v)*dsqrt(-2d0*log(u))
        exit
      end if
    end do

  end do


end subroutine

! **********************************************************************
!
!
!
! **********************************************************************

SUBROUTINE init_random_seed()

  implicit none

        INTEGER :: ii, kk, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
        CALL RANDOM_SEED(size = kk)
        ALLOCATE(seed(kk))
          
        CALL SYSTEM_CLOCK(COUNT=clock)
          
        seed = clock + 37 * (/ (ii - 1, ii = 1, kk) /)
        CALL RANDOM_SEED(PUT = seed)
          
  DEALLOCATE(seed)
END SUBROUTINE
