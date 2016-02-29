!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_singleshift
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file has been first published as part of 
! AVW --  Polynomial root finders 
! based on fast companion matrix eigensolvers on:
! https://github.com/jaurentz/avw
! It has been modified.
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
! 
! This subroutine computes three residuals for each computed root, 
! lambda, of a polynomial P(x). It is also capable of applying an 
! arbitrary number of Newton iterations to each computed root. 
!
! The residuals are |P(lambda)/P'(lambda)|, 
! |P(lambda)/P'(lambda)/lambda|, and ||Cv-lambda v||/||C||/||v||, 
! in the inifinity norm where C is the companion matrix and v is 
! the eigenvectro associated with lambda.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Input variables:
!
! N            degree of polynomial 
!
! K            Bandwidth of upper-triangular matrix 
!
! POLY         complex array of length N, containing the nonleading 
!              coefficients of P (ordered with decreasing N). P must
!              be monic
!
! COEFFS       complex array of dimension (N,K) containing the recursion
!              coefficients for the polynomial basis, COEFFS(1,1) must
!              be 1
!
! ROOTS        complex array of length N, containing the computed roots
!              of POLY
!
! NEWTNUM      a non-negative integer specifying the number of Newton
!              iterations zero is acceptable.
!
! Output variables:
!
! ALLROOTS     complex array of dimension (N,NEWTNUM+1) contains the 
!              original roots as well as any newton corrections
!
! RESIDUALS    double precision array of dimension (N,3*(NEWTNUM+1)).
!              each row of RESIDUALS corresponds to one root of POLY.
!              columns are in sets of three, with the columns 1, 2 
!              and 3 corresponding to the three residuals mentioned 
!              above respectively. Every set of three columns 
!              corresponds to a Newton iteration except for the 
!              first one.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_polyc_residuals(N,K,NEWTNUM,POLY,COEFFS,ROOTS,ALLROOTS,RESIDUALS)

  implicit none
  
  integer, intent(in) :: N,K,newtnum
  complex(8), intent(in) :: poly(N),roots(N),coeffs(N,K)
  real(8), intent(inout) :: residuals(N,3*(newtnum+1))
  complex(8), intent(inout) :: allroots(N,newtnum+1)
  
  integer ii,jj,kk,ll,length
  real(8) :: Cnorm,temp,Pnorm
  complex(8) :: f, fprime, lambda
  complex(8), allocatable :: P(:),Pprime(:)
  
  ! allocate memory
  allocate(P(N+1),Pprime(N+1))
  
  ! initialize residuals
  residuals = 10d0
  
  ! initialize allroots
  allroots = cmplx(0d0,0d0,kind=8)
  allroots(:,1) = roots
  
  ! Matrix infinity norm
  P(1:N) = poly
  P(1:(K-1)) = P(1:(K-1))- Coeffs(1,2:K)
  
  Cnorm = 0d0
  do ii=1,N
     Cnorm = Cnorm + abs(P(ii))
  end do
  
  do ii=2,N
     temp = 0d0
     length = min(K,N+2-ii)
     do jj=1,length
        temp = temp + abs(Coeffs(ii,jj))
     end do
     
     if(temp > Cnorm)then
        Cnorm = temp
     end if
  end do
  
  ! Function evaluations and Newton Corrections
  do ii=1,N    
     do jj=1,(newtnum+1)
        ! Roots inside or on the unit circle
        if(abs(allroots(ii,jj)) <= 1d0)then
           ! function evals
           lambda = allroots(ii,jj)
           
           ! P_k(lambda)
           P(1) = cmplx(1d0,0d0,kind=8)
           
           if(K == 1)then
              P(2) = lambda/Coeffs(N,1)
           else
              P(2) = (lambda - Coeffs(N,2))/Coeffs(N,1)
           end if
           
           do kk=1,(N-1)
              if(K == 1)then
                 P(kk+2) = lambda*P(kk+1)/Coeffs(N-kk,1)
              else
                 P(kk+2) = (lambda - Coeffs(N-kk,2))*P(kk+1)
                 
                 length=min(kk,K-2)
                 do ll=1,length
                    P(kk+2) = P(kk+2) - Coeffs(N-kk,2+ll)*P(kk+1-ll)
                 end do
                 
                 P(kk+2) = P(kk+2)/Coeffs(N-kk,1)
              end if
           end do                  
           
           ! P'_k(lambda)
           Pprime(1) = cmplx(0d0,0d0,kind=8)
           
           Pprime(2) = cmplx(1d0,0d0,kind=8)/Coeffs(N,1)
           
           do kk=1,(N-1)
              if(K == 1)then
                 Pprime(kk+2) = (lambda*Pprime(kk+1) + P(kk+1))/Coeffs(N-kk,1)
              else
                 Pprime(kk+2) = (lambda - Coeffs(N-kk,2))*Pprime(kk+1) + P(kk+1)
                 length=min(kk,K-2)
                 do ll=1,length
                    Pprime(kk+2) = Pprime(kk+2) - Coeffs(N-kk,2+ll)*Pprime(kk+1-ll)
                 end do
                 
                 Pprime(kk+2) = Pprime(kk+2)/Coeffs(N-kk,1)
                 
              end if
           end do
           
           ! compute vector norm
           Pnorm = abs(P(1))
           do kk=2,N+1
              temp = abs(P(kk))
              if(temp > Pnorm)then
                 Pnorm = temp
              end if
           end do
           
           !P(lambda) and P'(lambda)
           f = P(N+1)
           fprime = Pprime(N+1)
           do kk=1,N
              f = f + poly(kk)*P(N+1-kk)
              fprime = fprime + poly(kk)*Pprime(N+1-kk)
           end do
           
           ! Store residuals
           residuals(ii,3*(jj-1)+1) = abs(f/fprime)
           residuals(ii,3*(jj-1)+2) = abs(f/fprime/lambda)
           residuals(ii,3*(jj-1)+3) = abs(f)/Cnorm/Pnorm
           
           ! Newton correction
           if((newtnum+1-jj) > 0)then
              lambda = lambda - f/fprime
              allroots(ii,jj+1) = lambda
           end if
           ! Roots outside the unit circle
        else
           ! function evals
           lambda = allroots(ii,jj)
           
           ! P_k(lambda)
           P(1) = cmplx(1d0,0d0,kind=8)/lambda
           Pprime(1) = cmplx(0d0,0d0,kind=8)
           
           if(K == 1)then
              P(2) = P(1)/Coeffs(N,1)
              Pprime(2) = cmplx(1d0,0d0,kind=8)/Coeffs(N,1)/lambda/lambda
           else
              P(2) = (lambda - Coeffs(N,2))*P(1)/Coeffs(N,1)/lambda
              Pprime(2) = cmplx(1d0,0d0,kind=8)/Coeffs(N,1)/lambda/lambda
           end if
           
           P(1) = P(1)/lambda
           
           do kk=1,(N-1)
              if(K == 1)then
                 P(kk+2) = P(kk+1)/Coeffs(N-kk,1)
                 Pprime(kk+2) = (lambda*Pprime(kk+1) + P(kk+1))/Coeffs(N-kk,1)
              else
                 P(kk+2) = (lambda - Coeffs(N-kk,2))*P(kk+1)/lambda
                 Pprime(kk+2) = (lambda - Coeffs(N-kk,2))*Pprime(kk+1) + P(kk+1)
                 
                 length=min(kk,K-2)
                 do ll=1,length
                    P(kk+2) = P(kk+2) - Coeffs(N-kk,2+ll)*P(kk+1-ll)/lambda
                    Pprime(kk+2) = Pprime(kk+2) - Coeffs(N-kk,2+ll)*Pprime(kk+1-ll)
                 end do
                 
                 P(kk+2) = P(kk+2)/Coeffs(N-kk,1)
                 Pprime(kk+2) = Pprime(kk+2)/Coeffs(N-kk,1)
                 
              end if
              
              P(1:kk+1) = P(1:kk+1)/lambda
              Pprime(2:kk+2) = Pprime(2:kk+2)/lambda
              
           end do
           
           ! compute vector norm
           Pnorm = abs(P(1))
           do kk=2,N+1
              temp = abs(P(kk))
              if(temp > Pnorm)then
                 Pnorm = temp
              end if
           end do
           
           !P(lambda) and P'(lambda)
           f = P(N+1)
           fprime = Pprime(N+1)
           
           do kk=1,N
              f = f + poly(kk)*P(N+1-kk)
              fprime = fprime + poly(kk)*Pprime(N+1-kk)
           end do
           
           ! Store residuals
           residuals(ii,3*(jj-1)+1) = abs(f/fprime)
           residuals(ii,3*(jj-1)+2) = abs(f/fprime/lambda)
           residuals(ii,3*(jj-1)+3) = abs(f)/Cnorm/Pnorm
           
           ! Newton correction
           if((newtnum+1-jj) > 0)then
              lambda = lambda - f/fprime
              allroots(ii,jj+1) = lambda
           end if
        end if
     end do
     
     if (residuals(ii,3*(newtnum)+1).NE.residuals(ii,3*(newtnum)+1)) then                 
        call zq_polyc_residual(N,K,NEWTNUM,POLY,COEFFS,ROOTS(ii),&
             &ALLROOTS(ii,:),RESIDUALS(ii,:))
     end if     
     
  end do
  
  ! free memory
  deallocate(P,Pprime)
  
  
end subroutine z_polyc_residuals

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fact_singleshift
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Same function but with higher precision for internal operations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Input variables:
!
! N            degree of polynomial 
!
! K            Bandwidth of upper-triangular matrix 
!
! POLY         complex array of length N, containing the nonleading 
!              coefficients of P (ordered with decreasing N). P must be
!              monic
!
! COEFFS       complex array of dimension (N,K) containing the recursion
!              coefficients for the polynomial basis, COEFFS(1,1) must be 1
!
!
! ROOTS        complex number, containing the computed root
!              of POLY
!
! NEWTNUM      a non-negative integer specifying the number of Newton
!              iterations zero is acceptable.
!
! Output variables:
!
! ONEROOT      complex array of dimension (NEWTNUM+1) contains the original 
!              roots as well as any newton corrections
!
! RESIDUALS    double precision array of dimension (3*(NEWTNUM+1)).
!              each row of RESIDUALS corresponds to one root of POLY.
!              columns are in sets of three, with the columns 1, 2 and 3 
!              corresponding to the three residuals mentioned above 
!              respectively. Every set of three columns corresponds to a 
!              Newton iteration except for the first one.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zq_polyc_residual(N,K,NEWTNUM,POLY,COEFFS,ROOT,ONEROOT,RESIDUAL)

  implicit none
  
  integer, intent(in) :: N,K,newtnum
  complex(8), intent(in) :: poly(N),root,coeffs(N,K)
  real(8), intent(inout) :: residual(3*(newtnum+1))
  complex(8), intent(inout) :: oneroot(newtnum+1)
  
  integer ii,jj,kk,ll,length
  real(16) :: Cnorm,temp,Pnorm
  complex(16) :: f, fprime, lambda
  complex(16), allocatable :: P(:),Pprime(:)
  
  ! allocate memory
  allocate(P(N+1),Pprime(N+1))
  
  ! initialize residual
  residual = 10d0
  
  ! initialize allroots
  oneroot = cmplx(0d0,0d0,kind=8)
  oneroot(1) = root
  
  ! Matrix infinity norm
  P(1:N) = poly
  P(1:(K-1)) = P(1:(K-1))- Coeffs(1,2:K)
  
  Cnorm = 0d0
  do ii=1,N
     Cnorm = Cnorm + abs(P(ii))
  end do
  
  do ii=2,N
     temp = 0d0
     length = min(K,N+2-ii)
     do jj=1,length
        temp = temp + abs(Coeffs(ii,jj))
     end do
     
     if(temp > Cnorm)then
        Cnorm = temp
     end if
  end do
  
  ! Function evaluations and Newton Corrections
  do jj=1,(newtnum+1)
     ! Roots inside or on the unit circle
     if(abs(oneroot(jj)) <= 1d0)then
        ! function evals
        lambda = cmplx(oneroot(jj),kind=16)
        
        ! P_k(lambda)
        P(1) = cmplx(1d0,0d0,kind=16)
        
        if(K == 1)then
           P(2) = lambda/Coeffs(N,1)
        else
           P(2) = (lambda - Coeffs(N,2))/Coeffs(N,1)
        end if
           
        do kk=1,(N-1)
           if(K == 1)then
              P(kk+2) = lambda*P(kk+1)/Coeffs(N-kk,1)
           else
              P(kk+2) = (lambda - Coeffs(N-kk,2))*P(kk+1)
              
              length=min(kk,K-2)
              do ll=1,length
                 P(kk+2) = P(kk+2) - Coeffs(N-kk,2+ll)*P(kk+1-ll)
              end do
              
              P(kk+2) = P(kk+2)/Coeffs(N-kk,1)
           end if
        end do
        
        ! P'_k(lambda)
        Pprime(1) = cmplx(0d0,0d0,kind=16)
        
        Pprime(2) = cmplx(1d0,0d0,kind=16)/Coeffs(N,1)
        
        do kk=1,(N-1)
           if(K == 1)then
              Pprime(kk+2) = (lambda*Pprime(kk+1) + P(kk+1))/Coeffs(N-kk,1)
           else
              Pprime(kk+2) = (lambda - Coeffs(N-kk,2))*Pprime(kk+1) + P(kk+1)
              length=min(kk,K-2)
              do ll=1,length
                 Pprime(kk+2) = Pprime(kk+2) - Coeffs(N-kk,2+ll)*Pprime(kk+1-ll)
              end do
              
              Pprime(kk+2) = Pprime(kk+2)/Coeffs(N-kk,1)
              
           end if
        end do
        
        ! compute vector norm
        Pnorm = abs(P(1))
        do kk=2,N+1
           temp = abs(P(kk))
           if(temp > Pnorm)then
              Pnorm = temp
           end if
        end do
        
           !P(lambda) and P'(lambda)
        f = P(N+1)
        fprime = Pprime(N+1)
        do kk=1,N
           f = f + poly(kk)*P(N+1-kk)
           fprime = fprime + poly(kk)*Pprime(N+1-kk)
        end do
        
        ! Store residual
        residual(3*(jj-1)+1) = real(abs(f/fprime),kind=8)
        residual(3*(jj-1)+2) = real(abs(f/fprime/lambda),kind=8)
        residual(3*(jj-1)+3) = real(abs(f)/Cnorm/Pnorm,kind=8)
           
        ! Newton correction
        if((newtnum+1-jj) > 0)then
           lambda = lambda - f/fprime
           oneroot(jj+1) = cmplx(lambda,kind=8)
        end if
        ! Roots outside the unit circle
     else
        ! function evals
        lambda = oneroot(jj)
        
        ! P_k(lambda)
        P(1) = cmplx(1d0,0d0,kind=16)/lambda
        Pprime(1) = cmplx(0d0,0d0,kind=16)
        
        if(K == 1)then
           P(2) = P(1)/Coeffs(N,1)
           Pprime(2) = cmplx(1d0,0d0,kind=16)/Coeffs(N,1)/lambda/lambda
        else
           P(2) = (lambda - Coeffs(N,2))*P(1)/Coeffs(N,1)/lambda
           Pprime(2) = cmplx(1d0,0d0,kind=16)/Coeffs(N,1)/lambda/lambda
        end if
        
        P(1) = P(1)/lambda
        
        do kk=1,(N-1)
           if(K == 1)then
              P(kk+2) = P(kk+1)/Coeffs(N-kk,1)
              Pprime(kk+2) = (lambda*Pprime(kk+1) + P(kk+1))/Coeffs(N-kk,1)
           else
              P(kk+2) = (lambda - Coeffs(N-kk,2))*P(kk+1)/lambda
              Pprime(kk+2) = (lambda - Coeffs(N-kk,2))*Pprime(kk+1) + P(kk+1)
              
              length=min(kk,K-2)
              do ll=1,length
                 P(kk+2) = P(kk+2) - Coeffs(N-kk,2+ll)*P(kk+1-ll)/lambda
                 Pprime(kk+2) = Pprime(kk+2) - Coeffs(N-kk,2+ll)*Pprime(kk+1-ll)
              end do
              
              P(kk+2) = P(kk+2)/Coeffs(N-kk,1)
              Pprime(kk+2) = Pprime(kk+2)/Coeffs(N-kk,1)
              
           end if
           
           P(1:kk+1) = P(1:kk+1)/lambda
           Pprime(2:kk+2) = Pprime(2:kk+2)/lambda
           
        end do
        
        ! compute vector norm
        Pnorm = abs(P(1))
        do kk=2,N+1
           temp = abs(P(kk))
           if(temp > Pnorm)then
              Pnorm = temp
           end if
        end do
        
        !P(lambda) and P'(lambda)
        f = P(N+1)
        fprime = Pprime(N+1)
        
        do kk=1,N
           f = f + poly(kk)*P(N+1-kk)
           fprime = fprime + poly(kk)*Pprime(N+1-kk)
        end do
        
        ! Store residual
        residual(3*(jj-1)+1) = real(abs(f/fprime),kind=8)
        residual(3*(jj-1)+2) = real(abs(f/fprime/lambda),kind=8)
        residual(3*(jj-1)+3) = real(abs(f)/Cnorm/Pnorm,kind=8)
           
        ! Newton correction
        if((newtnum+1-jj) > 0)then
           lambda = lambda - f/fprime
           oneroot(jj+1) = cmplx(lambda,kind=8)
        end if
     end if
  end do
  
  ! free memory
  deallocate(P,Pprime)
    
end subroutine zq_polyc_residual
