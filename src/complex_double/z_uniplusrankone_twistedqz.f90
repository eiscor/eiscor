!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uniplusrankone_twistedqz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the generalized schur decomposition of an 
! extended upper-hessenberg, upper-triangular pencil. Both the hessenberg
! and triangular matrices are the sum of a unitary matrix and a rank 
! one matrix. These matrices are stored in 5 sequences of rotations.
! 
! This factorization is described in:
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  COMPZ           CHARACTER
!                    'N': no schurvectors
!                    'I': schuvectors, initializing V and W to the identity
!                    'V': schuvectors, assume V and W already initialized
!
!  RULE            CHARACTER(*)
!                    'periodic' : loop position vector periodically
!                    'upperhess': continue position vector as upper-hessenberg
!                    'lowerhess': continue position vector as lower-hessenberg
!                    'cmw'      : continue position vector as cmv
!                    'random'   : continue position vector randomly
!                    'user'     : continue position vector with user defined function
!
!  FUN             INTEGER PURE FUNCTION FUN(N,P)
!                    takes integer N and binary integer array P of 
!                    dimension N-1 and outputs a binary integer 
!                    if RULE = 'user' FUN is called on every iteration
!                    if RULE /= 'user' FUN is not used
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               INTEGER array of dimension (N-1)
!                    binary array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*N)
!                    array of generators for first sequence of rotations
!
!  D               REAL(8) array of dimension (2*(N+2))
!                    array of generators for complex diagonal matrix
!
!  R               REAL(8) array of dimension (3*N,2)
!                    array of generators for upper-triangular part
!                    of the upper-hessenberg matrix
!
!  S               REAL(8) array of dimension (3*N,2)
!                    array of generators for upper-triangular part
!                    of the pencil
!
! OUTPUT VARIABLES:
!
!  V              COMPLEX(8) array of dimension (N,N)
!                   right schur vectors
!                   if COMPZ = 'N' unused
!                   if COMPZ = 'I' stores right schurvectors in V 
!                   if COMPZ = 'V' update V to store right schurvectors 
!
!  W              COMPLEX(8) array of dimension (N,N)
!                   left schur vectors
!                   if COMPZ = 'N' unused
!                   if COMPZ = 'I' stores left schurvectors in W 
!                   if COMPZ = 'V' update W to store left schurvectors 
!
!  ITS            INTEGER array of dimension (N-1)
!                   Contains the number of iterations per deflation
!
!  INFO           INTEGER
!                   INFO = 0 implies successful computation
!                   INFO = 1 implies no convergence 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uniplusrankone_twistedqz(COMPZ,RULE,FUN,N,P,Q,D,R,S,V,W,ITS,INFO)

end subroutine z_uniplusrankone_twistedqz
