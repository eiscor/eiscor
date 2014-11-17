# ZUHFQR - Zomplex Unitary Hessenberg Fast QR eigensolver #

This routine computes the eigenvalues and optionally eigenvectors of a 
complex unitary upper-Hessenberg matrix using a fast QR algorithm. The output is the complex Schur factorization where the triangular part is stored in the same array as the input matrix.

## ZUHFQR(COMPZ,N,H,Z,ITS,WORK,INFO) ##

### INPUT VARIABLES: ###

__COMPZ__ - CHARACTER
 - 'N': do not compute eigenvectors
 - 'I': stores eigenvectors, initializing Z to the identity
 - 'V': stores eigenvectors, assume Z already initialized

__N__ - INTEGER
 - dimension of matrix

__H__ - COMPLEX(8) array of dimension (N,N)
 - unitary hessenberg matrix, assumed that H(ii,jj) = 0 for |ii-jj| > 0
 - on exit contains a diagonal matrix  

__WORK__ - REAL(8) array of dimension (5*N)
 - work space for eigensolver

### OUTPUT VARIABLES: ###

__Z__ - COMPLEX(8) array of dimension (N,N)
 - if COMPZ = 'N' unused
 - if COMPZ = 'I' stores eigenvectors, initializing to identity 
 - if COMPZ = 'V' stores eigenvectors, assume Z already initialized

__ITS__ - INTEGER array of dimension (N-1)
 - Contains the number of iterations per deflation

__INFO__ - INTEGER
 - INFO = 2 implies ZUFFQR failed
 - INFO = 1 implies ZUHRFF failed
 - INFO = 0 implies successful computation
 - INFO = -1 implies COMPZ is invalid
 - INFO = -2 implies N is invalid
 - INFO = -3 implies H is invalid
 - INFO = -4 implies Z is invalid

## Example call ##
```fortran
H(1,1) = cmplx(0d0,0d0,kind=8)
H(2,1) = cmplx(1d0,0d0,kind=8)
H(1,2) = cmplx(1d0,0d0,kind=8)
H(2,2) = cmplx(0d0,0d0,kind=8)

call ZUHFQR('N',2,H,Z,ITS,WORK,INFO)
```
