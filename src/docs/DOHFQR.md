# DOHFQR - Double Orthogonal Hessenberg Fast QR eigensolver #

This routine computes the eigenvalues and optionally eigenvectors of a 
real orthogonal upper-Hessenberg matrix using a fast QR algorithm. The output is the Winter-Murnaghan or real Schur factorization where the block triangular part is stored in the same array as the input matrix.

## DOHFQR(COMPZ,N,H,Z,ITS,WORK,INFO) ##

### INPUT VARIABLES: ###

__COMPZ__ - CHARACTER
 - 'N': do not compute eigenvectors
 - 'I': stores eigenvectors, initializing Z to the identity
 - 'V': stores eigenvectors, assume Z already initialized

__N__ - INTEGER
 - dimension of matrix

__H__ - REAL(8) array of dimension (N,N)
 -  orthogonal hessenberg matrix, assumed that H(ii,jj) = 0 for |ii-jj| > 0
 - on exit contains a block diagonal matrix with blocks no greater than 2x2. 

__WORK__ - REAL(8) array of dimension (4*N)
 - work space for eigensolver

### OUTPUT VARIABLES: ###

__Z__ - REAL(8) array of dimension (N,N)
 - if COMPZ = 'N' unused
 - if COMPZ = 'I' stores eigenvectors, initializing to identity 
 - if COMPZ = 'V' stores eigenvectors, assume Z already initialized

__ITS__ - INTEGER array of dimension (N-1)
 - Contains the number of iterations per deflation

__INFO__ - INTEGER
 - INFO = 2 implies DOFFQR failed
 - INFO = 1 implies DOHRFF failed
 - INFO = 0 implies successful computation
 - INFO = -1 implies COMPZ is invalid
 - INFO = -2 implies N is invalid
 - INFO = -3 implies H is invalid
 - INFO = -4 implies Z is invalid

## Example call ##
```fortran
H(1,1) = 0d0
H(2,1) = 1d0
H(1,2) = 1d0
H(2,2) = 0d0

call DOHFQR('N',2,H,Z,ITS,WORK,INFO)
```
