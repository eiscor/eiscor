# DOFFQR - Double Orthogonal hessenberg Factored Fast QR eigensolver #

This routine computes the eigenvalues and optionally eigenvectors of a 
real orthogonal upper-Hessenberg matrix that is stored as the product of N-1 Givens' rotations and a complex diagonal matrix with entries +/- 1 using a fast QR algorithm. On output the eigenvalues are stored in the diagonal matrix.

## DOFFQR(COMPZ,N,Q,D,Z,ITS,INFO) ##

### INPUT VARIABLES: ###

__COMPZ__ - CHARACTER
 - 'N': do not compute eigenvectors
 - 'I': stores eigenvectors, initializing Z to the identity
 - 'V': stores eigenvectors, assume Z already initialized

__N__ - INTEGER
 - dimension of matrix

__Q__ - REAL(8) array of dimension (2*(N-1))
 -  array of generators for givens rotations

__D__ - REAL(8) array of dimension (2*N)
 - array of generators for complex diagonal matrix
 - on output contains the eigenvalues

### OUTPUT VARIABLES: ###

__Z__ - REAL(8) array of dimension (N,N)
 - if COMPZ = 'N' unused
 - if COMPZ = 'I' stores eigenvectors, initializing to identity 
 - if COMPZ = 'V' stores eigenvectors, assume Z already initialized

__ITS__ - INTEGER array of dimension (N-1)
 - Contains the number of iterations per deflation

__INFO__ - INTEGER
 - INFO = 1 implies failure to converge
 - INFO = 0 implies successful computation
 - INFO = -1 implies COMPZ is invalid
 - INFO = -2 -2 implies N, Q, or D is invalid
 - INFO = -3 implies Z is invalid

## Example call ##
```fortran
Q(1) = 0d0
Q(2) = 1d0

D(1) = 1d0
D(2) = 0d0
D(3) = -1d0
D(4) = 0d0

call DOFFQR('N',2,Q,D,Z,ITS,INFO)
```
