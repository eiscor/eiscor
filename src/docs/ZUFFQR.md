# ZUFFQR - Zomplex Unitary hessenberg Factored Fast QR eigensolver #

This routine computes the eigenvalues and optionally eigenvectors of a 
complex unitary upper-Hessenberg matrix that is stored as the product of N-1 Givens' rotations 
and a complex diagonal matrix with entries +/- 1 using a fast QR algorithm. On output the eigenvalues 
are stored in the diagonal matrix.

Each Givens' rotation is stored as 3 real numbers _cr_, _ci_ and _s_ which correspond to the real part of a complex _cosine_,
the imaginary part of a complex _cosine_ and a strictly real _sine_ respectively. These numbers are stored sequentially in __Q__:
```fortran
Q(1) = cr1
Q(2) = ci1
Q(3) = s1
Q(4) = ...
```
The complex diagonal matrix has its real and imaginary entries stored sequentially in __D__. On input __D__ must have unimodular entries:
```fortran
D(1) = dr1
D(2) = di1
D(3) = dr2
D(4) = di2
D(5) = ...
```

## ZUFFQR(COMPZ,N,Q,D,Z,ITS,INFO) ##

### INPUT VARIABLES: ###

__COMPZ__ - CHARACTER
 - 'N': do not compute eigenvectors
 - 'I': stores eigenvectors, initializing Z to the identity
 - 'V': stores eigenvectors, assume Z already initialized

__N__ - INTEGER
 - dimension of matrix

__Q__ - REAL(8) array of dimension (3*(N-1))
 -  array of generators for givens rotations

__D__ - REAL(8) array of dimension (2*N)
 - array of generators for complex diagonal matrix
 - on output contains the eigenvalues

### OUTPUT VARIABLES: ###

__Z__ - COMPLEX(8) array of dimension (N,N)
 - if COMPZ = 'N' unused
 - if COMPZ = 'I' stores eigenvectors, initializing to identity 
 - if COMPZ = 'V' stores eigenvectors, assume Z already initialized

__ITS__ - INTEGER array of dimension (N-1)
 - Contains the number of iterations per deflation

__INFO__ - INTEGER
 - INFO = 1 implies failure to converge
 - INFO = 0 implies successful computation
 - INFO = -1 implies COMPZ is invalid
 - INFO = -2 implies N, Q, or D is invalid
 - INFO = -3 implies Z is invalid

## Example call ##
```fortran
Q(1) = 0d0
Q(2) = 0d0
Q(3) = 1d0

D(1) = 0d0
D(2) = 1d0
D(3) = 1d0
D(4) = 0d0

call ZUFFQR('N',2,Q,D,Z,ITS,INFO)
```
