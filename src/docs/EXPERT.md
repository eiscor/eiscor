# Expert routines #
This is a complete list of our currently supported _Expert_ routines sorted according to the types of problems they are associated with and their precision. Currently the supported precisions in __eiscor__ are double precision (__D__) and complex double precision (__Z__). Routines follow an [__LAPACK__](http://www.netlib.org/lapack/lug/node24.html)-style naming scheme where every subroutine consists of three blocks of letters __XYYZZZ__. 
- __X__ corresponds to the precision
- __YY__ corresponds to problem type
- __ZZZ__ corresponds to solver type

For example, __DOHFQR__ is a double precision (__D__) subroutine for solving orthogonal upper-Hessenberg (__OH__) eigenvalue problems using a fast QR (__FQR__) algorithm. 

## Unitary eigensolvers ##
Eigensolvers for unitary matrices can be interacted with at different levels in __eiscor__. At the highest level there are the routines:
- [__DOHFQR__](https://github.com/jaurentz/eiscor/blob/master/src/docs/DOHFQR.md)
- [__ZUHFQR__](https://github.com/jaurentz/eiscor/blob/master/src/docs/ZUHFQR.md)

These routines except unitary upper-Hessenberg matrices as input and solve for the eigenvalues and optionally the eigenvectors. The next level of interaction involves inputting the unitary upper-Hessenberg matrix in factored form:
- [__DOFFQR__](https://github.com/jaurentz/eiscor/blob/master/src/docs/DOFFQR.md)
- [__ZUFFQR__](https://github.com/jaurentz/eiscor/blob/master/src/docs/ZUFFQR.md)

These routines need less memory but require the user to have a good understanding of the underlying factorizations.
