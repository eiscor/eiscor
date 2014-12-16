# Expert routines #
This is a complete list of our currently supported _Expert_ routines sorted according to the types of problems they are associated with and their precision. Currently the supported precisions in __eiscor__ are double precision (__d__) and complex double precision (__z__). Routines names follow the scheme __prec_object_function__ 
- __prec__ corresponds to the precision 
- __object__ corresponds to problem type
- __function__ corresponds to solver type

For example, __d_orthhess_qr__ is a double precision (__d__) subroutine for solving orthogonal upper-Hessenberg (__orthhess__) eigenvalue problems using a QR (__qr__) algorithm. 

## Unitary eigensolvers ##
Eigensolvers for unitary matrices can be interacted with at different levels in __eiscor__. At the highest level there are the routines:
- [__d_orthhess_qr__](https://github.com/jaurentz/eiscor/blob/master/src/docs/d_orthhess_qr.md)
- [__z_unihess_qr__](https://github.com/jaurentz/eiscor/blob/master/src/docs/z_unihess_qr.md)

These routines except unitary upper-Hessenberg matrices as input and solve for the eigenvalues and optionally the eigenvectors. The next level of interaction involves inputting the unitary upper-Hessenberg matrix in factored form:
- [__d_orthfact_qr__](https://github.com/jaurentz/eiscor/blob/master/src/docs/d_orthfact_qr.md)
- [__z_unifact_qr__](https://github.com/jaurentz/eiscor/blob/master/src/docs/z_unifact_qr.md)

These routines need less memory but require the user to have a good understanding of the underlying factorizations.
