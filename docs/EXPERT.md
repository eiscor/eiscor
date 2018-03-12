# Expert routines #
This is a complete list of our currently supported _Expert_ routines sorted 
according to the types of problems they are associated with and their 
precision. Currently the supported precisions in __eiscor__ are double 
precision (__d__) and complex double precision (__z__). Routines names 
follow the scheme __prec_object_function__ 
- __prec__ corresponds to the precision (i.e. double or complex double)
- __object__ corresponds to the type of data it can act on 
(i.e. scalar or matrix)
- __function__ corresponds to the action performed 
(i.e. sort entries or compute eigenvalues)

For example, __d_orthhess_qr__ is a double precision (__d__) subroutine for 
solving orthogonal upper-Hessenberg (__orthhess__) eigenvalue problems using 
the QR (__qr__) algorithm. 

## Unitary eigensolvers ##
Eigensolvers for unitary matrices can be interacted with at different levels 
in __eiscor__. At the highest level there are the routines:
- [__d_orthhess_qr__](https://github.com/eiscor/eiscor/blob/master/src/double/d_orthhess_qr.f90)
- [__z_unihess_qr__](https://github.com/eiscor/eiscor/blob/master/src/complex_double/z_unihess_qr.f90)

These routines accept unitary upper-Hessenberg matrices as input and 
compute the (possibly real) Schur factorization. The next level of 
interaction involves inputting the unitary upper-Hessenberg matrix in 
factored form:
- [__d_orthfact_qr__](https://github.com/eiscor/eiscor/blob/master/src/double/d_orthfact_qr.f90)
- [__z_unifact_qr__](https://github.com/eiscor/eiscor/blob/master/src/complex_double/z_unifact_qr.f90)

These routines need less memory but require the user to have a good 
understanding of the underlying factorizations. For real matrices the
real Schur factorization is computed. In order to recover the eigenvalues
one should convert to the complex Schur factorization using the following
routines:
- [__d_orthhess_real2complex__](https://github.com/eiscor/eiscor/blob/master/src/double/d_orthhess_real2complex.f90)
- [__d_orthfact_real2complex__](https://github.com/eiscor/eiscor/blob/master/src/double/d_orthfact_real2complex.f90)