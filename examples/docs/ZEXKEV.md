# ZEXKEV #
_Zomplex EXample Known EigenValues_

In this example N random eigenvalues on the unit circle are chosen. Based on this eigenvalues a (unitary) inverse eigenvalue problem is solved providing a unitary matrix H with the prescribed eigenvalues. This example solves the unitary eigenvalue problem in two ways:
 1. Form H explicitly and then call [__ZUHFQR__]().
 2. Form H implicitly in factored form and then call [__ZUFFQR__]().

Based on the known eigenvalues the distance to the closest exact eigenvalue is computed.
