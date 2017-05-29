# Core transformations #

The underlying algorithms in __eiscor__ all rely on efficiently
manipulating essentially 2x2 matrices known as 
_core transformations_. A _core transformation_ Cₖ is any square matrix that is equal to the identity everywhere 
except the 2x2 submatrix at the intersection of rows and columns k and k+1.

# Rotations #

The implementation of __eiscor__ uses for efficiency reasons _rotations_. Rotations are a special case of unitary core transformations. Rotations have the form

```
       ⎡      ⎤
       ⎢ c -s ⎥
  Cₖ = ⎢    _ ⎥
       ⎢ s  c ⎥ ,
       ⎣      ⎦
```
where c is complex and s real. We store the rotations as vectors [ real(c) imag(c) s ].
