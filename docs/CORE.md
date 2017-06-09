# Core transformations #

The underlying algorithms in __eiscor__ all rely on efficiently
manipulating essentially 2x2 matrices known as 
_core transformations_. A _core transformation_ 
![C_i]
(http://www.sciweavers.org/tex2img.php?eq=C_i&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)
is any square matrix that is equal to the identity everywhere 
except the 2x2 submatrix at the intersection of rows and columns
![i](http://www.sciweavers.org/tex2img.php?eq=i&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)
and
![i+1](http://www.sciweavers.org/tex2img.php?eq=i%2B1&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0).
