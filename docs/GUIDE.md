# eiscor guide #
Jared L. Aurentz and Thomas Mach, May 2015

## Introduction ##
__eiscor__ is a collection of Fortran 90 subroutines for accurately and 
efficiently solving matrix eigenvalue problems using unitary core 
transformations. These algorithms are the result of an ongoing series of 
papers written by Jared L. Aurentz, Thomas Mach, Raf Vandebril and 
David S. Watkins. (See the [__README.md__](https://github.com/eiscor/eiscor/blob/master/README.md) for a list of related publications.) 

### Current features ###
__eiscor-v0.2.0__ has the following features:
 - double precision eigensolvers for real orthogonal upper-Hessenberg matrices
 - complex double precision eigensolvers for unitary upper-Hessenberg matrices

### Changes from previous version ###
 - changed to a more intuitive naming convention
 - added support for the computation of Gauss-Szegő quadrature rules

## Expert routines ##
Every subroutine in __eiscor__ contains a comment block that describes what 
the subroutine does and the basic interface. Most users will only need to 
interact with a handful of subroutines which we refer to as _Expert_. For 
more detailed descriptions of these _Expert_ routines please see 
[__EXPERT.md__](https://github.com/eiscor/eiscor/blob/master/docs/EXPERT.md).  

## Installation ##
__eiscor__ is entirely self contained which makes installation incredibly 
simple.

### Linux ###
To install on a Linux machine simply move into the __eiscor__ root directory, 
edit the file __make.inc__ to suit your system and type:
```
make install
```
This creates a shared object library __libeiscor.so._version___ and copies 
it into the user specified installation directory. The installation does not 
create any symbolic links or export any library paths. Example __make.inc__ 
files have been included for several common fortran compilers.

### Mac ###
To install on a Mac simply move into the __eiscor__ root directory, edit 
the file __make.inc__ to suit your system and type:
```
make install
```
This creates a shared object library __libeiscor.dylib._version___ and 
copies it into the user specified installation directory. 
The installation does not create any symbolic links or export any library 
paths. Example __make.inc__ files have been included for several common 
fortran compilers.

## Tests and Examples ##
Once __eiscor__ is installed building the tests and examples is straight 
forward.

### Linux ###
To run tests on a Linux machine after the library has been installed simply 
move into the __eiscor__ root directory and type:
```
make tests
```
This will automatically build and run all tests. (_Note: This requires the 
same_ __make.inc__ _file that was used to build the library._) 

To build the examples on a Linux machine type:
```
make examples
```
This will automatically build and run all the examples. (_Note: This 
requires the same_ __make.inc__ _file that was used to build the library._) 

### Mac ###
To run tests on a Mac machine after the library has been installed simply 
move into the __eiscor__ root directory and type:
```
make tests
```
This will automatically build and run all tests. (_Note: This requires the 
same_ __make.inc__ _file that was used to build the library._) 

To build the examples on a Mac machine type:
```
make examples
```
This will automatically build and run all the examples. (_Note: This requires 
the same_ __make.inc__ _file that was used to build the library._) 

## Removing eiscor ##
If the source directory has not been removed simply move into the __eiscor__ 
root directory and type:
```
make uninstall
```
If the source directory has been removed the install directory will have to 
be removed explicitly by the user.

## Questions and issues ##
If you have any questions or encounter any issues while using __eiscor__ 
please file an issue on the [__eiscor__ issues](https://github.com/eiscor/eiscor/issues) page of Github.
