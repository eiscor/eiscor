# eiscor guide #
Jared L. Aurentz, November 2014

## Installation ##
__eiscor__ is entirely self contained which makes installation incredibly simple.

### Linux ###
To install on a Linux machine simply move into the __eiscor__ root directory, edit the file __make.inc__ to 
suit your system and type:
```
make install
```
This creates a shared object library __libeiscor.so._version___ and copies it into the user specified 
installation directory. 
The installation does not create any symbolic links or export any library paths. Example __make.inc__ files have been included for several common fortran compilers.

### Windows ###
Not explicitly supported at this time. Stay tuned!

### Mac ###
Not explicitly supported at this time. Stay tuned!

## Tests and Examples ##
Once __eiscor__ is installed building the tests and examples is straight forward.

### Linux ###
To run tests on a Linux machine after the library has been installed simply move into the __eiscor__ root directory and type:
```
make tests
```
This will automatically build and run all tests. _Note: This requires the same_ __make.inc__ _file that was used to build the library._ 

To build the examples on a Linux machine type:
```
make examples
```
This will automatically build and run all the examples. _Note: This requires the same_ __make.inc__ _file that was used to build the library._ 

## Removing eiscor ##
If you decide to uninstall __eiscor__

## Reporting issues ##
If you encounter any issues while using __eiscor__ please file an issue on the [eiscor issues](https://github.com/jaurentz/eiscor/issues) page of Github.
