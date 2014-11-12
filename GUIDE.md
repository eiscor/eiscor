# EISCOR guide #
Jared L. Aurentz, November 2014.

## Installation ##
__EISCOR__ is entirely self contained which makes installation incredibly simple.

### Linux ###
To install on a Linux machine simply move into the __EISCOR__ root directory, edit the file __make.inc__ to 
suit your system and type:
```
make install
```
This creates a shared object library __libeiscor.so._version___ and copies it into the user specified 
installation directory. 
The installation does not create any symbolic links or export library paths. Example __make.inc__ files have been 
included for several common fortran compilers.

### Windows ###
Not supported at this time. Stay tuned!

### Mac ###
Not supported at this time. Stay tuned!

## Tests and Examples ##
Once __EISCOR__ is installed building the tests and examples is straight forward.

### Linux ###
To install on a Linux machine simply move into the __EISCOR__ root directory, edit the file __make.inc__ to 
suit your system and type:
```
make install
```
This creates a shared object library __libeiscor.so._version___ and copies it into the user specified 
installation directory. 
The installation does not create any symbolic links or export library paths. Example __make.inc__ files have been 
included for several common fortran compilers.
