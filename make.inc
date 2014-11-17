############################################################
#
# EISCOR (EIgenSolvers based on unitary CORe transformations)
#
############################################################
LIBNAME := eiscor
MAJOR := 0
MINOR := 1
PATCH := 0
VERSION := $(MAJOR).$(MINOR).$(PATCH)

# directory where the library will be installed
INSTALLDIR := $(HOME)

# fortran compiler and flags
FC := /opt/icc/bin/ifort 
FFLAGS := -std90 -heap-arrays

# libs required by compiler
LIBS := -lifport






